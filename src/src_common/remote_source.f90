!***************************************************************************
! remote_source.f90
! -----------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow®.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version. You should have received a copy
! of the GNU General Public License along with EddyFlow (R). If not,
! see <http://www.gnu.org/licenses/>.
!
! EddyFlow® contains additional Open Source Components. The licenses
! and/or notices these Components can be found in the file LIBRARIES.txt.
!
! EddyFlow® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
!***************************************************************************
!
! \brief       Read input files from a shared Google Drive or Dropbox link.
!
!              Every input - data_path, proj_file, dyn_metadata_file,
!              biom_file, biom_dir, pf_file, to_file, and FCC's ex_file,
!              sa_file, sa_bin_spectra and sa_full_spectra - may be a link to
!              a folder or file shared with "anyone with the link". Nothing
!              needs a login or an API key. The output folder must be local.
!
!              The small inputs are downloaded once, at startup, into TmpDir,
!              and the setting is repointed at the copy. Every reader
!              downstream is unchanged.
!
!              Raw files are not. A season of them does not fit on every disk,
!              so the folder is only LISTED at startup, and each file is given
!              the local path it WILL have, under TmpDir/remote/. The listing
!              is handed to FileListByExt in place of the `dir` it would have
!              run, so name matching, timestamps and ordering are exactly those
!              of a local folder. A file is downloaded when a reader asks for
!              it (RemoteEnsure), the next few are fetched in the background
!              while it is processed, and the least recently used are deleted
!              once more than a handful are on disk.
!
!              Eviction never touches the file being read nor the ones fetched
!              ahead of it, because a period may span two files and because the
!              GHG prefetch extracts the next archive while this one is read.
!
!              How each provider is listed, as found by probing it:
!               - Google Drive: https://drive.google.com/embeddedfolderview?id=
!                 is a plain HTML page, one `flip-entry` per item, linking to
!                 /drive/folders/<id> or /file/d/<id>. A file downloads from
!                 drive.usercontent.google.com/download?id=..&confirm=t, which
!                 skips the virus-scan page Drive puts before large files.
!               - Dropbox: the folder page is a script with no entries in it.
!                 The page loads them with a POST to
!                 /list_shared_link_folder_entries, which needs the CSRF cookie
!                 `t` the page sets. Each entry carries its own share URL, with
!                 its own secure hash; a subfolder is listed with THAT hash and
!                 its path without a leading slash, and a file is that URL with
!                 dl=1. This endpoint is not a published API and may change;
!                 when it does, the listing fails with Fatal error(120) naming
!                 Dropbox, rather than processing a partial folder.
!
!              EDDYFLOW_REMOTE_BASE, when set, replaces every provider host, so
!              the regression harness can serve a fixture from localhost.
!
! \author      Jonathan Muller
! \sa          ghg_prefetch.f90, filelist_by_ext.f90
!***************************************************************************
module m_remote_source
    use m_common_global_var
    use m_log
    implicit none
    private

    public :: IsRemotePath, RemoteRefuseOutput
    public :: RemoteResolveInputs, RemoteFetchFile, RemoteFetchFolder
    public :: RemoteOwnsDir, RemoteWriteFileList
    public :: RemoteEnsure, RemoteIsLocal, RemoteSizeOf
    public :: RemoteAdoptOrder, RemoteCleanup

    type :: RemoteEntry
        !> Full local path the file has, or will have, under the staging dir
        character(PathLen) :: local = ''
        !> Where to download it from
        character(PathLen) :: url = ''
        !> Size as listed, -1 where the provider does not say (Google)
        integer(kind = 8) :: bytes = -1
        !> 0 not on disk, 1 background download started, 2 on disk, 3 failed
        integer :: state = 0
        !> Position in processing order
        integer :: pos = 0
        !> When last asked for, for the least-recently-used eviction
        integer :: used = 0
    end type RemoteEntry

    integer, parameter :: stAbsent = 0
    integer, parameter :: stFetching = 1
    integer, parameter :: stLocal = 2
    integer, parameter :: stFailed = 3

    !> Seconds to wait for a background download before fetching it again
    integer, parameter :: MaxWait = 600

    !> Raw data source
    logical :: Active = .false.
    logical :: Listed = .false.
    logical :: Recurse = .true.
    character(PathLen) :: DataLink = ''
    character(PathLen) :: StagingDir = ''
    type(RemoteEntry), allocatable :: Raw(:)
    integer :: NumRaw = 0
    !> Raw entry indices sorted by local path, for lookup
    integer, allocatable :: ByPath(:)
    !> Raw entry index at each processing position
    integer, allocatable :: ByPos(:)
    integer :: Clock = 0
    integer :: NumFailed = 0
    !> Raw entries currently on disk - a handful, so eviction need not walk
    !> the whole listing for every file read
    integer, allocatable :: OnDisk(:)
    integer :: NumOnDisk = 0

    !> How many files are fetched ahead of the one being read
    integer :: Ahead = 2

    !> Directories already created, so each costs one mkdir
    character(PathLen), allocatable :: MadeDirs(:)
    integer :: NumMadeDirs = 0

    !> Dropbox session: the CSRF token, read once from the cookie jar
    character(256) :: DbxToken = ''
    logical :: CurlChecked = .false.

contains

    !***************************************************************************
    !> \brief Whether a setting holds a link rather than a local path.
    !>
    !> AdjDir and AdjFilePath have usually run over the value already, turning
    !> its slashes into backslashes on Windows, so both spellings are accepted.
    !***************************************************************************
    logical function IsRemotePath(path)
        character(*), intent(in) :: path
        character(8) :: head
        character(PathLen) :: low

        low = Lower(Unmangle(path))
        head = low(1:8)
        IsRemotePath = head(1:7) == 'http://' .or. head(1:8) == 'https://'
    end function IsRemotePath

    !***************************************************************************
    !> \brief Stop if an output location is a link.
    !>
    !> EddyFlow writes to these, and a shared link is read only.
    !***************************************************************************
    subroutine RemoteRefuseOutput(setting, path)
        character(*), intent(in) :: setting
        character(*), intent(in) :: path

        if (.not. IsRemotePath(path)) return
        call LogSayList('  Fatal error(122)> Setting "' // trim(setting) &
            // '" is a shared link:')
        call LogSayList('  Fatal error(122)>   ' // trim(Unmangle(path)))
        call ExceptionHandler(122)
    end subroutine RemoteRefuseOutput

    !***************************************************************************
    !> \brief If the setting links to a file, download it and point at the copy.
    !>
    !> For the files a user may supply rather than have EddyFlow compute:
    !> planar fit, time lag and spectral assessment files, and FCC's ex file.
    !***************************************************************************
    subroutine RemoteFetchFile(setting, path)
        character(*), intent(in) :: setting
        character(*), intent(inout) :: path

        if (.not. IsRemotePath(path)) return
        call FetchSingle(setting, path, setting)
    end subroutine RemoteFetchFile

    !***************************************************************************
    !> \brief If the setting links to a folder, download it whole.
    !>
    !> Files whose name ends with tail (without case; all when tail is empty),
    !> in subfolders too when recursive. For folders that are read as a whole -
    !> biomet, binned and full cospectra - which are small next to raw data.
    !***************************************************************************
    subroutine RemoteFetchFolder(setting, dir, tail, recursive)
        character(*), intent(in) :: setting
        character(*), intent(inout) :: dir
        character(*), intent(in) :: tail
        logical, intent(in) :: recursive

        if (.not. IsRemotePath(dir)) return
        call FetchFolder(setting, dir, tail, recursive)
    end subroutine RemoteFetchFolder

    !***************************************************************************
    !> \brief Download the small inputs and register the raw data link.
    !>
    !> Called once the project file has been read, before anything is opened.
    !> Every argument is repointed at a local copy when it holds a link;
    !> DataPath is repointed at the staging directory the raw files will be
    !> downloaded into, and is listed later, when FileListByExt asks.
    !***************************************************************************
    subroutine RemoteResolveInputs(DataPath, RawRecurse, MetaFile, DynMDFile, &
            BiometFile, BiometDir, BiometTail, BiometRecurse, PrefetchAhead)
        character(*), intent(inout) :: DataPath
        logical, intent(in) :: RawRecurse
        character(*), intent(inout) :: MetaFile
        character(*), intent(inout) :: DynMDFile
        character(*), intent(inout) :: BiometFile
        character(*), intent(inout) :: BiometDir
        character(*), intent(in) :: BiometTail
        logical, intent(in) :: BiometRecurse
        integer, intent(in) :: PrefetchAhead

        if (IsRemotePath(MetaFile)) &
            call FetchSingle('proj_file', MetaFile, 'metadata')
        if (IsRemotePath(DynMDFile)) &
            call FetchSingle('dyn_metadata_file', DynMDFile, 'dynamic_metadata')
        if (IsRemotePath(BiometFile)) &
            call FetchSingle('biom_file', BiometFile, 'biomet')
        if (IsRemotePath(BiometDir)) &
            call RemoteFetchFolder('biom_dir', BiometDir, BiometTail, BiometRecurse)

        if (IsRemotePath(DataPath)) then
            call CheckCurl()
            Active = .true.
            Recurse = RawRecurse
            if (PrefetchAhead >= 0) Ahead = PrefetchAhead
            DataLink = Unmangle(DataPath)
            StagingDir = trim(adjustl(TmpDir)) // 'remote' // slash
            call MakeDir(StagingDir)
            DataPath = StagingDir
            call LogSay(' Raw files are read from a shared link:')
            call LogSay('  ' // trim(DataLink))
            call LogSay('  Each file is downloaded when it is needed and deleted&
                & after use; nothing is kept once the run ends.')
        end if
    end subroutine RemoteResolveInputs

    !***************************************************************************
    !> \brief Whether this directory is the raw data staging directory.
    !***************************************************************************
    logical function RemoteOwnsDir(DirIn)
        character(*), intent(in) :: DirIn

        RemoteOwnsDir = Active
        if (.not. Active) return
        RemoteOwnsDir = trim(adjustl(DirIn)) == trim(StagingDir)
    end function RemoteOwnsDir

    !***************************************************************************
    !> \brief Write the raw listing, in the shape `dir /S /B` would have.
    !>
    !> One full local path per line, for the files whose name ends with Ext
    !> (compared without case, like `dir` and `find -iname`). The folder is
    !> listed on the first call only; FileListByExt and NumberOfFilesInDir
    !> both ask.
    !***************************************************************************
    subroutine RemoteWriteFileList(Ext, OutFile, status)
        character(*), intent(in) :: Ext
        character(*), intent(in) :: OutFile
        integer, intent(out) :: status
        integer :: u
        integer :: i
        integer :: n
        character(PathLen) :: name

        status = 0
        if (.not. Listed) call ListRaw()

        open(newunit = u, file = trim(OutFile), status = 'replace', &
            iostat = status)
        if (status /= 0) return
        n = len_trim(Ext)
        do i = 1, NumRaw
            name = Raw(ByPath(i))%local
            if (len_trim(name) < n) cycle
            if (Lower(name(len_trim(name) - n + 1:len_trim(name))) &
                /= Lower(Ext(1:n))) cycle
            write(u, '(a)') trim(name)
        end do
        close(u)
    end subroutine RemoteWriteFileList

    !***************************************************************************
    !> \brief Make sure a raw file is on disk before it is opened.
    !>
    !> Does nothing for a path that is not from a link, so every reader can
    !> call it unconditionally. A file that cannot be downloaded is reported
    !> once, as Warning(121), and left missing: the reader then fails to open
    !> it and skips the file as it would any unreadable one.
    !>
    !> FetchAhead starts the files that come next in processing order. The
    !> acquisition frequency survey reads files out of order and says no.
    !***************************************************************************
    subroutine RemoteEnsure(path, FetchAhead)
        character(*), intent(in) :: path
        logical, intent(in) :: FetchAhead
        integer :: e
        integer :: p

        if (.not. Active) return
        e = Lookup(path)
        if (e == 0) return

        Clock = Clock + 1
        Raw(e)%used = Clock

        if (Raw(e)%state == stFetching) call AwaitBackground(e)
        if (Raw(e)%state == stAbsent) then
            if (Download(Raw(e)%url, Raw(e)%local)) then
                call MarkLocal(e)
            else
                Raw(e)%state = stFailed
                NumFailed = NumFailed + 1
                call LogSayList(' Warning(121)> Could not download "' &
                    // trim(BaseName(Raw(e)%local)) // '" from the shared link.')
                call ExceptionHandler(121)
            end if
        end if

        if (FetchAhead .and. Ahead > 0) then
            do p = Raw(e)%pos + 1, min(Raw(e)%pos + Ahead, NumRaw)
                call StartBackground(ByPos(p))
            end do
        end if

        call Evict(e)
    end subroutine RemoteEnsure

    !***************************************************************************
    !> \brief Whether a file can be read right now without waiting.
    !>
    !> True for any path that is not from a link. The GHG prefetch uses this
    !> not to run 7-Zip over an archive that is still being downloaded.
    !***************************************************************************
    logical function RemoteIsLocal(path)
        character(*), intent(in) :: path
        integer :: e
        logical :: ex

        RemoteIsLocal = .true.
        if (.not. Active) return
        e = Lookup(path)
        if (e == 0) return
        if (Raw(e)%state == stFetching) then
            inquire(file = trim(Raw(e)%local), exist = ex)
            if (ex) call MarkLocal(e)
        end if
        RemoteIsLocal = Raw(e)%state == stLocal
    end function RemoteIsLocal

    !***************************************************************************
    !> \brief The listed size of a raw file; -1 when unknown or not remote.
    !***************************************************************************
    integer(kind = 8) function RemoteSizeOf(path)
        character(*), intent(in) :: path
        integer :: e

        RemoteSizeOf = -1
        if (.not. Active) return
        e = Lookup(path)
        if (e == 0) return
        RemoteSizeOf = Raw(e)%bytes
    end function RemoteSizeOf

    !***************************************************************************
    !> \brief Take the processing order from the sorted raw file list.
    !>
    !> Until this is called, processing order is listing order.
    !***************************************************************************
    subroutine RemoteAdoptOrder(FileList, nfiles)
        integer, intent(in) :: nfiles
        type(FileListType), intent(in) :: FileList(nfiles)
        integer :: i
        integer :: e
        integer :: p

        if (.not. Active) return
        do i = 1, NumRaw
            Raw(i)%pos = 0
        end do
        p = 0
        do i = 1, nfiles
            e = Lookup(FileList(i)%path)
            if (e == 0) cycle
            if (Raw(e)%pos /= 0) cycle
            p = p + 1
            Raw(e)%pos = p
            ByPos(p) = e
        end do
        !> Files the run will not read go last, never fetched ahead
        do i = 1, NumRaw
            if (Raw(i)%pos /= 0) cycle
            p = p + 1
            Raw(i)%pos = p
            ByPos(p) = i
        end do
    end subroutine RemoteAdoptOrder

    !***************************************************************************
    !> \brief Delete everything downloaded, at the end of a run.
    !***************************************************************************
    subroutine RemoteCleanup()
        character(16) :: count_str

        if (NumFailed > 0) then
            write(count_str, '(i0)') NumFailed
            call LogSayList(' Warning(121)> ' // trim(count_str) &
                // ' raw file(s) could not be downloaded from the shared link;&
                & their periods were skipped. See above for which.')
        end if
        if (len_trim(StagingDir) > 0) &
            call system(trim(comm_rmdir) // ' "' // trim(StagingDir) // '"' &
                // comm_err_redirect)
        call system(trim(comm_rmdir) // ' "' // trim(adjustl(TmpDir)) &
            // 'remote_aux' // slash // '"' // comm_err_redirect)
        Active = .false.
    end subroutine RemoteCleanup


    !===========================================================================
    ! Listing
    !===========================================================================

    !***************************************************************************
    !> \brief List the raw data folder into Raw(:).
    !***************************************************************************
    subroutine ListRaw()
        integer :: i
        character(16) :: count_str

        Listed = .true.
        NumRaw = 0
        allocate(Raw(256))
        call LogSayNoAdv('  Listing the shared folder..')
        call ListFolder(DataLink, StagingDir, Recurse, Raw, NumRaw)
        write(count_str, '(i0)') NumRaw
        call LogSay(' ' // trim(count_str) // ' files.')

        allocate(ByPath(NumRaw), ByPos(NumRaw))
        do i = 1, NumRaw
            ByPath(i) = i
        end do
        if (NumRaw > 1) call SortByPath(ByPath, NumRaw)
        !> Until the run says otherwise, processing order is name order
        do i = 1, NumRaw
            ByPos(i) = ByPath(i)
            Raw(ByPath(i))%pos = i
        end do
    end subroutine ListRaw

    !***************************************************************************
    !> \brief List a shared folder, files only, into List(1:n).
    !>
    !> LocalRoot is where the folder's files will go; subfolders keep their
    !> names below it, so a recursive local listing and this one agree.
    !***************************************************************************
    subroutine ListFolder(link, LocalRoot, recursive, List, n)
        character(*), intent(in) :: link
        character(*), intent(in) :: LocalRoot
        logical, intent(in) :: recursive
        type(RemoteEntry), allocatable, intent(inout) :: List(:)
        integer, intent(inout) :: n
        character(16) :: provider
        character(PathLen) :: id
        character(PathLen) :: key
        character(PathLen) :: hash
        character(PathLen) :: sub
        character(256) :: rlkey

        call CheckCurl()
        call ParseLink(link, provider, id, key, hash, sub, rlkey)
        select case (trim(provider))
            case ('gdrive')
                call ListGdrive(id, LocalRoot, recursive, List, n, 0)
            case ('dropbox')
                call DropboxSession(link)
                call ListDropbox(key, hash, sub, rlkey, LocalRoot, recursive, &
                    List, n, 0)
            case default
                call LinkNotUnderstood(link)
        end select
    end subroutine ListFolder

    !***************************************************************************
    !> \brief One Google Drive folder, and its subfolders if asked.
    !***************************************************************************
    recursive subroutine ListGdrive(id, LocalDir, recursive, List, n, depth)
        character(*), intent(in) :: id
        character(*), intent(in) :: LocalDir
        logical, intent(in) :: recursive
        type(RemoteEntry), allocatable, intent(inout) :: List(:)
        integer, intent(inout) :: n
        integer, intent(in) :: depth
        character(:), allocatable :: page
        character(PathLen) :: page_file
        character(PathLen) :: entry_id
        character(PathLen) :: name
        integer :: at
        integer :: next
        integer :: q
        integer :: stop_at
        logical :: is_folder
        logical :: ok

        if (depth > 32) return
        page_file = trim(adjustl(TmpDir)) // 'remote_page.tmp'
        ok = Curl('-o "' // trim(page_file) // '" "' &
            // trim(Host('https://drive.google.com')) &
            // '/embeddedfolderview?id=' // trim(id) // '"')
        if (ok) call ReadWhole(page_file, page, ok)
        if (.not. ok) call ListingFailed('Google Drive', id)
        if (index(page, 'flip-entries') == 0) call ListingFailed('Google Drive', id)

        at = index(page, 'class="flip-entry" id="entry-')
        do while (at > 0)
            at = at + len('class="flip-entry" id="entry-')
            q = index(page(at:), '"')
            if (q <= 1) exit
            entry_id = page(at:at + q - 2)
            next = index(page(at:), 'class="flip-entry" id="entry-')
            if (next > 0) then
                stop_at = at + next - 2
            else
                stop_at = len(page)
            end if

            is_folder = index(page(at:stop_at), '/drive/folders/') > 0
            q = index(page(at:stop_at), 'flip-entry-title">')
            if (q > 0) then
                name = Between(page(at + q - 1 + len('flip-entry-title">'):stop_at), '<')
                name = DecodeHtml(name)
                if (is_folder) then
                    if (recursive) call ListGdrive(entry_id, &
                        trim(LocalDir) // trim(SafeName(name)) // slash, &
                        recursive, List, n, depth + 1)
                else
                    call Append(List, n, trim(LocalDir) // trim(SafeName(name)), &
                        trim(Host('https://drive.usercontent.google.com')) &
                        // '/download?id=' // trim(entry_id) &
                        // '&export=download&confirm=t', -1_8)
                end if
            end if

            if (next == 0) exit
            at = at + next - 1
        end do
        deallocate(page)
    end subroutine ListGdrive

    !***************************************************************************
    !> \brief Get the Dropbox CSRF cookie from the link's own page, once.
    !***************************************************************************
    subroutine DropboxSession(link)
        character(*), intent(in) :: link
        character(PathLen) :: jar
        character(1024) :: line
        character(256) :: fields(8)
        integer :: u
        integer :: io_status
        integer :: nf
        logical :: ok

        if (len_trim(DbxToken) > 0) return
        jar = trim(adjustl(TmpDir)) // 'remote_cookies.txt'
        ok = Curl('-c "' // trim(jar) // '" -o "' // trim(adjustl(TmpDir)) &
            // 'remote_page.tmp" "' // trim(HostSwap(StripFragment(link))) // '"')
        if (.not. ok) call ListingFailed('Dropbox', link)

        open(newunit = u, file = trim(jar), status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status /= 0) call ListingFailed('Dropbox', link)
        do
            read(u, '(a)', iostat = io_status) line
            if (io_status /= 0) exit
            !> Netscape format: domain, flag, path, secure, expiry, name, value
            call SplitTabs(line, fields, nf)
            if (nf >= 7) then
                if (trim(fields(6)) == 't') DbxToken = fields(7)
            end if
        end do
        close(u)
        if (len_trim(DbxToken) == 0) call ListingFailed('Dropbox', link)
    end subroutine DropboxSession

    !***************************************************************************
    !> \brief One Dropbox folder, and its subfolders if asked.
    !***************************************************************************
    recursive subroutine ListDropbox(key, hash, sub, rlkey, LocalDir, &
            recursive, List, n, depth)
        character(*), intent(in) :: key
        character(*), intent(in) :: hash
        character(*), intent(in) :: sub
        character(*), intent(in) :: rlkey
        character(*), intent(in) :: LocalDir
        logical, intent(in) :: recursive
        type(RemoteEntry), allocatable, intent(inout) :: List(:)
        integer, intent(inout) :: n
        integer, intent(in) :: depth
        character(:), allocatable :: reply
        character(PathLen) :: reply_file
        character(PathLen) :: voucher
        character(PathLen) :: href
        character(PathLen) :: name
        character(PathLen) :: s_key, s_hash, s_sub, s_id
        character(256) :: s_rlkey
        character(16) :: s_provider
        character(2048) :: args
        integer :: obj_start, obj_end, cursor
        integer(kind = 8) :: bytes
        integer :: pages
        logical :: ok
        logical :: is_dir
        logical :: found
        logical :: more

        if (depth > 32) return
        reply_file = trim(adjustl(TmpDir)) // 'remote_page.tmp'
        voucher = ''
        pages = 0
        do
            pages = pages + 1
            args = '-b "' // trim(adjustl(TmpDir)) // 'remote_cookies.txt" -o "' &
                // trim(reply_file) // '"' &
                // ' --data-urlencode "is_xhr=true"' &
                // ' --data-urlencode "t=' // trim(DbxToken) // '"' &
                // ' --data-urlencode "link_key=' // trim(key) // '"' &
                // ' --data-urlencode "link_type=c"' &
                // ' --data-urlencode "secure_hash=' // trim(hash) // '"' &
                // ' --data-urlencode "sub_path=' // trim(sub) // '"' &
                // ' --data-urlencode "rlkey=' // trim(rlkey) // '"'
            if (len_trim(voucher) > 0) args = trim(args) &
                // ' --data-urlencode "voucher=' // trim(voucher) // '"'
            ok = Curl(trim(args) // ' "' // trim(Host('https://www.dropbox.com')) &
                // '/list_shared_link_folder_entries"')
            if (ok) call ReadWhole(reply_file, reply, ok)
            if (.not. ok) call ListingFailed('Dropbox', trim(key) // '/' // trim(sub))

            cursor = index(reply, '"entries"')
            if (cursor == 0) call ListingFailed('Dropbox', trim(key) // '/' // trim(sub))
            cursor = cursor + index(reply(cursor:), '[')
            do
                call NextJsonObject(reply, cursor, obj_start, obj_end)
                if (obj_start == 0) exit
                cursor = obj_end + 1
                name = JsonString(reply(obj_start:obj_end), 'filename', found)
                if (.not. found) cycle
                href = JsonString(reply(obj_start:obj_end), 'href', found)
                if (.not. found) cycle
                is_dir = JsonBool(reply(obj_start:obj_end), 'is_dir')
                if (is_dir) then
                    if (.not. recursive) cycle
                    call ParseLink(href, s_provider, s_id, s_key, s_hash, &
                        s_sub, s_rlkey)
                    if (trim(s_provider) /= 'dropbox') cycle
                    call ListDropbox(s_key, s_hash, s_sub, rlkey, &
                        trim(LocalDir) // trim(SafeName(name)) // slash, &
                        recursive, List, n, depth + 1)
                else
                    bytes = JsonInt(reply(obj_start:obj_end), 'bytes')
                    call Append(List, n, trim(LocalDir) // trim(SafeName(name)), &
                        trim(DropboxDownloadUrl(href)), bytes)
                end if
            end do

            more = JsonBool(reply, 'has_more_entries')
            if (.not. more) exit
            voucher = JsonString(reply, 'next_request_voucher', found)
            if (.not. found .or. len_trim(voucher) == 0 .or. pages > 1000) &
                call ListingFailed('Dropbox', trim(key) // '/' // trim(sub))
        end do
        deallocate(reply)
    end subroutine ListDropbox


    !===========================================================================
    ! Small inputs
    !===========================================================================

    !***************************************************************************
    !> \brief Download one linked file and repoint the setting at the copy.
    !***************************************************************************
    subroutine FetchSingle(setting, path, fallback_name)
        character(*), intent(in) :: setting
        character(*), intent(inout) :: path
        character(*), intent(in) :: fallback_name
        character(PathLen) :: link
        character(PathLen) :: url
        character(PathLen) :: name
        character(PathLen) :: dir

        call CheckCurl()
        link = Unmangle(path)
        call FileLinkTarget(link, url, name)
        if (len_trim(url) == 0) then
            call LogSayList('  Fatal error(120)> "' // trim(setting) &
                // '" must link to a file, not a folder:')
            call LogSayList('  Fatal error(120)>   ' // trim(link))
            call ExceptionHandler(120)
        end if
        if (len_trim(name) == 0) name = fallback_name

        dir = trim(adjustl(TmpDir)) // 'remote_aux' // slash &
            // trim(setting) // slash
        call MakeDir(dir)
        path = trim(dir) // trim(SafeName(name))
        call LogSayNoAdv(' Downloading "' // trim(setting) &
            // '" from the shared link..')
        if (.not. Download(url, path)) then
            call LogSay('')
            call LogSayList('  Fatal error(120)> Could not download "' &
                // trim(setting) // '" from:')
            call LogSayList('  Fatal error(120)>   ' // trim(link))
            call ExceptionHandler(120)
        end if
        call LogSay(' Done.')
    end subroutine FetchSingle

    !***************************************************************************
    !> \brief Download a linked folder whole, and point the setting at the copy.
    !***************************************************************************
    subroutine FetchFolder(setting, dir, tail, recursive)
        character(*), intent(in) :: setting
        character(*), intent(inout) :: dir
        character(*), intent(in) :: tail
        logical, intent(in) :: recursive
        type(RemoteEntry), allocatable :: List(:)
        character(PathLen) :: link
        character(PathLen) :: local_root
        integer :: n
        integer :: i
        integer :: nt
        character(PathLen) :: name

        link = Unmangle(dir)
        local_root = trim(adjustl(TmpDir)) // 'remote_aux' // slash &
            // trim(setting) // slash
        call MakeDir(local_root)
        allocate(List(64))
        n = 0
        call LogSayNoAdv(' Downloading "' // trim(setting) &
            // '" from the shared link..')
        call ListFolder(link, local_root, recursive, List, n)
        nt = len_trim(tail)
        do i = 1, n
            name = List(i)%local
            if (nt > 0) then
                if (len_trim(name) < nt) cycle
                if (Lower(name(len_trim(name) - nt + 1:len_trim(name))) &
                    /= Lower(tail(1:nt))) cycle
            end if
            if (.not. Download(List(i)%url, List(i)%local)) then
                call LogSay('')
                call LogSayList('  Fatal error(120)> Could not download "' &
                    // trim(BaseName(List(i)%local)) // '" of "' // trim(setting) &
                    // '" from:')
                call LogSayList('  Fatal error(120)>   ' // trim(link))
                call ExceptionHandler(120)
            end if
        end do
        call LogSay(' Done.')
        dir = local_root
    end subroutine FetchFolder


    !===========================================================================
    ! Downloads
    !===========================================================================

    !***************************************************************************
    !> \brief Download url to dest, waiting for it. False if it failed.
    !>
    !> Written to dest.part and renamed, so dest never exists half written. A
    !> page of HTML where data was expected - a quota, permission or sign-in
    !> page - counts as a failure, since curl sees HTTP 200 for those.
    !***************************************************************************
    logical function Download(url, dest)
        character(*), intent(in) :: url
        character(*), intent(in) :: dest
        character(PathLen) :: part
        integer :: attempt
        integer :: rename_status

        Download = .false.
        call MakeDir(DirName(dest))
        part = trim(dest) // '.part'
        do attempt = 1, 3
            call DeleteFile(part)
            if (.not. Curl('-o "' // trim(part) // '" "' // trim(url) // '"')) cycle
            if (LooksLikeHtml(part)) cycle
            call DeleteFile(dest)
            rename_status = RenameFile(part, dest)
            if (rename_status == 0) then
                Download = .true.
                return
            end if
        end do
        call DeleteFile(part)
    end function Download

    !***************************************************************************
    !> \brief Start downloading raw entry e without waiting for it.
    !>
    !> The launcher writes dest.failed when curl fails and otherwise renames
    !> dest.part to dest, so dest appearing means the download is complete.
    !***************************************************************************
    subroutine StartBackground(e)
        integer, intent(in) :: e
        character(PathLen) :: script
        character(PathLen) :: part
        character(PathLen) :: failed
        character(16) :: tag
        integer :: u
        integer :: io_status

        if (Raw(e)%state /= stAbsent) return
        call MakeDir(DirName(Raw(e)%local))
        part = trim(Raw(e)%local) // '.part'
        failed = trim(Raw(e)%local) // '.failed'
        call DeleteFile(part)
        call DeleteFile(failed)

        write(tag, '(i0)') e
        if (OS == 'win') then
            script = trim(adjustl(TmpDir)) // 'remote_fetch_' // trim(tag) // '.bat'
        else
            script = trim(adjustl(TmpDir)) // 'remote_fetch_' // trim(tag) // '.sh'
        end if
        open(newunit = u, file = trim(script), status = 'replace', &
            iostat = io_status)
        if (io_status /= 0) return
        if (OS == 'win') then
            write(u, '(a)') '@echo off'
            !> % is special in a batch file, and links are percent-encoded
            write(u, '(a)') trim(CurlLine('-o "' // trim(part) // '" "' &
                // trim(Percents(Raw(e)%url)) // '"'))
            write(u, '(a)') 'if errorlevel 1 (echo failed> "' // trim(failed) &
                // '") else (move /y "' // trim(part) // '" "' &
                // trim(Raw(e)%local) // '" >nul 2>nul)'
            !> The launcher removes itself; one is written per file
            write(u, '(a)') '(goto) 2>nul & del "%~f0"'
        else
            write(u, '(a)') '#!/bin/sh'
            write(u, '(a)') trim(CurlLine('-o "' // trim(part) // '" ''' &
                // trim(Raw(e)%url) // '''')) // ' && mv -f "' // trim(part) &
                // '" "' // trim(Raw(e)%local) // '" || echo failed > "' &
                // trim(failed) // '"'
            write(u, '(a)') 'rm -f "$0"'
        end if
        close(u)

        if (OS == 'win') then
            call system('start "" /B cmd /c "' // trim(script) // '"')
        else
            call system('sh "' // trim(script) // '" &')
        end if
        Raw(e)%state = stFetching
    end subroutine StartBackground

    !***************************************************************************
    !> \brief Wait for a background download to finish, then check it.
    !>
    !> Leaves the entry on disk, or absent so the caller downloads it again.
    !***************************************************************************
    subroutine AwaitBackground(e)
        integer, intent(in) :: e
        integer :: waited
        logical :: done
        logical :: failed

        waited = 0
        do
            inquire(file = trim(Raw(e)%local), exist = done)
            inquire(file = trim(Raw(e)%local) // '.failed', exist = failed)
            if (done .or. failed .or. waited >= MaxWait) exit
            call sleep(1)
            waited = waited + 1
        end do

        call MarkGone(e)
        call DeleteFile(trim(Raw(e)%local) // '.failed')
        if (done) then
            if (LooksLikeHtml(Raw(e)%local)) then
                call DeleteFile(Raw(e)%local)
                call MarkGone(e)
            else
                call MarkLocal(e)
            end if
        end if
    end subroutine AwaitBackground

    !***************************************************************************
    !> \brief Delete the least recently used raw files beyond the budget.
    !>
    !> Never the one just asked for, nor those fetched ahead of it.
    !***************************************************************************
    subroutine Evict(e)
        integer, intent(in) :: e
        integer :: k
        integer :: i
        integer :: victim

        do while (NumOnDisk > Ahead + 3)
            victim = 0
            do k = 1, NumOnDisk
                i = OnDisk(k)
                if (i == e) cycle
                if (Raw(i)%pos > Raw(e)%pos .and. &
                    Raw(i)%pos <= Raw(e)%pos + Ahead) cycle
                if (victim == 0) then
                    victim = i
                else if (Raw(i)%used < Raw(victim)%used) then
                    victim = i
                end if
            end do
            if (victim == 0) exit
            call DeleteFile(Raw(victim)%local)
            call MarkGone(victim)
        end do
    end subroutine Evict

    !> Entry e is on disk
    subroutine MarkLocal(e)
        integer, intent(in) :: e
        integer, allocatable :: grown(:)

        if (Raw(e)%state == stLocal) return
        Raw(e)%state = stLocal
        if (.not. allocated(OnDisk)) allocate(OnDisk(16))
        if (NumOnDisk >= size(OnDisk)) then
            allocate(grown(2 * size(OnDisk)))
            grown(1:NumOnDisk) = OnDisk(1:NumOnDisk)
            call move_alloc(grown, OnDisk)
        end if
        NumOnDisk = NumOnDisk + 1
        OnDisk(NumOnDisk) = e
    end subroutine MarkLocal

    !> Entry e is not on disk (any more)
    subroutine MarkGone(e)
        integer, intent(in) :: e
        integer :: k

        Raw(e)%state = stAbsent
        do k = 1, NumOnDisk
            if (OnDisk(k) /= e) cycle
            OnDisk(k) = OnDisk(NumOnDisk)
            NumOnDisk = NumOnDisk - 1
            exit
        end do
    end subroutine MarkGone

    !***************************************************************************
    !> \brief Run curl with the common options; true if it succeeded.
    !***************************************************************************
    logical function Curl(args)
        character(*), intent(in) :: args

        Curl = system(trim(CurlLine(args))) == 0
    end function Curl

    character(4096) function CurlLine(args)
        character(*), intent(in) :: args
        character(16) :: exe

        exe = 'curl'
        if (OS == 'win') exe = 'curl.exe'
        !> --fail: an HTTP error is an error. -L: shares redirect.
        CurlLine = trim(exe) // ' -sS -L --fail --retry 2 -A "Mozilla/5.0" ' &
            // trim(args) // ' ' // comm_out_redirect // comm_err_redirect
    end function CurlLine

    !***************************************************************************
    !> \brief Stop, once, if there is no curl to download with.
    !***************************************************************************
    subroutine CheckCurl()
        if (CurlChecked) return
        CurlChecked = .true.
        if (system(trim(CurlLine('--version'))) == 0) return
        call LogSayList('  Fatal error(120)> Reading from a shared link needs curl,&
            & which was not found.')
        call LogSayList('  Fatal error(120)> It ships with Windows 10 and later,&
            & macOS and most Linux distributions.')
        call ExceptionHandler(120)
    end subroutine CheckCurl

    !***************************************************************************
    !> \brief True if a downloaded file starts like an HTML page.
    !***************************************************************************
    logical function LooksLikeHtml(path)
        character(*), intent(in) :: path
        character(64) :: head
        character(PathLen) :: low
        integer :: u
        integer :: io_status
        integer :: fsize

        LooksLikeHtml = .false.
        inquire(file = trim(path), size = fsize)
        if (fsize <= 0) then
            !> An empty download is no better than an error page
            LooksLikeHtml = fsize == 0
            return
        end if
        head = ''
        open(newunit = u, file = trim(path), access = 'stream', &
            form = 'unformatted', status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status /= 0) return
        read(u, iostat = io_status) head(1:min(64, fsize))
        close(u)
        low = Lower(adjustl(head))
        head = low(1:64)
        LooksLikeHtml = head(1:14) == '<!doctype html' .or. head(1:5) == '<html'
    end function LooksLikeHtml


    !===========================================================================
    ! Links
    !===========================================================================

    !***************************************************************************
    !> \brief Pick a share link apart.
    !>
    !> Google Drive: id is the folder or file id.
    !> Dropbox (/scl/fo/<key>/<hash>[/<sub path>]?rlkey=..): key, hash, the
    !> decoded sub path with no leading slash, and rlkey.
    !> provider is left empty for a link neither understands.
    !***************************************************************************
    subroutine ParseLink(link_in, provider, id, key, hash, sub, rlkey)
        character(*), intent(in) :: link_in
        character(*), intent(out) :: provider
        character(*), intent(out) :: id
        character(*), intent(out) :: key
        character(*), intent(out) :: hash
        character(*), intent(out) :: sub
        character(*), intent(out) :: rlkey
        character(PathLen) :: link
        character(PathLen) :: path_part
        integer :: at
        integer :: q
        integer :: slash_at

        provider = ''
        id = ''
        key = ''
        hash = ''
        sub = ''
        rlkey = ''
        link = StripFragment(Unmangle(link_in))

        if (index(link, 'google.com') > 0 .or. index(link, '/drive/folders/') > 0 &
            .or. index(link, '/file/d/') > 0 .or. index(link, 'embeddedfolderview') > 0) then
            at = index(link, '/folders/')
            if (at > 0) then
                id = UpTo(link(at + len('/folders/'):), '/?&#')
            else
                at = index(link, '/file/d/')
                if (at > 0) then
                    id = UpTo(link(at + len('/file/d/'):), '/?&#')
                else
                    at = index(link, 'id=')
                    if (at > 0) id = UpTo(link(at + 3:), '&#')
                end if
            end if
            if (len_trim(id) > 0) provider = 'gdrive'
            return
        end if

        at = index(link, '/scl/fo/')
        if (at > 0) then
            path_part = UpTo(link(at + len('/scl/fo/'):), '?#')
            slash_at = index(path_part, '/')
            if (slash_at <= 1) return
            key = path_part(1:slash_at - 1)
            path_part = path_part(slash_at + 1:)
            slash_at = index(path_part, '/')
            if (slash_at == 0) then
                hash = path_part
            else
                hash = path_part(1:slash_at - 1)
                sub = UrlDecode(path_part(slash_at + 1:))
            end if
            q = index(link, 'rlkey=')
            if (q > 0) rlkey = UpTo(link(q + len('rlkey='):), '&#')
            if (len_trim(hash) > 0) provider = 'dropbox'
        end if
    end subroutine ParseLink

    !***************************************************************************
    !> \brief Where a linked FILE downloads from, and its name if the link says.
    !>
    !> url is empty when the link names a folder.
    !***************************************************************************
    subroutine FileLinkTarget(link, url, name)
        character(*), intent(in) :: link
        character(*), intent(out) :: url
        character(*), intent(out) :: name
        character(16) :: provider
        character(PathLen) :: id, key, hash, sub
        character(256) :: rlkey
        character(PathLen) :: path_part
        integer :: at

        url = ''
        name = FragmentName(link)
        if (index(link, '/drive/folders/') > 0) return

        if (index(link, '/scl/fi/') > 0) then
            !> A single shared file: /scl/fi/<id>/<name>?rlkey=..
            at = index(link, '/scl/fi/')
            path_part = UpTo(link(at + len('/scl/fi/'):), '?#')
            if (len_trim(name) == 0) name = UrlDecode(BaseNameOf(path_part, '/'))
            url = DropboxDownloadUrl(StripFragment(link))
            return
        end if

        call ParseLink(link, provider, id, key, hash, sub, rlkey)
        select case (trim(provider))
            case ('gdrive')
                url = trim(Host('https://drive.usercontent.google.com')) &
                    // '/download?id=' // trim(id) // '&export=download&confirm=t'
            case ('dropbox')
                !> A file inside a shared folder has a sub path; the folder's own
                !> link has none, and is a folder
                if (len_trim(sub) == 0) return
                if (len_trim(name) == 0) name = BaseNameOf(sub, '/')
                url = DropboxDownloadUrl(StripFragment(link))
        end select
    end subroutine FileLinkTarget

    !***************************************************************************
    !> \brief A Dropbox share URL made to download: dl=1, st= dropped.
    !***************************************************************************
    character(PathLen) function DropboxDownloadUrl(href)
        character(*), intent(in) :: href
        character(PathLen) :: base
        character(PathLen) :: query
        character(PathLen) :: kept
        character(PathLen) :: item
        integer :: q
        integer :: amp

        q = index(href, '?')
        if (q == 0) then
            DropboxDownloadUrl = trim(HostSwap(href)) // '?dl=1'
            return
        end if
        base = href(1:q - 1)
        query = href(q + 1:)
        kept = ''
        do while (len_trim(query) > 0)
            amp = index(query, '&')
            if (amp == 0) then
                item = query
                query = ''
            else
                item = query(1:amp - 1)
                query = query(amp + 1:)
            end if
            if (item(1:3) == 'dl=' .or. item(1:3) == 'st=') cycle
            if (len_trim(item) == 0) cycle
            kept = trim(kept) // trim(item) // '&'
        end do
        DropboxDownloadUrl = trim(HostSwap(base)) // '?' // trim(kept) // 'dl=1'
    end function DropboxDownloadUrl

    !***************************************************************************
    !> \brief The provider host, or EDDYFLOW_REMOTE_BASE when that is set.
    !***************************************************************************
    character(PathLen) function Host(default)
        character(*), intent(in) :: default
        integer :: n
        integer :: env_status

        call get_environment_variable('EDDYFLOW_REMOTE_BASE', Host, n, env_status)
        if (env_status /= 0 .or. n == 0) then
            Host = default
        else if (Host(n:n) == '/') then
            Host(n:n) = ' '
        end if
    end function Host

    !***************************************************************************
    !> \brief A provider URL with its host replaced by the test base, if set.
    !***************************************************************************
    character(PathLen) function HostSwap(url)
        character(*), intent(in) :: url
        character(PathLen) :: base
        integer :: scheme_end
        integer :: path_start

        HostSwap = url
        base = Host('')
        if (len_trim(base) == 0) return
        scheme_end = index(url, '://')
        if (scheme_end == 0) return
        path_start = index(url(scheme_end + 3:), '/')
        if (path_start == 0) return
        HostSwap = trim(base) // url(scheme_end + 2 + path_start:)
    end function HostSwap

    subroutine LinkNotUnderstood(link)
        character(*), intent(in) :: link

        call LogSayList('  Fatal error(120)> Not a Google Drive or Dropbox&
            & folder link that EddyFlow can read:')
        call LogSayList('  Fatal error(120)>   ' // trim(link))
        call ExceptionHandler(120)
    end subroutine LinkNotUnderstood

    subroutine ListingFailed(provider, what)
        character(*), intent(in) :: provider
        character(*), intent(in) :: what

        call LogSay('')
        call LogSayList('  Fatal error(120)> Could not list the ' // provider &
            // ' folder:')
        call LogSayList('  Fatal error(120)>   ' // trim(what))
        if (provider == 'Dropbox') then
            call LogSayList('  Fatal error(120)> Either the link is not shared&
                & with "anyone with the link", or Dropbox')
            call LogSayList('  Fatal error(120)> changed its listing format&
                & (not recognised).')
        else
            call LogSayList('  Fatal error(120)> Check that it is shared with&
                & "anyone with the link".')
        end if
        call ExceptionHandler(120)
    end subroutine ListingFailed


    !===========================================================================
    ! JSON, just enough for the Dropbox listing
    !===========================================================================

    !***************************************************************************
    !> \brief The next {...} at array level from cursor, skipping strings.
    !>
    !> obj_start is 0 when the array ends first.
    !***************************************************************************
    subroutine NextJsonObject(text, cursor, obj_start, obj_end)
        character(*), intent(in) :: text
        integer, intent(in) :: cursor
        integer, intent(out) :: obj_start
        integer, intent(out) :: obj_end
        integer :: i
        integer :: depth
        logical :: in_string

        obj_start = 0
        obj_end = 0
        i = cursor
        do while (i <= len(text))
            if (text(i:i) == '{') exit
            if (text(i:i) == ']') return
            i = i + 1
        end do
        if (i > len(text)) return
        obj_start = i
        depth = 0
        in_string = .false.
        do while (i <= len(text))
            if (in_string) then
                if (text(i:i) == achar(92)) then
                    i = i + 1
                else if (text(i:i) == '"') then
                    in_string = .false.
                end if
            else
                select case (text(i:i))
                    case ('"')
                        in_string = .true.
                    case ('{', '[')
                        depth = depth + 1
                    case ('}', ']')
                        depth = depth - 1
                        if (depth == 0) then
                            obj_end = i
                            return
                        end if
                end select
            end if
            i = i + 1
        end do
        obj_start = 0
    end subroutine NextJsonObject

    !> Position just after `"key":` and any blanks, 0 if absent
    integer function JsonValueAt(text, key)
        character(*), intent(in) :: text
        character(*), intent(in) :: key
        integer :: at

        JsonValueAt = 0
        at = index(text, '"' // key // '"')
        if (at == 0) return
        at = at + len(key) + 2
        do while (at <= len(text))
            if (text(at:at) /= ' ' .and. text(at:at) /= ':') exit
            at = at + 1
        end do
        if (at <= len(text)) JsonValueAt = at
    end function JsonValueAt

    character(PathLen) function JsonString(text, key, found)
        character(*), intent(in) :: text
        character(*), intent(in) :: key
        logical, intent(out) :: found
        integer :: at
        integer :: i
        integer :: o
        integer :: code
        integer :: io_status

        JsonString = ''
        found = .false.
        at = JsonValueAt(text, key)
        if (at == 0) return
        if (text(at:at) /= '"') return
        found = .true.
        o = 0
        i = at + 1
        do while (i <= len(text) .and. o < PathLen - 4)
            if (text(i:i) == '"') exit
            if (text(i:i) == achar(92) .and. i < len(text)) then
                i = i + 1
                select case (text(i:i))
                    case ('u')
                        read(text(i + 1:i + 4), '(z4)', iostat = io_status) code
                        if (io_status /= 0) code = iachar('?')
                        call PutUtf8(JsonString, o, code)
                        i = i + 4
                    case ('n')
                        o = o + 1
                        JsonString(o:o) = ' '
                    case ('t')
                        o = o + 1
                        JsonString(o:o) = ' '
                    case default
                        o = o + 1
                        JsonString(o:o) = text(i:i)
                end select
            else
                o = o + 1
                JsonString(o:o) = text(i:i)
            end if
            i = i + 1
        end do
    end function JsonString

    logical function JsonBool(text, key)
        character(*), intent(in) :: text
        character(*), intent(in) :: key
        integer :: at

        JsonBool = .false.
        at = JsonValueAt(text, key)
        if (at == 0) return
        if (at + 3 > len(text)) return
        JsonBool = text(at:at + 3) == 'true'
    end function JsonBool

    integer(kind = 8) function JsonInt(text, key)
        character(*), intent(in) :: text
        character(*), intent(in) :: key
        integer :: at
        integer :: e
        integer :: io_status

        JsonInt = -1
        at = JsonValueAt(text, key)
        if (at == 0) return
        e = at
        do while (e <= len(text))
            if (index('0123456789', text(e:e)) == 0) exit
            e = e + 1
        end do
        if (e == at) return
        read(text(at:e - 1), *, iostat = io_status) JsonInt
        if (io_status /= 0) JsonInt = -1
    end function JsonInt

    !> Append a code point as UTF-8
    subroutine PutUtf8(s, o, code)
        character(*), intent(inout) :: s
        integer, intent(inout) :: o
        integer, intent(in) :: code

        if (code < 128) then
            o = o + 1
            s(o:o) = achar(code)
        else if (code < 2048) then
            s(o + 1:o + 1) = achar(192 + code / 64)
            s(o + 2:o + 2) = achar(128 + mod(code, 64))
            o = o + 2
        else
            s(o + 1:o + 1) = achar(224 + code / 4096)
            s(o + 2:o + 2) = achar(128 + mod(code / 64, 64))
            s(o + 3:o + 3) = achar(128 + mod(code, 64))
            o = o + 3
        end if
    end subroutine PutUtf8


    !===========================================================================
    ! Strings and files
    !===========================================================================

    subroutine Append(List, n, local, url, bytes)
        type(RemoteEntry), allocatable, intent(inout) :: List(:)
        integer, intent(inout) :: n
        character(*), intent(in) :: local
        character(*), intent(in) :: url
        integer(kind = 8), intent(in) :: bytes
        type(RemoteEntry), allocatable :: grown(:)

        if (n >= size(List)) then
            allocate(grown(2 * size(List)))
            grown(1:n) = List(1:n)
            call move_alloc(grown, List)
        end if
        n = n + 1
        List(n)%local = local
        List(n)%url = url
        List(n)%bytes = bytes
        List(n)%state = stAbsent
        List(n)%pos = 0
        List(n)%used = 0
    end subroutine Append

    !> Raw entry with this local path, by binary search; 0 if none
    integer function Lookup(path)
        character(*), intent(in) :: path
        integer :: lo, hi, mid

        Lookup = 0
        if (.not. Listed) return
        lo = 1
        hi = NumRaw
        do while (lo <= hi)
            mid = (lo + hi) / 2
            if (Raw(ByPath(mid))%local == path) then
                Lookup = ByPath(mid)
                return
            else if (Raw(ByPath(mid))%local < path) then
                lo = mid + 1
            else
                hi = mid - 1
            end if
        end do
    end function Lookup

    !> Sort entry indices by local path (bottom-up merge sort)
    subroutine SortByPath(idx, n)
        integer, intent(in) :: n
        integer, intent(inout) :: idx(n)
        integer, allocatable :: tmp(:)
        integer :: width, lo, mid, hi, i, j, k

        allocate(tmp(n))
        width = 1
        do while (width < n)
            lo = 1
            do while (lo <= n)
                mid = min(lo + width - 1, n)
                hi = min(lo + 2 * width - 1, n)
                i = lo
                j = mid + 1
                k = lo
                do while (i <= mid .and. j <= hi)
                    if (Raw(idx(j))%local < Raw(idx(i))%local) then
                        tmp(k) = idx(j)
                        j = j + 1
                    else
                        tmp(k) = idx(i)
                        i = i + 1
                    end if
                    k = k + 1
                end do
                do while (i <= mid)
                    tmp(k) = idx(i)
                    i = i + 1
                    k = k + 1
                end do
                do while (j <= hi)
                    tmp(k) = idx(j)
                    j = j + 1
                    k = k + 1
                end do
                lo = lo + 2 * width
            end do
            idx(1:n) = tmp(1:n)
            width = 2 * width
        end do
    end subroutine SortByPath

    !> A whole file into one string
    subroutine ReadWhole(path, text, ok)
        character(*), intent(in) :: path
        character(:), allocatable, intent(out) :: text
        logical, intent(out) :: ok
        integer :: u
        integer :: fsize
        integer :: io_status

        ok = .false.
        inquire(file = trim(path), size = fsize)
        if (fsize <= 0) then
            text = ''
            return
        end if
        allocate(character(fsize) :: text)
        open(newunit = u, file = trim(path), access = 'stream', &
            form = 'unformatted', status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status /= 0) return
        read(u, iostat = io_status) text
        close(u)
        ok = io_status == 0
    end subroutine ReadWhole

    subroutine MakeDir(dir)
        character(*), intent(in) :: dir
        character(PathLen) :: d
        character(PathLen), allocatable :: grown(:)
        integer :: i

        d = dir
        do i = 1, NumMadeDirs
            if (MadeDirs(i) == d) return
        end do
        if (.not. allocated(MadeDirs)) allocate(MadeDirs(16))
        if (NumMadeDirs >= size(MadeDirs)) then
            allocate(grown(2 * size(MadeDirs)))
            grown(1:NumMadeDirs) = MadeDirs(1:NumMadeDirs)
            call move_alloc(grown, MadeDirs)
        end if
        NumMadeDirs = NumMadeDirs + 1
        MadeDirs(NumMadeDirs) = d
        if (OS == 'win') then
            call system('mkdir "' // trim(d) // '"' // comm_err_redirect)
        else
            call system('mkdir -p "' // trim(d) // '"' // comm_err_redirect)
        end if
    end subroutine MakeDir

    subroutine DeleteFile(path)
        character(*), intent(in) :: path
        integer :: u
        integer :: io_status
        logical :: ex

        inquire(file = trim(path), exist = ex)
        if (.not. ex) return
        open(newunit = u, file = trim(path), status = 'old', iostat = io_status)
        if (io_status == 0) close(u, status = 'delete', iostat = io_status)
    end subroutine DeleteFile

    integer function RenameFile(from, to)
        character(*), intent(in) :: from
        character(*), intent(in) :: to

        call rename(trim(from), trim(to), RenameFile)
    end function RenameFile

    !> The setting as it was typed: AdjDir/AdjFilePath may have turned its
    !> slashes into backslashes and added a trailing one. And the interface
    !> saves projects with QSettings, which puts a value in double quotes when
    !> it contains = ; or , - as every Dropbox link does, in ?rlkey= - and the
    !> engine's ini reader keeps those quotes.
    character(PathLen) function Unmangle(path)
        character(*), intent(in) :: path
        integer :: i
        integer :: n

        Unmangle = adjustl(path)
        do i = 1, len_trim(Unmangle)
            if (Unmangle(i:i) == achar(92)) Unmangle(i:i) = '/'
        end do
        call DropLast(Unmangle, '/')
        if (Unmangle(1:1) == '"') then
            Unmangle = Unmangle(2:)
            call DropLast(Unmangle, '"')
            call DropLast(Unmangle, '/')
        end if
    contains
        subroutine DropLast(s, c)
            character(*), intent(inout) :: s
            character(1), intent(in) :: c

            n = len_trim(s)
            if (n == 0) return
            if (s(n:n) == c) s(n:n) = ' '
        end subroutine DropLast
    end function Unmangle

    character(PathLen) function StripFragment(link)
        character(*), intent(in) :: link
        integer :: h

        StripFragment = link
        h = index(link, '#')
        if (h > 0) StripFragment = link(1:h - 1)
    end function StripFragment

    !> Last element of the `path=` in the fragment the GUI adds, if any
    character(PathLen) function FragmentName(link)
        character(*), intent(in) :: link
        character(PathLen) :: frag
        integer :: h
        integer :: p

        FragmentName = ''
        h = index(link, '#')
        if (h == 0) return
        frag = link(h + 1:)
        p = index(frag, 'path=')
        if (p == 0) return
        frag = UpTo(frag(p + 5:), '&')
        FragmentName = UrlDecode(BaseNameOf(frag, '/'))
    end function FragmentName

    !> s up to (not including) the first of any character in stops
    character(PathLen) function UpTo(s, stops)
        character(*), intent(in) :: s
        character(*), intent(in) :: stops
        integer :: e

        e = scan(s, stops)
        if (e == 0) then
            UpTo = s
        else
            UpTo = s(1:e - 1)
        end if
    end function UpTo

    character(PathLen) function Between(s, stop_char)
        character(*), intent(in) :: s
        character(*), intent(in) :: stop_char
        integer :: e

        e = index(s, stop_char)
        if (e == 0) then
            Between = s
        else
            Between = s(1:e - 1)
        end if
    end function Between

    character(PathLen) function BaseNameOf(path, sep)
        character(*), intent(in) :: path
        character(*), intent(in) :: sep

        BaseNameOf = path(index(trim(path), sep, .true.) + 1:)
    end function BaseNameOf

    character(PathLen) function BaseName(path)
        character(*), intent(in) :: path

        BaseName = BaseNameOf(path, slash)
    end function BaseName

    character(PathLen) function DirName(path)
        character(*), intent(in) :: path

        DirName = path(1:index(trim(path), slash, .true.))
    end function DirName

    !> A remote name that is safe as a file name here
    character(PathLen) function SafeName(name)
        character(*), intent(in) :: name
        integer :: i

        SafeName = name
        do i = 1, len_trim(SafeName)
            if (index('<>:"/|?*' // achar(92), SafeName(i:i)) > 0) &
                SafeName(i:i) = '_'
        end do
    end function SafeName

    character(PathLen) function DecodeHtml(s)
        character(*), intent(in) :: s

        DecodeHtml = s
        call Replace(DecodeHtml, '&amp;', '&')
        call Replace(DecodeHtml, '&#39;', "'")
        call Replace(DecodeHtml, '&quot;', '"')
        call Replace(DecodeHtml, '&lt;', '<')
        call Replace(DecodeHtml, '&gt;', '>')
    end function DecodeHtml

    subroutine Replace(s, what, with)
        character(*), intent(inout) :: s
        character(*), intent(in) :: what
        character(*), intent(in) :: with
        integer :: at
        integer :: from

        from = 1
        do
            at = index(s(from:), what)
            if (at == 0) exit
            at = from + at - 1
            s = s(1:at - 1) // with // s(at + len(what):)
            from = at + len(with)
        end do
    end subroutine Replace

    character(PathLen) function UrlDecode(s)
        character(*), intent(in) :: s
        integer :: i
        integer :: o
        integer :: code
        integer :: io_status

        UrlDecode = ''
        o = 0
        i = 1
        do while (i <= len_trim(s))
            if (s(i:i) == '%' .and. i + 2 <= len_trim(s)) then
                read(s(i + 1:i + 2), '(z2)', iostat = io_status) code
                if (io_status == 0) then
                    o = o + 1
                    UrlDecode(o:o) = achar(code)
                    i = i + 3
                    cycle
                end if
            end if
            o = o + 1
            if (s(i:i) == '+') then
                UrlDecode(o:o) = ' '
            else
                UrlDecode(o:o) = s(i:i)
            end if
            i = i + 1
        end do
    end function UrlDecode

    !> Double every %, for a batch file
    character(PathLen) function Percents(s)
        character(*), intent(in) :: s
        integer :: i
        integer :: o

        Percents = ''
        o = 0
        do i = 1, len_trim(s)
            o = o + 1
            Percents(o:o) = s(i:i)
            if (s(i:i) == '%') then
                o = o + 1
                Percents(o:o) = '%'
            end if
        end do
    end function Percents

    subroutine SplitTabs(line, fields, nf)
        character(*), intent(in) :: line
        character(*), intent(out) :: fields(:)
        integer, intent(out) :: nf
        integer :: from
        integer :: t

        fields = ''
        nf = 0
        from = 1
        do while (nf < size(fields))
            t = index(line(from:), achar(9))
            nf = nf + 1
            if (t == 0) then
                fields(nf) = line(from:)
                exit
            end if
            fields(nf) = line(from:from + t - 2)
            from = from + t
        end do
    end subroutine SplitTabs

    character(PathLen) function Lower(s)
        character(*), intent(in) :: s
        integer :: i
        integer :: c

        Lower = s
        do i = 1, len_trim(Lower)
            c = iachar(Lower(i:i))
            if (c >= 65 .and. c <= 90) Lower(i:i) = achar(c + 32)
        end do
    end function Lower

end module m_remote_source
