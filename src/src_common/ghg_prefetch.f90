!***************************************************************************
! ghg_prefetch.f90
! ----------------
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
! \brief       Decompress the next GHG archive while the current one is
!              being processed.
!
!              A LI-COR archive costs about 170 ms more than the same data as
!              a plain file - measured over the base_ghg fixture against
!              base_tlag_opt on identical data - and almost none of that is
!              work this program does. UnZipArchive shells out five times per
!              file, once of them to run 7-Zip over the whole archive, and
!              ReadLicorGhgArchive a sixth to clean up afterwards. That is
!              time spent waiting on another process, so it can be spent while
!              this one computes a flux instead.
!
!              The arrangement is deliberately the simplest that overlaps
!              anything. One archive is prefetched at a time, into one
!              directory, and the request is made only AFTER the current
!              archive's files have been read and deleted - so nothing is ever
!              extracted into a directory something is still reading.
!
!              Nothing here is allowed to change what is processed. A prefetch
!              is claimed only when the archive path matches exactly and the
!              extraction has signalled that it finished; anything else falls
!              through to extracting synchronously, exactly as before. So the
!              worst a failed, slow or stale prefetch can do is waste itself.
!
!              A run that asks for one process - `-j 1` - gets no prefetch,
!              since that is the switch for "do not start anything alongside".
!
! \author      Jonathan Muller
! \note        The sentinel is written by the launcher script as a separate
!              command AFTER the extraction, so its existence means the
!              extraction returned. Testing for the extracted files instead
!              would race: they appear while 7-Zip is still writing them.
! \sa          unzip_archive.f90, read_licor_ghg_archive.f90
!***************************************************************************
module m_ghg_prefetch
    use m_common_global_var
    implicit none
    private

    public :: GhgPrefetchDir, GhgPrefetchClaim, GhgPrefetchStart
    public :: GhgPrefetchRelease, GhgPrefetchCleanup

    !> The archive the outstanding prefetch is for, empty when there is none.
    character(PathLen) :: Pending = ''
    !> Set while the caller is reading out of the prefetch directory, so a new
    !> request cannot extract into it underneath them.
    logical :: InUse = .false.

contains

    !***************************************************************************
    !> \brief Where a prefetched archive is extracted to.
    !***************************************************************************
    character(PathLen) function GhgPrefetchDir()
        GhgPrefetchDir = trim(adjustl(TmpDir)) // 'ghg_prefetch' // slash
    end function GhgPrefetchDir

    !***************************************************************************
    !> \brief Take the prefetch of this archive if there is a finished one.
    !>
    !> ready is .false. for every reason there might be - no request, a request
    !> for a different archive, one that has not finished, one whose extraction
    !> failed - and the caller then does what it always did.
    !***************************************************************************
    subroutine GhgPrefetchClaim(ZipFile, ready)
        character(*), intent(in) :: ZipFile
        logical, intent(out) :: ready
        logical :: ex

        ready = .false.
        if (len_trim(Pending) == 0) return
        if (trim(Pending) /= trim(ZipFile)) return

        inquire(file = trim(GhgPrefetchDir()) // 'ready.flag', exist = ex)
        if (.not. ex) return

        ready = .true.
        InUse = .true.
        Pending = ''
    end subroutine GhgPrefetchClaim

    !***************************************************************************
    !> \brief Done reading the prefetch directory; it may be written again.
    !***************************************************************************
    subroutine GhgPrefetchRelease()
        InUse = .false.
    end subroutine GhgPrefetchRelease

    !***************************************************************************
    !> \brief Ask for an archive to be extracted in the background.
    !>
    !> Silently does nothing when there is nothing to fetch, when a fetch is
    !> already outstanding, or when the caller is still reading the directory.
    !***************************************************************************
    subroutine GhgPrefetchStart(ZipFile)
        character(*), intent(in) :: ZipFile
        integer :: u
        integer :: io_status
        character(PathLen) :: dir
        character(PathLen) :: script
        character(2048) :: comm

        if (len_trim(ZipFile) == 0) return
        if (InUse) return
        if (len_trim(Pending) /= 0) return
        !> -j 1 is the switch for "start nothing alongside this".
        if (NumJobs == 1) return

        dir = GhgPrefetchDir()
        if (OS == 'win') then
            script = trim(adjustl(TmpDir)) // 'ghg_prefetch.bat'
        else
            script = trim(adjustl(TmpDir)) // 'ghg_prefetch.sh'
        end if

        open(newunit = u, file = trim(script), status = 'replace', &
            iostat = io_status)
        if (io_status /= 0) return

        if (OS == 'win') then
            write(u, '(a)') '@echo off'
            write(u, '(a)') 'if exist "' // trim(dir) // '" rmdir /s /q "' &
                // trim(dir) // '"'
            write(u, '(a)') 'mkdir "' // trim(dir) // '"'
            write(u, '(a)') trim(comm_7zip) // ' ' // trim(comm_7zip_x_opt) &
                // ' "' // trim(ZipFile) // '" -o"' // trim(dir) // '"' &
                // comm_out_redirect // comm_err_redirect
            !> Separate command, so it runs only once the extraction returned.
            write(u, '(a)') 'echo ready > "' // trim(dir) // 'ready.flag"'
        else
            write(u, '(a)') '#!/bin/sh'
            write(u, '(a)') 'rm -rf "' // trim(dir) // '"'
            write(u, '(a)') 'mkdir -p "' // trim(dir) // '"'
            write(u, '(a)') trim(comm_7zip) // ' ' // trim(comm_7zip_x_opt) &
                // ' "' // trim(ZipFile) // '" -o"' // trim(dir) // '"' &
                // comm_out_redirect // comm_err_redirect
            write(u, '(a)') 'echo ready > "' // trim(dir) // 'ready.flag"'
        end if
        close(u)

        if (OS == 'win') then
            comm = 'start "" /B cmd /c "' // trim(script) // '"'
        else
            comm = 'sh "' // trim(script) // '" &'
        end if
        call system(trim(comm))

        Pending = ZipFile
    end subroutine GhgPrefetchStart

    !***************************************************************************
    !> \brief Remove the prefetch directory at the end of a run.
    !***************************************************************************
    subroutine GhgPrefetchCleanup()
        character(2048) :: comm

        Pending = ''
        InUse = .false.
        comm = trim(comm_rmdir) // ' "' // trim(GhgPrefetchDir()) // '"' &
            // comm_err_redirect
        call system(trim(comm))
    end subroutine GhgPrefetchCleanup

end module m_ghg_prefetch
