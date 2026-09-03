!***************************************************************************
! unzip_archive.f90
! -----------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
! Copyright © 2011-2026, LI-COR Biosciences, Gerardo Fratini
! Copyright © 2026-    , ETH Zurich, Jonathan Muller
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
! \brief       Unzip archive, containing at max a status, a metadata and a data file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UnZipArchive(ZipFile, MetaExt, DataExt, MetaFile, DataFile, &
    BiometFile, BiometMetaFile, skip_file, WorkDir, prefetched)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: ZipFile
    character(*), intent(in) :: MetaExt
    character(*), intent(in) :: DataExt
    character(*), intent(out) :: MetaFile
    character(*), intent(out) :: DataFile
    character(*), intent(out) :: BiometFile
    character(*), intent(out) :: BiometMetaFile
    logical, intent(out) :: skip_file
    !> The directory to work in, and whether the archive is already extracted
    !> there. Both are required rather than optional because this is an
    !> external procedure with no explicit interface, and an optional argument
    !> without one is not something the standard defines - gfortran infers the
    !> interface from whichever call it compiles first and rejects the other.
    character(*), intent(in) :: WorkDir
    logical, intent(in) :: prefetched
    !> local variables

    integer :: io_status
    integer :: del_status
    integer :: unzip_status
    integer :: dir_status
    character(CommLen) :: comm
    character(128) :: dataline
    character(128) :: TmpString


    call clearstr(MetaFile)
    call clearstr(DataFile)
    call clearstr(BiometFile)
    call clearstr(BiometMetaFile)

    !> Already extracted by something else: skip the deleting and 7-Zip and
    !> go straight to reading back what is there.
    if (prefetched) go to 100

    !> Delete residual files in tmp folder
    comm = trim(comm_del) // ' "' // trim(adjustl(TmpDir)) &
        // '"*.*' // trim(adjustl(DataExt)) &
        // ' ' // comm_err_redirect
    del_status = system(comm)
    comm = trim(comm_del) // ' "' // trim(adjustl(TmpDir)) &
        // '"*.tmp' // ' ' // comm_err_redirect
    del_status = system(comm)

    !> Extract files from archive
    comm = trim(comm_7zip) // ' ' // trim(comm_7zip_x_opt) &
        // ' "' // ZipFile(1:len_trim(ZipFile)) // '" -o"' &
        // trim(adjustl(TmpDir)) // '"' &
        // comm_out_redirect // comm_err_redirect

    unzip_status = system(comm)
    if (unzip_status /= 0) then
        call ExceptionHandler(14)
        skip_file = .true.
        return
    end if
    call clearstr(comm)

    !> Everything below reads back what is now in WorkDir, and is the same
    !> whether this routine extracted it or something else did.
100 continue

    !> One listing of everything the archive held, rather than one per
    !> extension. Both were only ever asking "what is in here", and each cost a
    !> shell of its own - about 50 ms, against the 100 ms the extraction itself
    !> takes. Classifying the names here instead is free.
    comm = trim(adjustl(comm_dir)) // ' "' &
        // trim(adjustl(WorkDir)) // '"*.*' &
        // ' > "' // trim(adjustl(WorkDir)) // 'arch_flist.tmp" ' &
        // comm_err_redirect
    dir_status = system(comm)

    MetaFile = 'none'
    DataFile = 'none'
    BiometMetaFile = 'none'
    BiometFile = 'none'

    open(udf, file = trim(adjustl(WorkDir)) // 'arch_flist.tmp', &
        iostat = io_status)
    if (io_status == 0) then
        do
            read(udf, '(a256)', iostat = io_status) dataline
            if (io_status /= 0) exit
            if (len_trim(dataline) == 0) cycle
            !> The listing names itself, and the ready flag a prefetch leaves.
            if (index(dataline, '.tmp') /= 0) cycle
            if (index(dataline, 'ready.flag') /= 0) cycle

            if (index(dataline, '-biomet.' // trim(adjustl(MetaExt))) /= 0) then
                BiometMetaFile = dataline(1:len_trim(dataline))
                call StripFileName(BiometMetaFile)
            else if (index(dataline, '-biomet.' // trim(adjustl(DataExt))) /= 0) then
                BiometFile = dataline(1:len_trim(dataline))
                call StripFileName(BiometFile)
            else if (index(dataline, '.' // trim(adjustl(MetaExt))) /= 0) then
                MetaFile = dataline(1:len_trim(dataline))
                call StripFileName(MetaFile)
            else if (index(dataline, '.' // trim(adjustl(DataExt))) /= 0) then
                DataFile = dataline(1:len_trim(dataline))
                call StripFileName(DataFile)
            end if
        end do
    end if
    close(udf, status = 'delete')

    TmpString = MetaFile
    call basename(TmpString, MetaFile, slash)
    TmpString = BiometMetaFile
    call basename(TmpString, BiometMetaFile, slash)
    TmpString = DataFile
    call basename(TmpString, DataFile, slash)
    TmpString = BiometFile
    call basename(TmpString, BiometFile, slash)
end subroutine UnZipArchive
