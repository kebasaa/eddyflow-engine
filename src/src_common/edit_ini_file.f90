!***************************************************************************
! edit_ini_file.f90
! -----------------
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
! \brief       Modify INI file, changing the value for the specified tag
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EditIniFile(fname, tag, newval)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: fname
    character(*), intent(in) :: tag
    character(*), intent(in) :: newval
    !> Local variables
    integer :: ierr
    integer :: sepa
    integer :: io_error
    character(PathLen) :: tfname
    character(ShortInstringLen) :: dataline
    character(64) :: currtag


    !> Open file to be edited and temp file.
    !>
    !> Both opens are checked before anything is written, because the close
    !> below deletes the input. Ignoring the status meant a temp file left by
    !> an earlier run - status='new' fails on an existing one - was enough to
    !> destroy the file this routine was asked to edit, having written the
    !> replacement nowhere.
    open(10, file=trim(fname), iostat=io_error)
    if (io_error /= 0) return
    tfname = trim(fname) // '.tmp'
    open(11, file=trim(tfname), status='new', iostat=io_error)
    if (io_error /= 0) then
        close(10)
        return
    end if

    !> Copy whole file into temp file including modification
    do
        read(10, '(a)', iostat=ierr) dataline
        if (ierr < 0) exit
        sepa = index(dataline, '=')
        if (sepa /= 0) then
            currtag = dataline(1:sepa-1)
            if (currtag == trim(tag)) dataline = trim(tag) // '='// trim(newval)
        end if
        write(11, '(a)') trim(adjustl(dataline))
    end do
    close(10, status='DELETE')
    close(11)

    !> Input file has been deleted, now change name of tmp file into old input file name
    call rename(trim(adjustl(tfname)), trim(adjustl(fname)), status = io_error)
end subroutine EditIniFile

