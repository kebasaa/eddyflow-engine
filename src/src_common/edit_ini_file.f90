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


!***************************************************************************
!
! \brief       Set a tag in one [section] of an INI file, adding it at the end
!              of that section when the file does not have it yet.
! \author      Jonathan Muller
! \note        EditIniFile only rewrites a line that is already there, which
!              is right for a setting every project writes. A per-gas record
!              key whose absence means something - gas_<i>_sa_source absent
!              is "automatic" - is missing exactly when it would be set, so
!              EditIniFile would drop the change without a word. The section
!              matters for the same reason: a tag added after the last group
!              would be read as part of whichever group came last.
!
!              Same temp-file discipline as EditIniFile: both files are
!              checked open before anything is written. Nothing is added when
!              the section itself is absent.
!***************************************************************************
subroutine SetIniKeyInSection(fname, section, tag, newval)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: fname
    character(*), intent(in) :: section
    character(*), intent(in) :: tag
    character(*), intent(in) :: newval
    !> Local variables
    integer :: ierr
    integer :: sepa
    integer :: io_error
    integer :: ini_in
    integer :: ini_out
    logical :: in_section
    logical :: done
    character(PathLen) :: tfname
    character(ShortInstringLen) :: dataline
    character(ShortInstringLen) :: trimmed
    character(64) :: currtag

    open(newunit = ini_in, file = trim(fname), status = 'old', action = 'read', &
        iostat = io_error)
    if (io_error /= 0) return
    tfname = trim(fname) // '.tmp'
    open(newunit = ini_out, file = trim(tfname), status = 'new', iostat = io_error)
    if (io_error /= 0) then
        close(ini_in)
        return
    end if

    in_section = .false.
    done = .false.
    do
        read(ini_in, '(a)', iostat = ierr) dataline
        if (ierr < 0) exit
        trimmed = adjustl(dataline)
        if (trimmed(1:1) == '[') then
            !> Leaving the section without having met the tag: add it as the
            !> section's last line, before the next header.
            if (in_section .and. .not. done) then
                write(ini_out, '(a)') trim(tag) // '=' // trim(newval)
                done = .true.
            end if
            in_section = trim(trimmed) == '[' // trim(section) // ']'
        else if (in_section .and. .not. done) then
            sepa = index(trimmed, '=')
            if (sepa > 1) then
                currtag = trimmed(1:sepa - 1)
                if (trim(currtag) == trim(tag)) then
                    trimmed = trim(tag) // '=' // trim(newval)
                    done = .true.
                end if
            end if
        end if
        write(ini_out, '(a)') trim(trimmed)
    end do
    !> The section ran to the end of the file
    if (in_section .and. .not. done) write(ini_out, '(a)') trim(tag) // '=' // trim(newval)
    close(ini_in, status = 'DELETE')
    close(ini_out)
    call rename(trim(adjustl(tfname)), trim(adjustl(fname)), status = io_error)
end subroutine SetIniKeyInSection
