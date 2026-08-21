!***************************************************************************
! abort_on_missing_path.f90
! -------------------------
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
! \brief       Stop, because the project names a path that is not there.
! \author      Jonathan Muller
! \note        A project that states a path which does not exist used to be
!              answered by quietly computing something else: a missing planar
!              fit file became double rotation, a missing time-lag file became
!              covariance maximisation, and a missing assessment file or
!              co-spectra directory became Moncrieff 1997. The run finished and
!              the numbers were labelled as the method that had been asked for,
!              which is the part that makes it dangerous - a whole dataset can
!              be processed the wrong way and look right.
!
!              The engine already stops for three other missing inputs - the
!              project file, the alternative metadata file and the raw data
!              directory - and this is the same rule applied to the rest.
!
!              Only a path the project *states* reaches here. Nothing changes
!              for a method that was never selected: each caller sits inside a
!              branch that already means "the user chose the file option and
!              named a file".
!
!              The remedy is passed in rather than composed here because only
!              the caller knows which method was selected and what the
!              alternative is. Stopping without saying what to change just
!              moves the dead end.
!
!              The remedy arrives as one string and is wrapped here. It began
!              as an array of lines, which compiled cleanly and printed
!              nothing: an assumed-shape dummy needs an explicit interface, and
!              these callers reach an external procedure without one, so the
!              loop over size() ran zero times and the advice - the whole point
!              of the change - vanished silently.
!
! \sa          exception_handler.f90, cases 21, 22 and 86
!***************************************************************************
subroutine AbortOnMissingPath(setting, path, remedy)
    use m_index_parameters
    use m_log
    implicit none
    !> in/out variables
    !> The key as it is spelled in the project file, so it can be searched for
    character(*), intent(in) :: setting
    !> The path as the project gave it, unmodified - a trailing space or a
    !> stale drive letter is exactly what the reader needs to see
    character(*), intent(in) :: path
    !> What to do about it, as one sentence or several
    character(*), intent(in) :: remedy
    !> local variables
    integer, parameter :: wrap = 66
    integer :: first
    integer :: last
    integer :: brk
    character(*), parameter :: tag = ' Fatal error(111)> '

    call LogSayList(tag // 'The project setting "' // &
        trim(adjustl(setting)) // '" names a path that does not exist:')
    call LogSayList(tag // '  ' // trim(path))
    call LogSayList(trim(tag))

    !> Wrap on spaces so a long remedy stays inside a terminal, and so the
    !> callers can write one readable sentence instead of counting columns.
    first = 1
    do while (first <= len_trim(remedy))
        if (len_trim(remedy) - first + 1 <= wrap) then
            last = len_trim(remedy)
        else
            brk = index(remedy(first : first + wrap), ' ', back = .true.)
            if (brk == 0) then
                last = first + wrap - 1
            else
                last = first + brk - 2
            end if
        end if
        call LogSayList(tag // remedy(first : last))
        first = last + 2
    end do

    call LogSayList(tag // 'Program execution aborted.')
    stop 1
end subroutine AbortOnMissingPath
