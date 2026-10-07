!***************************************************************************
! m_log.f90
! ---------
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
! \brief       Everything the engine prints, to the console and to a file.
! \author      Jonathan Muller, ETH Zurich
! \note        A run leaves a dozen output files behind and no record of what
!              the engine said while producing them - which warnings fired,
!              which periods were skipped and why, what the settings resolved
!              to. That record only existed in whatever terminal the run
!              happened in, and under the interface it did not exist at all.
!
!              Fortran cannot duplicate the default output unit, so this is a
!              tee rather than a redirect: every message goes to the console
!              exactly as it did before AND to unit ulog. Reopening unit 6 onto
!              a file would have been one line, but it would have silenced the
!              console, and the interface reads the engine's standard output to
!              drive its progress display.
!
!              **ulog is connected from the program's first statement**, to an
!              unnamed scratch file, because output starts long before the log
!              can have a name: the banner, the project file being read, and any
!              exception raised while reading it all happen before the output
!              directory is known. LogInit copies the scratch across once the
!              path exists and carries on in the real file. Connecting it late
!              instead would have left every one of those early lines writing to
!              a stray fort.163 in the working directory - which is exactly what
!              it did, until this was fixed.
! \sa          init_run_log.f90
!***************************************************************************
module m_log
    use m_index_parameters
    use iso_fortran_env, only: iostat_end, output_unit
    implicit none
    save
    private

    public :: LogStart, LogInit, LogSay, LogSayList, LogSayNoAdv, LogClose, &
              LogIsOpen, LogFlush, LogOpenLine, LogEndLine

    !> Long enough for the widest message the engine emits.
    integer, parameter :: LogLineLen = 2048
    !> A spare unit for the copy out of the scratch file. Only ever open inside
    !> LogInit, so it cannot collide with anything the engine holds.
    integer, parameter :: uLogCopy = 164

    logical :: Connected = .false.
    logical :: Named = .false.
    !> A progress line opened by LogOpenLine and not yet ended. A message
    !> said meanwhile ends it first, so the message keeps its own line - it
    !> used to be written onto the progress text, where the interface took
    !> the whole line for progress and dropped the message's first line.
    logical :: LineOpen = .false.

contains

    !> Whether the log has reached its real file yet.
    logical function LogIsOpen()
        LogIsOpen = Named
    end function LogIsOpen

    !***********************************************************************
    !> Connect the log before anything is said. The first executable statement
    !> of each binary.
    !***********************************************************************
    subroutine LogStart()
        integer :: open_status

        if (Connected) return
        open(ulog, status = 'scratch', iostat = open_status)
        Connected = open_status == 0
    end subroutine LogStart

    !***********************************************************************
    !> Give the log its name, once there is an output folder to put it in.
    !>
    !> Everything said so far is in the scratch file, in order; it is copied
    !> across and the scratch discarded. Silent on failure - an unwritable
    !> output folder is a problem the rest of the run reports far better than
    !> this can, and losing the log is not a reason to stop processing.
    !***********************************************************************
    subroutine LogInit(path)
        character(*), intent(in) :: path
        character(LogLineLen) :: line
        integer :: open_status
        integer :: read_status
        integer :: nchars

        if (Named .or. .not. Connected) return

        open(uLogCopy, file = path, status = 'replace', iostat = open_status, &
            encoding = 'utf-8')
        if (open_status /= 0) return

        !> Read non-advancing with size=, so the record's own length comes back
        !> rather than being guessed at. trim() would have been wrong: reading
        !> into a fixed buffer pads with blanks, and trimming cannot tell that
        !> padding from a message that really ends in a space - " Executing
        !> EddyFlow " does, and lost it.
        !>
        !> End-of-record and end-of-file are both negative, so they have to be
        !> told apart by name: treating any negative as the end would stop at
        !> the first blank line, of which the preamble has several.
        rewind(ulog)
        do
            read(ulog, '(a)', advance = 'no', size = nchars, &
                iostat = read_status) line
            if (read_status > 0) exit
            if (read_status == iostat_end .and. nchars == 0) exit
            write(uLogCopy, '(a)') line(1:nchars)
            if (read_status == iostat_end) exit
        end do
        close(uLogCopy)

        !> Closing a scratch file deletes it, which is the whole point of using
        !> one: nothing is left behind in whatever directory the run started in.
        close(ulog)
        open(ulog, file = path, status = 'old', position = 'append', &
            iostat = open_status, encoding = 'utf-8')
        if (open_status /= 0) then
            Connected = .false.
            return
        end if
        Named = .true.
    end subroutine LogInit

    !***********************************************************************
    !> One line, to the console and to the log.
    !>
    !> The console write is the one that was here before, unchanged, because
    !> the interface parses this stream and the regression harness compares it
    !> line for line.
    !***********************************************************************
    subroutine LogSay(text)
        character(*), intent(in) :: text

        !> LogSay('') is how a caller ends an open line before a message; on
        !> a line opened here that is all it has to do.
        if (LineOpen .and. len_trim(text) == 0) then
            call EndOpenLine()
            call LogFlush()
            return
        end if
        call EndOpenLine()
        write(*, '(a)') text
        if (Connected) write(ulog, '(a)') text
        call LogFlush()
    end subroutine LogSay

    !***********************************************************************
    !> One line, list-directed.
    !>
    !> Kept apart from LogSay because the two do NOT render alike: list-directed
    !> output opens a record with a blank, so `write(*,*) 'x'` prints " x" where
    !> `write(*,'(a)') 'x'` prints "x". Most of the engine's messages - every
    !> line of the exception handler among them - were written the first way,
    !> and quietly losing that leading blank would change every one of them.
    !***********************************************************************
    subroutine LogSayList(text)
        character(*), intent(in) :: text

        call EndOpenLine()
        write(*, *) text
        if (Connected) write(ulog, *) text
        call LogFlush()
    end subroutine LogSayList

    !***********************************************************************
    !> The same, without ending the line - the progress messages that print a
    !> label, work, and then print " Done." onto the same line.
    !***********************************************************************
    subroutine LogSayNoAdv(text)
        character(*), intent(in) :: text

        write(*, '(a)', advance = 'no') text
        if (Connected) write(ulog, '(a)', advance = 'no') text
        call LogFlush()
    end subroutine LogSayNoAdv

    !> Flushed per line once the log has a name, because the run that most
    !> needs reading back is the one that died, and a buffered tail is exactly
    !> what would be missing from it. Pointless on the scratch file, which
    !> LogInit reads back in full anyway.
    !>
    !> The console too. When it is not a terminal - the interface reading a
    !> pipe, a pre-pass worker writing to its .out file - gfortran buffers it
    !> in blocks, and a worker's file stayed empty for an hour of a parallel
    !> run while its parent's progress arrived in bursts. A handful of lines a
    !> period costs nothing to flush.
    !> Starts a progress line - "  Calculating storage terms.." - that
    !> LogEndLine finishes. Only lines a message can interrupt need it; the
    !> rest keep LogSayNoAdv.
    subroutine LogOpenLine(text)
        character(*), intent(in) :: text

        call EndOpenLine()
        write(*, '(a)', advance = 'no') text
        if (Connected) write(ulog, '(a)', advance = 'no') text
        LineOpen = .true.
        call LogFlush()
    end subroutine LogOpenLine

    !> Finishes the progress line with `text` - " Done." - or, when a message
    !> has ended it meanwhile, writes `text` on a line of its own.
    subroutine LogEndLine(text)
        character(*), intent(in) :: text

        write(*, '(a)') text
        if (Connected) write(ulog, '(a)') text
        LineOpen = .false.
        call LogFlush()
    end subroutine LogEndLine

    subroutine EndOpenLine()
        if (.not. LineOpen) return
        write(*, '(a)') ''
        if (Connected) write(ulog, '(a)') ''
        LineOpen = .false.
    end subroutine EndOpenLine

    subroutine LogFlush()
        flush(output_unit)
        if (Named) flush(ulog)
    end subroutine LogFlush

    subroutine LogClose()
        if (.not. Connected) return
        close(ulog)
        Connected = .false.
        Named = .false.
    end subroutine LogClose

end module m_log
