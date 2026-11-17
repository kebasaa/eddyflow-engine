!***************************************************************************
! read_planar_fit_file.f90
! ------------------------
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
! \brief       Read planar fit file and import rotation matrices
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadPlanarFitFile()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: io_status
    integer :: open_status
    integer :: sec
    integer :: i
    integer :: zero
    integer :: j
    integer :: start
    character(ShortInstringLen) :: dataline
    character(64) :: strg


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading planar-fit file: ' // AuxFile%pf(1:len_trim(AuxFile%pf))
    write(ulog,'(a)') ' Reading planar-fit file: ' // AuxFile%pf(1:len_trim(AuxFile%pf))
    open(udf, file = AuxFile%pf, status = 'old', iostat = open_status)

    if (open_status == 0) then
        call LogSay(' planar fit file found, reading rotation matrices..')

        !> A moving-window file (after innFLUX) carries its own layout: the
        !> windows it was fitted with are applied, not the project's.
        if (ReadMovingWindows()) then
            close(udf)
            call LogSay(' Done.')
            return
        end if
        rewind(udf)

        !> Read number of sectors from relevant line
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'Number_of_selected_wind_sectors') == 0) cycle
            read(dataline(index(dataline, ':') + 1: len_trim(dataline)), *) PFSetup%num_sec
            exit
        end do

        !> Skip remaining lines until beginning of fitting plane coefficients
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'WindSector') == 0) cycle
            exit
        end do

        !> Read fitting plane coefficients
        do sec = 1, PFSetup%num_sec
            call clearstr(strg)
            read(udf, '(a)', iostat = io_status) strg
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            start = index(strg, '-') + 6
            read(strg(start:), *) (PFb(i, sec), i = 1, 3)
        end do

        !> Skip remaining lines until beginning of rotation matrices
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'Rotation matrices') == 0) cycle
            exit
        end do

        !> Read rotation matrices
        do sec = 1, PFSetup%num_sec
            call clearstr(strg)
            read(udf, '(a)', iostat = io_status) strg
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if

            if(sec == 1) then
                !> From first sector retrieve north offset
                read(strg(26:30), '(i5)', iostat = io_status) zero
                if (zero >= 0) then
                    PFSetup%north_offset = zero
                else
                    PFSetup%north_offset = zero - 360d0
                end if
            end if

            !> Read end-of-sector angle
            read(strg(32:35), '(i4)', iostat = io_status) PFSetup%wsect_end(sec)
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            !> Back to the frame the run compares in. The file states sector
            !> ends as true bearings - the writer adds the north offset - while
            !> WindSector compares them with the wind direction minus that
            !> offset, as a run that fits the sectors itself stores them. Read
            !> as written, every boundary moved by the offset again and the
            !> wind in that band got its neighbour's rotation; with a negative
            !> offset some directions fell in no sector at all. The last
            !> sector closes the circle at 360.
            PFSetup%wsect_end(sec) = modulo(PFSetup%wsect_end(sec) - zero, 360)
            if (PFSetup%wsect_end(sec) == 0) PFSetup%wsect_end(sec) = 360
            !> Determine whether sector is to be excluded
            PFSetup%wsect_exclude(sec) = index(strg,'excluded') /= 0

            inloop: do i = 1, 3
                read(udf, *, iostat = io_status) (PFMat(i, j, sec), j = 1, 3)
                if (io_status /= 0) then
                    Meth%rot = 'double_rotation'
                    call ExceptionHandler(29)
                    return
                end if
            enddo inloop
        end do
        write(LogInteger, '(i6)') PFSetup%num_sec
        write(*,'(a)') '  ' // trim(adjustl(LogInteger)) // ' sector(s) found.'
        write(ulog,'(a)') '  ' // trim(adjustl(LogInteger)) // ' sector(s) found.'
    else
        !> The project asked for planar fit and named a file that is not there.
        !> This used to fall back to double rotation and carry on, which moves
        !> every flux in the run and says so only in an alert nobody reads.
        call AbortOnMissingPath('pf_file', AuxFile%pf, &
            'Correct the path to the planar fit file, or choose ' &
            // '"Planar fit file not available" so the rotation matrices are ' &
            // 'computed from this dataset instead of being read.')
    end if
    call LogSay(' Done.')
contains

    !> True when the file is a moving-window file and was read; the
    !> project's own pf_mw_* settings then give way to the file's.
    logical function ReadMovingWindows()
        integer :: n
        integer :: k
        integer :: ios
        integer :: nmw
        integer :: nb
        integer :: seedf
        integer :: ne
        real(kind = dbl) :: width
        real(kind = dbl) :: centre
        real(kind = dbl) :: vals(14)
        character(ShortInstringLen) :: line
        character(32) :: word

        ReadMovingWindows = .false.
        nmw = 0
        do
            read(udf, '(a)', iostat = ios) line
            if (ios /= 0) return
            if (index(line, 'Moving_window_sectors:') /= 0) then
                read(line(index(line, ':') + 1:), *, iostat = ios) nmw, word, width
                if (ios /= 0 .or. nmw <= 0 .or. nmw > MaxNumPfSect) return
                exit
            end if
            if (index(line, 'Rotation matrices') /= 0) return
        end do
        do
            read(udf, '(a)', iostat = ios) line
            if (ios /= 0) return
            if (index(line, 'North_offset:') /= 0) &
                read(line(index(line, ':') + 1:), *, iostat = ios) PFSetup%north_offset
            if (index(line, 'Bootstrap:') /= 0) then
                read(line(index(line, ':') + 1:), *, iostat = ios) nb, word, seedf
                if (ios == 0) then
                    PFSetup%boot_n = nb
                    PFSetup%boot_seed = seedf
                end if
            end if
            if (index(line, 'window,centre_deg') /= 0) exit
        end do
        if (PFSetup%mw_num > 0 .and. (PFSetup%mw_num /= nmw &
            .or. abs(PFSetup%mw_width - width) > 1d-6)) &
            call LogSay('  Note: the planar-fit file''s moving windows differ from ' &
                // 'the project''s pf_mw_* settings; the file''s are applied.')
        PFSetup%mw_num = nmw
        PFSetup%mw_width = width
        PFSetup%num_sec = nmw
        do n = 1, nmw
            read(udf, *, iostat = ios) k, centre, ne, vals
            if (ios /= 0 .or. k /= n) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                ReadMovingWindows = .true.
                return
            end if
            PFb(1:3, n) = vals(1:3)
            PFStdAlpha(n) = vals(4)
            PFStdBeta(n) = vals(5)
            PFMat(1, 1:3, n) = vals(6:8)
            PFMat(2, 1:3, n) = vals(9:11)
            PFMat(3, 1:3, n) = vals(12:14)
        end do
        write(LogInteger, '(i6)') nmw
        call LogSay('  ' // trim(adjustl(LogInteger)) // ' moving window(s) found.')
        ReadMovingWindows = .true.
    end function ReadMovingWindows
end subroutine ReadPlanarFitFile
