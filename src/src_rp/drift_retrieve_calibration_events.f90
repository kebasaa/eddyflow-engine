!***************************************************************************
! drift_retrieve_calibration_events.f90
! -------------------------------------
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
!*******************************************************************************
!
! \brief       Reads external dynamic (time-varying) metadata file and retrieve
!              cleaning and calibration events information
!              alternative metadata file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine driftRetrieveCalibrationEvents(nCalibEvents)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: nCalibEvents
    !> local variables
    !> Four scalar columns, and two per-gas column families. The families
    !> used to share one array with the gases packed at 5..8 and their
    !> references at 9..12, read back through `Calib(i)%ref(j - 4)`. That
    !> arithmetic is why there could be exactly four: the two index spaces are
    !> separate now, so no offset relates them and none is needed.
    integer :: cal_date
    integer :: cal_time
    integer :: cal_t
    integer :: cleaning
    integer :: offcol(GHGNumVar)  !< header column of <gas>_offset, 0 = absent
    integer :: refcol(GHGNumVar)  !< header column of <gas>_ref,    0 = absent
    integer :: open_status
    integer :: read_status
    integer :: cnt
    integer :: i
    integer :: gas
    integer :: slot
    integer :: var_num
    integer :: sepa
    real(kind = dbl) :: biased_abs
    real(kind = dbl) :: unbiased_abs
    real(kind = dbl) :: abs_scale
    character(32) :: text_vars(128)
    character(64) :: field
    character(LongInstringLen) :: dataline
    logical :: bias_is_negative
    integer, external :: GasSlotFromDynMDTag
    logical, external :: GasSlotIsWater


    !> Open dynamic metadata file
    write(*,'(a)') ' Looking for calibration information &
        &in dynamic metadata file:  '
    write(*,'(a)') '  ' // AuxFile%DynMD(1:len_trim(AuxFile%DynMD))
    open(udf, file = AuxFile%DynMD, status = 'old', iostat = open_status)

    Calib = CalibType('', '', tsNull, .false., nint(error), &
        error, error, error, error, error, error)

    nCalibEvents = 0
    cal_date = 0
    cal_time = 0
    cal_t    = 0
    cleaning = 0
    offcol = 0
    refcol = 0
    if (open_status == 0) then
        !> Look in the header for calibration information. Every configured
        !> gas is offered its own pair of columns, named as the record names
        !> it; the four historical spellings stay valid as aliases for the
        !> first four slots, since that is what files in the wild carry.
        read(udf, '(a)', iostat = read_status) dataline
        cnt = 0
        do
            sepa = index(dataline, ',')
            if (sepa == 0) sepa = len_trim(dataline) + 1
            if (len_trim(dataline) == 0) exit
            cnt = cnt + 1
            field = dataline(1:sepa - 1)
            !> Date and time
            if (trim(field) &
                == trim(adjustl(StdDynMDVars(dynmd_date)))) cal_date = cnt
            if (trim(field) &
                == trim(adjustl(StdDynMDVars(dynmd_time)))) cal_time = cnt
            !> calibration data columns
            slot = GasSlotFromDynMDTag(field, '_offset')
            if (slot > 0) offcol(slot) = cnt
            slot = GasSlotFromDynMDTag(field, '_ref')
            if (slot > 0) refcol(slot) = cnt
            if (trim(field) == 'calibration_temperature') cal_t = cnt
            if (trim(field) == 'cleaning') cleaning = cnt

            dataline = dataline(sepa + 1: len_trim(dataline))
        end do

        if (all(offcol == 0)) then
            close(udf)
            return
        end if

        !> Reads file and store available calibration data
        i = 0
        do
            read(udf, '(a)', iostat = read_status) dataline
            if (read_status /= 0) exit
            i = i + 1

            !> For each data line, store data in a temporary array as text
            text_vars = 'none'
            var_num = 0
            do
                sepa = index(dataline, separator)
                if (sepa == 0) sepa = len_trim(dataline) + 1
                if (len_trim(dataline) == 0) exit
                var_num = var_num + 1
                text_vars(var_num) = dataline(1:sepa - 1)
                dataline = dataline(sepa + 1: len_trim(dataline))
            end do

            !> Associate stored data to relevant variables, for current dataline
            if (cleaning /= 0 .and. text_vars(cleaning) == '1') &
                Calib(i)%cleaning = .true.
            if (cal_date /= 0) read(text_vars(cal_date), *) Calib(i)%date
            if (cal_time /= 0) read(text_vars(cal_time), *) Calib(i)%time
            if (cal_date /= 0 .and. cal_time /= 0) &
                call DateTimeToDateType(Calib(i)%date, &
                    Calib(i)%time, Calib(i)%ts)

            if (cal_t /= 0) read(text_vars(cal_t), *) Calib(i)%Tcell
            do gas = firstGas, lastGas
                if (offcol(gas) /= 0) &
                    read(text_vars(offcol(gas)), *) Calib(i)%offset(gas)
                if (refcol(gas) /= 0) &
                    read(text_vars(refcol(gas)), *) Calib(i)%ref(gas)
            end do
        end do
        nCalibEvents = i
    end if
    close(udf)


    !> Convert offsets into absorptance offsets, considering the error as a
    !> span error thus, evaluating the abs_offset on the cal curve starting
    !> from the reference concentration indicated by the user
    !> Two unrolled blocks before, co2 and h2o, differing only in the 1d3 that
    !> puts water on the mmol basis. That is a property of the species, so it
    !> is asked of the record; every other gas takes the umol arm, which is
    !> what the co2 block was.
    !>
    !> A gas with no inverse calibration polynomial is left alone. The
    !> coefficients are per-channel instrument calibrations - there is no
    !> general form for an arbitrary species on an arbitrary analyser - so
    !> drift is opt-in per gas, and `error` is the opt-out DriftCorr is
    !> initialised to.
    do i = 1, nCalibEvents
        do gas = firstGas, lastGas
            if (offcol(gas) == 0 .or. refcol(gas) == 0) cycle
            if (DriftCorr%inv_cal(0, gas) == error) cycle
            if (Calib(i)%Tcell == error .or. Calib(i)%Tcell <= 0d0) cycle

            abs_scale = 1d0
            if (GasSlotIsWater(gas)) abs_scale = 1d3

            bias_is_negative = &
                Calib(i)%ref(gas) + Calib(i)%offset(gas) < 0d0

            call PolyVal(DriftCorr%inv_cal(0:6, gas), 6, &
                dabs(Calib(i)%ref(gas) + Calib(i)%offset(gas)) &
                / Calib(i)%Tcell / Ru * abs_scale, &
                1, biased_abs)
            if (bias_is_negative) biased_abs = - biased_abs

            call PolyVal(DriftCorr%inv_cal(0:6, gas), 6, &
                Calib(i)%ref(gas) / Calib(i)%Tcell / Ru * abs_scale, &
                1, unbiased_abs)

            Calib(i)%offset(gas) = biased_abs - unbiased_abs
        end do
    end do
    write(*,'(a)') ' Done.'
end subroutine driftRetrieveCalibrationEvents
