!***************************************************************************
! write_out_fluxnet_only_biomet.f90
! ---------------------------------
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
! \brief       Write line to FLUXNET output with only biomet data, if available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnetOnlyBiomet()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    integer :: indx
    integer :: int_doy
    real(kind = dbl) :: float_doy
    character(32) :: char_doy
    character(LongOutstringLen) :: csv_row
    character(14) :: tsIso
    real(kind = dbl), allocatable :: bAggrOut(:)
    real(kind = dbl) :: lrad
    include '../src_common/interfaces.inc'

    call clearstr(csv_row)

    !> Start/end imestamps
    tsIso = Stats%start_date(1:4) // Stats%start_date(6:7) // Stats%start_date(9:10) &
                // Stats%start_time(1:2) // Stats%start_time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)
    tsIso = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)

    !> DOYs
    !>  Start
    call DateTimeToDOY(Stats%start_date, Stats%start_time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)
    !>  End
    call DateTimeToDOY(Stats%date, Stats%time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)

    !> Not enough data
    call AddDatum(csv_row, 'not_enough_data', separator)

    !> Potential Radiations
    indx = DateTimeToHalfHourNumber(Stats%date, Stats%time) - 1
    indx = max(indx, 2)
    lrad = (PotRad(indx) + PotRad(indx - 1)) / 2
    call AddFloatDatumToDataline(lrad, csv_row, EddyFlowProj%err_label)

    !> Daytime
    if (Stats%daytime) then
        call AddDatum(csv_row, '0', separator)
    else
        call AddDatum(csv_row, '1', separator)
    endif

    !> Write error codes in place of fixed columns.
    !>
    !> Padded to the width InitFluxnetFile_rp measured off the header it wrote,
    !> less the seven columns already written above. This was a literal 465,
    !> which had drifted four columns short of the header - every skipped
    !> period emitted a row that did not line up with its own header from the
    !> custom variables onward.
    do i = 1, nFluxnetFixedCols - 7
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end do

    !> Write error codes in place of custom variables
    do i = 1, NumUserVar + 1
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end do

    !> Per-gas moisture and analyser blocks, same position and width as in
    !> WriteOutFluxnet: the slot is real so the row stays parseable, the
    !> values are error because this period was skipped.
    call AddIntDatumToDataline(nFluxnetGasSlots, csv_row, EddyFlowProj%err_label)
    do i = 1, nFluxnetGasSlots
        call AddIntDatumToDataline(FluxnetGasSlots(i), csv_row, EddyFlowProj%err_label)
        do indx = 1, 2
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end do
    call AddIntDatumToDataline(nFluxnetInstrSlots, csv_row, EddyFlowProj%err_label)
    do i = 1, nFluxnetInstrSlots
        call AddIntDatumToDataline(FluxnetInstrSlots(i), csv_row, EddyFlowProj%err_label)
        do indx = 1, 13
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end do

    !> write all aggregated biomet values in FLUXNET units
    call AddIntDatumToDataline(nbVars, csv_row, EddyFlowProj%err_label)
    if (nbVars > 0) then
        if (.not. allocated(bAggrOut)) allocate(bAggrOut(size(bAggr)))
        if (EddyFlowProj%fluxnet_standardize_biomet) then
            bAggrOut = bAggrFluxnet
        else
            bAggrOut = bAggr
        end if

        do i = 1, nbVars
            call AddFloatDatumToDataline(bAggrOut(i), csv_row, EddyFlowProj%err_label)
        end do
    end if
    !> CEC ratios — error for skipped periods
    call AddFloatDatumToDataline(CECDescriptor%r_ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%r_Fc, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_valid, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_O1, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_O2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%frac_O1, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%frac_O2, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, CECDescriptor%h2o_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, CECDescriptor%co2_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%h2o_status, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%co2_status, csv_row, EddyFlowProj%err_label)

    write(uflxnt, '(a)') csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteOutFluxnetOnlyBiomet
