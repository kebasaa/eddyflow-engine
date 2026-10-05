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
    integer :: n
    real(kind = dbl), allocatable :: bAggrOut(:)
    type(DeferredFluxnetRowType), allocatable :: grown(:)


    !> This period's biomet values, in the units the header names them in.
    if (nbVars > 0) then
        if (EddyFlowProj%fluxnet_standardize_biomet) then
            if (allocated(bAggrFluxnet)) bAggrOut = bAggrFluxnet
        else
            if (allocated(bAggr)) bAggrOut = bAggr
        end if
    end if
    if (.not. allocated(bAggrOut)) allocate(bAggrOut(0))

    !> This period's own day or night. Stats%daytime is otherwise set only for
    !> a period that gets as far as AssessDaytime, so a skipped one carried the
    !> NIGHT flag of the last period processed before it. Its biomet values are
    !> this period's (error where there are none), and without them the
    !> timestamp's potential radiation decides, as for any period.
    call AssessDaytime(Stats%date, Stats%time)

    if (FluxnetFileOpen) then
        call WriteFluxnetOnlyBiometRow(Stats%start_date, Stats%start_time, &
            Stats%date, Stats%time, Stats%daytime, bAggrOut, size(bAggrOut))
        return
    end if

    !> No header yet, so no row width either: keep what is particular to this
    !> period and let InitFluxnetFile_rp write the row once it has written the
    !> header. Writing it now would go to a unit nobody has opened, which
    !> gfortran quietly turns into fort.132 in the working directory.
    if (.not. allocated(DeferredFluxnetRows)) allocate(DeferredFluxnetRows(64))
    if (nDeferredFluxnetRows == size(DeferredFluxnetRows)) then
        allocate(grown(2 * size(DeferredFluxnetRows)))
        grown(1:nDeferredFluxnetRows) = DeferredFluxnetRows
        call move_alloc(grown, DeferredFluxnetRows)
    end if
    n = nDeferredFluxnetRows + 1
    DeferredFluxnetRows(n)%start_date = Stats%start_date
    DeferredFluxnetRows(n)%start_time = Stats%start_time
    DeferredFluxnetRows(n)%date = Stats%date
    DeferredFluxnetRows(n)%time = Stats%time
    DeferredFluxnetRows(n)%daytime = Stats%daytime
    DeferredFluxnetRows(n)%biomet = bAggrOut
    nDeferredFluxnetRows = n

end subroutine WriteOutFluxnetOnlyBiomet

!***************************************************************************
!
! \brief       Write the rows of periods skipped before the FLUXNET header
!              existed, in the order they were skipped
! \note        Called by InitFluxnetFile_rp straight after the header, so the
!              rows land ahead of the first period that had data.
!***************************************************************************
subroutine FlushDeferredFluxnetRows()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i


    do i = 1, nDeferredFluxnetRows
        call WriteFluxnetOnlyBiometRow(DeferredFluxnetRows(i)%start_date, &
            DeferredFluxnetRows(i)%start_time, DeferredFluxnetRows(i)%date, &
            DeferredFluxnetRows(i)%time, DeferredFluxnetRows(i)%daytime, &
            DeferredFluxnetRows(i)%biomet, size(DeferredFluxnetRows(i)%biomet))
    end do
    nDeferredFluxnetRows = 0
    if (allocated(DeferredFluxnetRows)) deallocate(DeferredFluxnetRows)

end subroutine FlushDeferredFluxnetRows

!***************************************************************************
!
! \brief       Build and write one skipped period's FLUXNET row
! \note        Everything that sets the row's shape - the fixed width, the
!              custom variables, the gas, analyser, hygrometer and CEC blocks,
!              the biomet count - is read here, as the header has it. Only
!              what belongs to the period itself comes in as arguments, which
!              is what lets a row be written later than its period.
!              bvals holds nb values; the header's nbVars columns past those
!              are written as errors.
!***************************************************************************
subroutine WriteFluxnetOnlyBiometRow(start_date, start_time, end_date, &
    end_time, is_daytime, bvals, nb)
    use m_rp_global_var
    use m_cec
    implicit none
    !> in/out variables
    character(*), intent(in) :: start_date
    character(*), intent(in) :: start_time
    character(*), intent(in) :: end_date
    character(*), intent(in) :: end_time
    logical, intent(in) :: is_daytime
    integer, intent(in) :: nb
    real(kind = dbl), intent(in) :: bvals(nb)
    !> local variables
    integer :: i
    integer :: indx
    integer :: int_doy
    real(kind = dbl) :: float_doy
    character(32) :: char_doy
    character(LongOutstringLen) :: csv_row
    character(14) :: tsIso
    real(kind = dbl) :: lrad
    integer :: j
    integer :: cec_p
    integer :: cec_k
    integer :: n_cec_pairs
    integer :: n_cec_fields
    integer :: n_w_slots
    integer :: w_slots(GHGNumVar)
    character(32) :: w_tags(GHGNumVar)
    type(CECResolvedPairType) :: cec_pairs(MaxNumCecPairs)
    type(CECDescriptorType) :: cec_blank
    real(kind = dbl) :: cec_values(MaxNumCecTargets * 6 + 9)
    logical :: cec_is_int(MaxNumCecTargets * 6 + 9)
    include '../src_common/interfaces.inc'

    call clearstr(csv_row)

    !> Start/end imestamps
    tsIso = start_date(1:4) // start_date(6:7) // start_date(9:10) &
                // start_time(1:2) // start_time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)
    tsIso = end_date(1:4) // end_date(6:7) // end_date(9:10) &
                // end_time(1:2) // end_time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)

    !> DOYs
    !>  Start
    call DateTimeToDOY(start_date, start_time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)
    !>  End
    call DateTimeToDOY(end_date, end_time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)

    !> Not enough data
    call AddDatum(csv_row, 'not_enough_data', separator)

    !> Potential Radiations
    indx = DateTimeToHalfHourNumber(end_date, end_time) - 1
    indx = max(indx, 2)
    lrad = (PotRad(indx) + PotRad(indx - 1)) / 2
    call AddFloatDatumToDataline(lrad, csv_row, EddyFlowProj%err_label)

    !> Daytime
    if (is_daytime) then
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

    !> Error codes in place of the custom variables and their count: one more
    !> than the header names, not than this period declares. Where every raw
    !> file brings its own metadata, as GHG archives do, a skipped period's file
    !> can declare fewer variables than the one the header was written from:
    !> base_ghg_mixed_60's 02:00 period declared 13 against the header's 16,
    !> and its row came out three columns short, every column from there on
    !> under the wrong name.
    do i = 1, nFluxnetCustomVars + 1
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end do

    !> Per-gas moisture, analyser and hygrometer blocks, same position and
    !> width as in WriteOutFluxnet: the slot is real so the row stays
    !> parseable, the values are error because this period was skipped.
    !>
    !> The widths here were three, fourteen and nothing at all, against the
    !> seven, fifteen and one-plus-seven-per-hygrometer the header names and
    !> ReadExRecord steps over. It went unnoticed because the CEC descriptor
    !> used to be found by counting back from the end of the row, so the
    !> shortfall was quietly absorbed into the biomet chunk. Nothing is
    !> end-anchored now, and a short row is a misread row.
    call AddIntDatumToDataline(nFluxnetGasSlots, csv_row, EddyFlowProj%err_label)
    do i = 1, nFluxnetGasSlots
        call AddIntDatumToDataline(FluxnetGasSlots(i), csv_row, EddyFlowProj%err_label)
        do indx = 1, 6
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end do
    call AddIntDatumToDataline(nFluxnetInstrSlots, csv_row, EddyFlowProj%err_label)
    do i = 1, nFluxnetInstrSlots
        call AddIntDatumToDataline(FluxnetInstrSlots(i), csv_row, EddyFlowProj%err_label)
        do indx = 1, 14
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end do

    call WaterOutSlots(w_slots, w_tags, n_w_slots)
    j = 0
    do i = 1, n_w_slots
        if (len_trim(w_tags(i)) > 0) j = j + 1
    end do
    call AddIntDatumToDataline(j, csv_row, EddyFlowProj%err_label)
    do i = 1, n_w_slots
        if (len_trim(w_tags(i)) == 0) cycle
        call AddIntDatumToDataline(w_slots(i), csv_row, EddyFlowProj%err_label)
        do indx = 1, 6
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end do

    !> The CEC descriptors, at their full width and all error: a skipped period
    !> has no partition, but it still has to occupy its own columns.
    call CecPairs(cec_pairs, n_cec_pairs)
    call AddIntDatumToDataline(n_cec_pairs, csv_row, EddyFlowProj%err_label)
    call ResetCecDescriptor(cec_blank)
    do cec_p = 1, n_cec_pairs
        call CecExRowValues(cec_pairs(cec_p), cec_blank, cec_values, &
            cec_is_int, n_cec_fields)
        do cec_k = 1, n_cec_fields
            if (cec_is_int(cec_k)) then
                call AddIntDatumToDataline(nint(cec_values(cec_k)), csv_row, &
                    EddyFlowProj%err_label)
            else
                call AddFloatDatumToDataline(cec_values(cec_k), csv_row, &
                    EddyFlowProj%err_label)
            end if
        end do
    end do

    !> write all aggregated biomet values in FLUXNET units
    !> A period skipped before embedded biomet was first read has no values;
    !> it still takes the header's count, with errors in its columns.
    call AddIntDatumToDataline(nbVars, csv_row, EddyFlowProj%err_label)
    do i = 1, nbVars
        if (i <= nb) then
            call AddFloatDatumToDataline(bvals(i), csv_row, EddyFlowProj%err_label)
        else
            call AddFloatDatumToDataline(error, csv_row, EddyFlowProj%err_label)
        end if
    end do
    write(uflxnt, '(a)') csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteFluxnetOnlyBiometRow
