!***************************************************************************
! writeout_timelag_optimization.f90
! ---------------------------------
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
! \brief       Write time-lag optimization results on output file \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutTimelagOptimization(actn, M, h2o_n, ncls, cls_size)
    use m_rp_global_var
    use m_pwb_timelag, only: TimelagOptGasLabel
    implicit none
    !> in/out variables
    integer, intent(in) :: M
    integer, intent(in) :: ncls
    integer, intent(in) :: actn(M)
    real(kind = dbl), intent(in) :: cls_size
    integer, intent(out) :: h2o_n(ncls)
    integer, external :: CreateDir
    !> local variables
    integer :: cls
    integer :: gas
    integer :: open_status = 1
    real(kind = dbl) :: tl_def, tl_min, tl_max
    character(32) :: gname
    logical, external :: GasSlotIsWater
    character(4) :: min
    character(4) :: max
    character(9) :: txt


    !> Create output file
    TimelagOpt_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
              // EddyFlowProj%id(1:len_trim(EddyFlowProj%id))
    if (PwbAggregateSummary) then
        TimelagOpt_Path = trim(TimelagOpt_Path) // '_pwb_timelag_opt' // Timestamp_FilePadding // TxtExt
    else
        TimelagOpt_Path = trim(TimelagOpt_Path) // TimelagOpt_FilePadding // Timestamp_FilePadding // TxtExt
    end if
    open(uto, file = TimelagOpt_Path, iostat = open_status, encoding = 'utf-8')

    !> Write on output file time-lag optimization results
    write(uto, '(a)') 'Time-lag_optimisation_results'
    if (PwbAggregateSummary) write(uto, '(a)') 'PWB_aggregate_summary: true'
    write(uto, '(a, f7.2)') 'Plausibility_range_[timefolds_standard_deviation]:',TOSetup%pg_range
    write(uto, '(a, a)') 'Beginning_of_timelag_optimization_period: ', TOSetup%start_date
    write(uto, '(a, a)') 'End_of_timelag_optimization_period: ', TOSetup%end_date
    write(uto, '(a)')

    !> One block per configured gas, named from TimelagOptGasLabel so the
    !> reader looks each one up under the name it was written with. These were
    !> four hand-written blocks, so a gas past the fourth had no entry at all
    !> and its optimised window was lost between runs.
    !>
    !> A gas the optimiser found no determinations for is written as the error
    !> code, not as 0.00. Zero is a window - [0, 0] - and SetTimelags would
    !> take it in place of the metadata's declared one, detecting every lag as
    !> zero. Its guard is `max > min`, which the sentinel fails, so the
    !> declared window survives. Absent means not performed, as everywhere
    !> else in this file format.
    !>
    !> Water is skipped when it is classed by relative humidity; the per-class
    !> table below carries it instead, which is what the original `ncls <= 1`
    !> guard on the h2o block said.
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        if (GasSlotIsWater(gas) .and. ncls > 1) cycle
        gname = TimelagOptGasLabel(gas)
        if (actn(gas) > 0) then
            tl_def = toPasGas(gas)%def
            tl_min = toPasGas(gas)%min
            tl_max = toPasGas(gas)%max
        else
            tl_def = error
            tl_min = error
            tl_max = error
        end if
        write(uto, '(a, i5)') 'Number_of_timelags_used_for_' &
            // trim(gname) // ':', actn(gas)
        if (PwbAggregateSummary) call WritePwbProvenance(uto, gas)
        write(uto, '(a, f9.2)') 'Median_' // trim(gname) // '_timelag_[s]:', tl_def
        write(uto, '(a, f9.2)') 'Mimimum_' // trim(gname) // '_timelag_[s]:', tl_min
        write(uto, '(a, f9.2)') 'Maximum_' // trim(gname) // '_timelag_[s]:', tl_max
        write(uto, '(a)')
    end do

    if (E2Col(h2o)%present .and. ncls > 1) then
        if (PwbAggregateSummary) call WritePwbProvenance(uto, h2o)
        write(uto, '(a, i4)') 'H2O_timelag_determinations_as_a_function_of_relative_humidity'
        write(uto, '(a, i4)') 'Classes with numerosity < 30 are inferred (see software documentation)'
        write(uto,'(a)')             'class     RH-range       med_h2o       min_h2o       max_h2o     class_num'
        do cls = 1, ncls
            write(min, '(i4)') nint((cls - 1) * cls_size)
            call ShrinkString(min)
            write(max, '(i4)') nint(cls * cls_size)
            call ShrinkString(max)
            txt = min(1:len_trim(min)) // ' - ' // max(1:len_trim(max)) // '%'
            write(uto,'(i5, 5x, a9, 3(f13.2,1x), i13)') cls,  txt, toH2O(cls)%def, toH2O(cls)%min, toH2O(cls)%max, h2o_n(cls)
        end do
    end if
    close(uto)
    write(*,'(a)') '  Results written on file: ' &
        // TimelagOpt_Path(1:len_trim(TimelagOpt_Path))
contains

subroutine WritePwbProvenance(unit, gas)
    integer, intent(in) :: unit, gas
    character(32) :: source
    character(32) :: name

    !> Had no default arm, so a gas past the fourth printed whatever `name`
    !> held. One helper, shared with the writer above and with the reader.
    name = TimelagOptGasLabel(gas)
    if (PwbSummarySource(gas) == gas) then
        source = 'native'
    elseif (PwbSummarySource(gas) > 0) then
        select case (PwbSummarySource(gas))
            case (co2)
                source = 'inferred_from_co2'
            case (h2o)
                source = 'inferred_from_h2o'
            case (ch4)
                source = 'inferred_from_ch4'
            case default
                source = 'inferred'
        end select
    else
        source = 'unavailable'
    end if
    write(unit, '(a,a,a,i0,a)') 'PWB_summary_source_for_' // trim(name) // ': ', &
        trim(source), ' (donor evidence=', PwbSummaryEvidence(gas), ')'
end subroutine WritePwbProvenance
end subroutine WriteOutTimelagOptimization

