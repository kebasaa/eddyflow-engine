!***************************************************************************
! write_out_full.f90
! ------------------
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
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFull(init_string, PeriodRecords, PeriodActualRecords)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    integer, intent(in) :: PeriodRecords
    integer, intent(in) :: PeriodActualRecords
    !> local variables
    integer :: var
    integer :: gas
!    integer :: prof
    character(LongOutstringLen) :: csv_row
    character(DatumLen) :: field_val
    include '../src_common/interfaces.inc'

    !> Preliminary file and timestamp information
    call clearstr(csv_row)
    call AddDatum(csv_row, trim(adjustl(init_string)), separator)
    if (Stats%daytime) then
        call AddDatum(csv_row, '1', separator)
    else
        call AddDatum(csv_row, '0', separator)
    endif
    call WriteDatumInt(PeriodRecords, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(PeriodActualRecords, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> Corrected fluxes (Level 3)
    !> Tau
    call WriteDatumFloat(Flux3%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(QCFlag%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (RUsetup%meth /= 'none') then
        call WriteDatumFloat(Essentials%rand_uncer(u), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> H
    call WriteDatumFloat(Flux3%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(QCFlag%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (RUsetup%meth /= 'none') then
        call WriteDatumFloat(Essentials%rand_uncer(ts), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> LE
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux3%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer_LE, field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Gases
    if(OutVarPresent(co2)) then
        call WriteDatumFloat(Flux3%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(co2), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux3%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(h2o), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(OutVarPresent(ch4)) then
        call WriteDatumFloat(Flux3%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(ch4), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(OutVarPresent(gas4)) then
        call WriteDatumFloat(Flux3%gas4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%gas4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(gas4), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> storage fluxes
    call WriteDatumFloat(Stor%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Stor%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stor%of(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> vertical advection fluxes
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            if (Stats5%Mean(w) /= error .and. Stats%d(gas) >= 0d0) then
                if (gas /= h2o) then
                    call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas) * 1d3, field_val, EddyFlowProj%err_label)
                    call AddDatum(csv_row, field_val, separator)
                else
                    call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas), field_val, EddyFlowProj%err_label)
                    call AddDatum(csv_row, field_val, separator)
                end if
            else
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end if
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> Gas concentrations, densities and timelags
    do gas = co2, gas4
        if (OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%d(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            call WriteDatumFloat(Stats%chi(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            call WriteDatumFloat(Stats%r(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            call WriteDatumFloat(Essentials%used_timelag(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            if (Essentials%def_tlag(gas)) then
                call AddDatum(csv_row, '1', separator)
            else
                call AddDatum(csv_row, '0', separator)
            endif
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, '9', separator)
        end if
    end do

    !> Air properties
    call WriteDatumFloat(Stats7%Mean(ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%Ta, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats%Pr, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(RHO%a, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (RHO%a /= 0d0 .and. RHO%a /= error) then
        call WriteDatumFloat(Ambient%RhoCp / RHO%a, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    call WriteDatumFloat(Ambient%Va, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (Flux3%h2o /= error) then
        call WriteDatumFloat(Flux3%h2o * h2o_to_ET, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    call WriteDatumFloat(RHO%w, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%e, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%es, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%Q, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats%RH, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%VPD, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%Td, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> Unrotated and rotated wind components
    call WriteDatumFloat(Stats4%Mean(u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats4%Mean(v), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats4%Mean(w), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats5%Mean(u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats5%Mean(v), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats5%Mean(w), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%WS, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%MWS, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats4%wind_dir, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> rotation angles
    call WriteDatumFloat(Essentials%yaw, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Essentials%pitch, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Essentials%roll, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> turbulence
    call WriteDatumFloat(Ambient%us, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Stats7%TKE, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%L, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%zL, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%bowen, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(Ambient%Ts, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> footprint
    if (Meth%foot /= 'none') then
        select case(foot_model_used(1:len_trim(foot_model_used)))
            case('kljun_04')
            call AddDatum(csv_row, '0', separator)
            case('kormann_meixner_01')
            call AddDatum(csv_row, '1', separator)
            case('hsieh_00')
            call AddDatum(csv_row, '2', separator)
        end select
        call WriteDatumFloat(Foot%peak, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%offset, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%x10, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%x30, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%x50, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%x70, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(Foot%x90, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, EddyFlowProj%err_label, separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Uncorrected fluxes (Level 0)
    !> Tau
    call WriteDatumFloat(Flux0%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(BPCF%of(w_u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> H
    call WriteDatumFloat(Flux0%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(BPCF%of(w_ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> LE
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux0%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> Gases
    if(OutVarPresent(co2)) then
        call WriteDatumFloat(Flux0%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_co2), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux0%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(OutVarPresent(ch4)) then
        call WriteDatumFloat(Flux0%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_ch4), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(OutVarPresent(gas4)) then
        call WriteDatumFloat(Flux0%gas4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_gas4), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Vickers and Mahrt 97 hard flags
    call AddDatum(csv_row, '8'//CharHF%sr(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%ar(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%do(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%al(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%sk(2:9), separator)
    call AddDatum(csv_row, '8'//CharSF%sk(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%ds(2:9), separator)
    call AddDatum(csv_row, '8'//CharSF%ds(2:9), separator)
    call AddDatum(csv_row, '8'//CharHF%tl(6:9), separator)
    call AddDatum(csv_row, '8'//CharSF%tl(6:9), separator)
    call AddDatum(csv_row, '8'//CharHF%aa(9:9), separator)
    call AddDatum(csv_row, '8'//CharHF%ns(9:9), separator)

    !> Spikes for EddyFlow variables
    call WriteDatumInt(Essentials%e2spikes(u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(Essentials%e2spikes(v), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(Essentials%e2spikes(w), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(Essentials%e2spikes(ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    do var = co2, gas4
        if(OutVarPresent(var)) then
            call WriteDatumInt(Essentials%e2spikes(var), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> LI-COR's diagnostic flags
    if (Diag7200%present) then
        call WriteDatumInt(Diag7200%head_detect, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%t_out, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%t_in, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%aux_in, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%delta_p, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%chopper, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%detector, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%pll, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7200%sync, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (Diag7500%present) then
        call WriteDatumInt(Diag7500%chopper, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7500%detector, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7500%pll, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7500%sync, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (Diag7700%present) then
        call WriteDatumInt(Diag7700%not_ready, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%no_signal, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%re_unlocked, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%bad_temp, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%laser_temp_unregulated, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%block_temp_unregulated, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%motor_spinning, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%pump_on, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%top_heater_on, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%bottom_heater_on, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%calibrating, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%motor_failure, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc1, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc3, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(Diag7700%box_connected, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> AGCs and RSSIs for LI-7200 and LI-7500
    if (Diag7200%present) then
        call WriteDatumInt(nint(Essentials%AGC72), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (Diag7500%present) then
        call WriteDatumInt(nint(Essentials%AGC75), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Variances
    do var = u, ts
        call WriteDatumFloat(Stats%Cov(var, var), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end do
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%Cov(gas, gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do
    !> w-covariances
    call WriteDatumFloat(Stats%Cov(w, ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%Cov(w, gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, &
            trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    enddo

    !> Mean values of user variables
    if (NumUserVar > 0) then
        do var = 1, NumUserVar
            call WriteDatumFloat(UserStats%Mean(var), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end do
    end if

    !> Conditional Eddy Covariance outputs (Zahn et al. 2022)
    if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 2) then
        call WriteDatumFloat(CECFlux%E_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%Tr_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%r_ET_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if
    if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 3) then
        call WriteDatumFloat(CECFlux%Reco_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%GPP_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%NEE_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%r_Fc_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if

    write(uflx, '(a)') csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteOutFull
