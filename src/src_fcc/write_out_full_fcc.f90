!***************************************************************************
! write_out_full_fcc.f90
! ----------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFullFcc(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(16000) :: csv_row

    !> local variables
    integer :: var
    integer :: i
    integer :: gas
    character(DatumLen) :: field_val
    include '../src_common/interfaces_1.inc'


    call clearstr(csv_row)
    !> Preliminary file and timestamp information
    call AddDatum(csv_row, trim(lEx%fname), separator)
    call AddDatum(csv_row, lEx%end_date(1:10), separator)
    call AddDatum(csv_row, lEx%end_time(1:5), separator)
    call WriteDatumFloat(float_doy, field_val, EddyFlowProj%err_label)
    call stripstr(field_val)  !< Added to fix a strange behaviour
    call AddDatum(csv_row, field_val(1:index(field_val, '.') + 3), separator)
    if (lEx%daytime) then
        call AddDatum(csv_row, '1', separator)
    else
        call AddDatum(csv_row, '0', separator)
    endif
    call WriteDatumInt(lEx%file_records, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(lEx%nr_after_wdf, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> Corrected fluxes (Level 3)
    !> Tau
    call WriteDatumFloat(Flux3%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(QCFlag%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if(RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
        call WriteDatumFloat(lEx%rand_uncer(u), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if

    !> H
    call WriteDatumFloat(Flux3%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(QCFlag%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
        call WriteDatumFloat(lEx%rand_uncer(ts), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if

    !> LE
    if(fcc_var_present(h2o)) then
        call WriteDatumFloat(Flux3%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
            call WriteDatumFloat(lEx%rand_uncer_LE, field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Gases
    if(fcc_var_present(co2)) then
        call WriteDatumFloat(Flux3%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
            call WriteDatumFloat(lEx%rand_uncer(co2), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(fcc_var_present(h2o)) then
        call WriteDatumFloat(Flux3%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
            call WriteDatumFloat(lEx%rand_uncer(h2o), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(fcc_var_present(ch4)) then
        call WriteDatumFloat(Flux3%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
            call WriteDatumFloat(lEx%rand_uncer(ch4), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    if(fcc_var_present(gas4)) then
        call WriteDatumFloat(merge(Flux3%gas4 * gas4_full_flux_sc, error, &
            Flux3%gas4 /= error), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumInt(QCFlag%gas4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        if (RUsetup%meth /= 'none' .or. EddyFlowProj%fix_out_format) then
            call WriteDatumFloat(merge(lEx%rand_uncer(gas4) * gas4_full_flux_sc, error, &
                lEx%rand_uncer(gas4) /= error), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end if
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> storage
    call WriteDatumFloat(lEx%Stor%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if(fcc_var_present(h2o)) then
        call WriteDatumFloat(lEx%Stor%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    do gas = co2, h2o
        if(fcc_var_present(gas)) then
            call WriteDatumFloat(lEx%Stor%of(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do
    do gas = ch4, gas4
        if(fcc_var_present(gas)) then
            if (lEx%Stor%of(gas) /= error) then
                if (gas == gas4) then
                    call WriteDatumFloat(lEx%Stor%of(gas) * 1d-3 * gas4_full_flux_sc, &
                        field_val, EddyFlowProj%err_label)
                else
                    call WriteDatumFloat(lEx%Stor%of(gas) * 1d-3, field_val, EddyFlowProj%err_label)
                end if
                call AddDatum(csv_row, field_val, separator)
            else
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end if 
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> vertical advection fluxes
    do gas = co2, n2o
        if(fcc_var_present(gas)) then
            if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
                if (gas == gas4) then
                    call WriteDatumFloat(lEx%rot_w * lEx%d(gas) * 1d3 * gas4_full_flux_sc, &
                        field_val, EddyFlowProj%err_label)
                    call AddDatum(csv_row, field_val, separator)
                else if (gas /= h2o) then
                    call WriteDatumFloat(lEx%rot_w * lEx%d(gas) * 1d3, field_val, EddyFlowProj%err_label)
                    call AddDatum(csv_row, field_val, separator)
                else
                    call WriteDatumFloat(lEx%rot_w * lEx%d(gas), field_val, EddyFlowProj%err_label)
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
    do gas = co2, n2o
        if (fcc_var_present(gas)) then
            if (gas == gas4) then
                call WriteDatumFloat(merge(lEx%d(gas) * gas4_full_dens_sc, error, &
                    lEx%d(gas) /= error), field_val, EddyFlowProj%err_label)
            else
                call WriteDatumFloat(lEx%d(gas), field_val, EddyFlowProj%err_label)
            end if
            call AddDatum(csv_row, field_val, separator)
            if (gas == gas4) then
                call WriteDatumFloat(merge(lEx%chi(gas) * gas4_full_flux_sc, error, &
                    lEx%chi(gas) /= error), field_val, EddyFlowProj%err_label)
            else
                call WriteDatumFloat(lEx%chi(gas), field_val, EddyFlowProj%err_label)
            end if
            call AddDatum(csv_row, field_val, separator)
            if (gas == gas4) then
                call WriteDatumFloat(merge(lEx%r(gas) * gas4_full_flux_sc, error, &
                    lEx%r(gas) /= error), field_val, EddyFlowProj%err_label)
            else
                call WriteDatumFloat(lEx%r(gas), field_val, EddyFlowProj%err_label)
            end if
            call AddDatum(csv_row, field_val, separator)
            call WriteDatumFloat(lEx%tlag(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
            if (lEx%def_tlag(gas)) then
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
    call WriteDatumFloat(lEx%Ts, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%Ta, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%Pa, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%RHO%a, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (lEx%RHO%a /= 0d0 .and. lEx%RHO%a /= error) then
        call WriteDatumFloat(lEx%RhoCp /lEx%RHO%a, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    call WriteDatumFloat(lEx%Va, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    if (Flux3%h2o /= error) then
        call WriteDatumFloat(Flux3%h2o * 0.0648d0, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    call WriteDatumFloat(lEx%RHO%w, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%e, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%es, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%Q, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%RH, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%VPD, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%Tdew, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> Unrotated and rotated wind components
    call WriteDatumFloat(lEx%unrot_u, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%unrot_v, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%unrot_w, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%rot_u, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%rot_v, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%rot_w, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%WS, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%MWS, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%WD, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> rotation angles
    call WriteDatumFloat(lEx%yaw, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%pitch, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%roll, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)

    !> turbulence
    call WriteDatumFloat(Flux3%ustar, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%TKE, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%L, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%zL, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%bowen, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(lEx%Tstar, field_val, EddyFlowProj%err_label)
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
    call WriteDatumFloat(lEx%Flux0%tau, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(BPCF%of(w_u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> H
    call WriteDatumFloat(lEx%Flux0%H, field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumFloat(BPCF%of(w_ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    !> LE
    if(fcc_var_present(h2o)) then
        call WriteDatumFloat(lEx%Flux0%LE, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> Gases
    if(fcc_var_present(co2)) then
        call WriteDatumFloat(lEx%Flux0%co2, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_co2), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(fcc_var_present(h2o)) then
        call WriteDatumFloat(lEx%Flux0%h2o, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(fcc_var_present(ch4)) then
        call WriteDatumFloat(lEx%Flux0%ch4, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_ch4), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if(fcc_var_present(gas4)) then
        call WriteDatumFloat(merge(lEx%Flux0%gas4 * gas4_full_flux_sc, error, &
            lEx%Flux0%gas4 /= error), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(BPCF%of(w_gas4), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Vickers and Mahrt 97 flags
    do i = 1, 8
        write(field_val, *) lEx%vm_flags(i)
        call AddDatum(csv_row, field_val, separator)
    end do
    call AddDatum(csv_row, lEx%vm_tlag_hf, separator)
    call AddDatum(csv_row, lEx%vm_tlag_sf, separator)
    call AddDatum(csv_row, lEx%vm_aoa_hf, separator)
    call AddDatum(csv_row, lEx%vm_nshw_hf, separator)

    !> Spikes for EddyFlow variables
    call WriteDatumInt(lEx%spikes(u), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(lEx%spikes(v), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(lEx%spikes(w), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    call WriteDatumInt(lEx%spikes(ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    do var = co2, gas4
        if(fcc_var_present(var)) then
            call WriteDatumInt(lEx%spikes(var), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> LI-COR's diagnostic flags
    if (Diag7200%present) then
        do i = 1, 9
            call WriteDatumInt(nint(lEx%licor_flags(i)), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end do
    elseif(EddyFlowProj%fix_out_format) then
        do i = 1, 9
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if
    if (Diag7500%present) then
        do i = 10, 13
            call WriteDatumInt(nint(lEx%licor_flags(i)), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end do
    elseif(EddyFlowProj%fix_out_format) then
        do i = 1, 4
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if
    if (Diag7700%present) then
        do i = 14, 29
            call WriteDatumInt(nint(lEx%licor_flags(i)), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end do
    elseif(EddyFlowProj%fix_out_format) then
        do i = 1, 16
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if

    !> AGCs
    if (Diag7200%present) then
        call WriteDatumInt(nint(lEx%agc72), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (Diag7500%present) then
        call WriteDatumInt(nint(lEx%agc75), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    elseif(EddyFlowProj%fix_out_format) then
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
!        if (Diag7700%present) then
!            call WriteDatumInt(nint(lEx%rssi77), field_val, EddyFlowProj%err_label)
!            call AddDatum(csv_row, field_val, separator)
!        elseif(EddyFlowProj%fix_out_format) then
!            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
!        end if

    !> Variances
    do var = u, ts
        call WriteDatumFloat(lEx%var(var), field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end do
    do gas = co2, gas4
        if(fcc_var_present(gas)) then
            call WriteDatumFloat(lEx%var(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do
    !> w-covariances
    call WriteDatumFloat(lEx%cov_w(ts), field_val, EddyFlowProj%err_label)
    call AddDatum(csv_row, field_val, separator)
    do gas = co2, gas4
        if(fcc_var_present(gas)) then
            call WriteDatumFloat(lEx%cov_w(gas), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        elseif(EddyFlowProj%fix_out_format) then
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    enddo

    !> Mean values of user variables
    if (lEx%ncustom > 0) then
        do var = 1, lEx%ncustom
            call WriteDatumFloat(lEx%user_var(var), field_val, EddyFlowProj%err_label)
            call AddDatum(csv_row, field_val, separator)
        end do
    end if

    !> Conditional Eddy Covariance outputs (Zahn et al. 2022)
    if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 2) then
        call WriteDatumFloat(CECFlux%E_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%Tr_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%E_cec_ET, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%Tr_cec_ET, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%r_ET_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if
    if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 3) then
        call WriteDatumFloat(CECFlux%Reco_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%P_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%NEE_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
        call WriteDatumFloat(CECFlux%r_Fc_cec, field_val, EddyFlowProj%err_label)
        call AddDatum(csv_row, field_val, separator)
    end if

    write(uflx, '(a)')   csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteOutFullFcc
