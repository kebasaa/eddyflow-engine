!***************************************************************************
! configure_for_express.f90
! -------------------------
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
! \brief       Set all entries to predefined values, valid for Express
!              This bypasses all user settings in terms of processing choices
! \author      Gerardo Fratini
! \note        Angle of attack correction is set later in main,
!              because it requires master_sonic
!              which is not know at this stage, when running in embedded mode
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForExpress
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    character(32) :: species
    include '../src_common/interfaces.inc'


    !> raw data processing methods
    Meth%tlag = 'maxcov&default'
    Meth%det = 'ba'
    Meth%rot = 'double_rotation'
    Meth%qcflag = 'mauder_foken_04'
    Meth%foot = 'kljun_04'
    EddyFlowProj%hf_meth = 'moncrieff_97'
    EddyFlowProj%wpl  = .true.
    RUsetup%meth = 'none'
    RPsetup%bu_corr = 'none'
    RPsetup%bu_multi = .false.
    RPsetup%filter_sr = .true.
    RPsetup%filter_al = .true.
    RPsetup%calib_cw = .false.
    RPSetup%despike_vickers97 = .true.
    RPsetup%calib_aoa = 'automatic'

    !> Raw statistical tests
    Test%sr = .true.
    Test%ar = .true.
    Test%do = .true.
    Test%al = .true.
    Test%sk = .true.
    Test%ds = .false.
    Test%tl = .false.
    Test%aa = .false.
    Test%ns = .false.

    !> test parameters
    sr%num_spk = 3
    sr%lim_u   = 3.5d0
    sr%hf_lim  = 1d0
    sr%lim_w   = 5.0d0
    !> Per-gas defaults, chosen by species rather than by slot. Written as
    !> co2/h2o/ch4/gas4 they described positions: a fifth gas got whichever
    !> value the blanket assignment left, and a project ordering its records
    !> differently got CO2's plausibility band applied to something else.
    !>
    !> The fall-through is the fourth slot's old set - a wide band and the
    !> looser spike limit - which is exactly what "any other trace gas" meant
    !> when it was spelled gas4.
    sr%lim_gas = 8d0
    do gas = firstGas, lastGas
        species = GasOutputLabel(gas)
        call lowercase(species)
        select case (trim(adjustl(species)))
            case ('co2', 'h2o'); sr%lim_gas(gas) = 3.5d0
        end select
    end do
    ar%lim     = 7
    ar%bins    = 100
    ar%hf_lim  = 70
    do%extlim_dw = 10d0
    do%hf1_lim = 10d0
    do%hf2_lim = 6d0
    al%u_max   = 30d0
    al%w_max   = 5d0
    al%t_min   = -40d0
    al%t_max   = 50d0
    do gas = firstGas, lastGas
        species = GasOutputLabel(gas)
        call lowercase(species)
        select case (trim(adjustl(species)))
            case ('co2')
                al%gas_min(gas) = 200d0
                al%gas_max(gas) = 900d0
            case ('h2o')
                al%gas_min(gas) = 0d0
                al%gas_max(gas) = 40d0
            case ('ch4')
                !< 1.7 ppm is minimum in unpolluted troposphere, 0.1 is safety factor
                al%gas_min(gas) = 1.7d0 * 0.1d0
                al%gas_max(gas) = 1000d0    !< to be better assessed
            case default
                !< no default lower bound - any trace gas at any concentration
                al%gas_min(gas) = 0d0
                al%gas_max(gas) = 1000d0    !< to be better assessed
        end select
    end do
    sk%hf_skmin = -2d0
    sk%hf_skmax = 2d0
    sk%sf_skmin = -1d0
    sk%sf_skmax = 1d0
    sk%hf_kumin = 1d0
    sk%hf_kumax = 8d0
    sk%sf_kumin = 2d0
    sk%sf_kumax = 5d0
    BurbaPar%l(daytime, bot, 1)     =  0.944d0
    BurbaPar%l(daytime, bot, 2)     =  2.57d0
    BurbaPar%l(daytime, top, 1)     =  1.005d0
    BurbaPar%l(daytime, top, 2)     =  0.24d0
    BurbaPar%l(daytime, spar, 1)    =  1.01d0
    BurbaPar%l(daytime, spar, 2)    =  0.36d0
    BurbaPar%l(nighttime, bot, 1)   =  0.883d0
    BurbaPar%l(nighttime, bot, 2)   =  2.17d0
    BurbaPar%l(nighttime, top, 1)   =  1.008d0
    BurbaPar%l(nighttime, top, 2)   =  -0.41d0
    BurbaPar%l(nighttime, spar, 1)  =  1.01d0
    BurbaPar%l(nighttime, spar, 2)  =  -0.17d0

    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0

    !> Output files and other settings
    ! EddyFlowProj%out_fluxnet  = .false.
    EddyFlowProj%out_full     = .true.
    EddyFlowProj%out_md       = .true.
    RPsetup%out_st        = .true.
    RPsetup%out_qc_details = .false.
    RPsetup%out_raw        = .false.
    RPsetup%out_bin_sp     = .false.
    RPsetup%out_bin_og     = .false.
    RPsetup%out_full_sp    = .false.
    RPsetup%out_full_cosp  = .false.
    EddyFlowProj%out_avrg_cosp = .false.
    EddyFlowProj%out_avrg_spec = .false.
    EddyFlowProj%fcc_follows  = .false.
    EddyFlowProj%make_dataset = .true.

    if (EddyFlowProj%biomet_data /= 'none') then
        EddyFlowProj%out_biomet = .true.
    else
        EddyFlowProj%out_biomet = .false.
    end if

end subroutine ConfigureForExpress
