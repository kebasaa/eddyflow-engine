!***************************************************************************
! configure_project_for_express.f90
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
! \brief       The part of the Express configuration that both RP and FCC
!              read: the [Project] settings that decide what the final
!              fluxes are and which files carry them.
! \author      Gerardo Fratini
! \note        Express used to be applied by RP alone (ConfigureForExpress),
!              while FCC read the same project file and took its values as
!              written. RP then ran with WPL on and FCC computed the final
!              fluxes with whatever wpl_meth the file held - so an express
!              project carrying wpl_meth=0 had RP convert to mixing ratios
!              point by point and FCC leave the fluxes uncorrected. LE came
!              out ~22x off EddyPro's, which ignores the 0 in express mode.
!              Called from WriteProcessingProjectVariables, so it takes
!              effect in both programs before either derives anything from
!              these settings (FCC's spectral-assessment switches in
!              WriteVariablesFCC are computed from hf_meth and out_avrg_*).
!              RP's ConfigureForExpress calls it again, which is harmless:
!              it only assigns constants.
!              The rest of Express - raw-data methods, tests, despiking,
!              Burba - lives in RP's types and is RP's alone; FCC takes
!              their effect from the ex file. So do out_md and out_biomet:
!              FCC never reads out_biomet, and nothing in the metadata file
!              depends on what FCC computes, while FCC's metadata writer
!              lays it out differently from RP's - forcing it on here would
!              have FCC overwrite the file RP already wrote.
! \sa          ConfigureForExpress
! \bug
! \deprecated
! \test        static_checks/test_express_shared_static.py
! \todo
!***************************************************************************
subroutine ConfigureProjectForExpress
    use m_common_global_var
    implicit none


    !> Flux computation and corrections
    EddyFlowProj%hf_meth = 'moncrieff_97'
    EddyFlowProj%hf_meth_in_situ = .false.
    !> Express is a fixed configuration, so it states the cospectral
    !> shape rather than inheriting whatever the project asked for.
    EddyFlowProj%cosp_model = 'moncrieff_97'
    EddyFlowProj%wpl  = .true.
    Meth%qcflag = 'mauder_foken_04'
    Meth%foot = 'kljun_04'
    RUsetup%meth = 'none'

    !> Output files that carry, or decide, the fluxes FCC computes
    EddyFlowProj%out_full     = .true.
    EddyFlowProj%out_avrg_cosp = .false.
    EddyFlowProj%out_avrg_spec = .false.
    EddyFlowProj%make_dataset = .true.
end subroutine ConfigureProjectForExpress
