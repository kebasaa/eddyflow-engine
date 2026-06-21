!***************************************************************************
! configure_for_md_retrieval.f90
! ------------------------------
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
! \brief       Shut down everything except production of metadata file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ConfigureForMdRetrieval()
    use m_rp_global_var
    implicit none


    !> raw data processing methods
    Meth%tlag = 'none'
    Meth%det = 'ba'
    Meth%rot = 'none'
    Meth%qcflag = 'none'
    Meth%foot = 'none'
    RUsetup%meth = 'none'
    RPsetup%bu_corr = 'none'
    RPsetup%calib_aoa = 'none'
    RPsetup%bu_multi = .false.
    RPsetup%calib_cw = .false.
    RPsetup%filter_by_raw_flags = .false.
    EddyFlowProj%use_extmd_file = .false.
    EddyFlowProj%biomet_data = 'none'
    EddyFlowProj%wpl = .false.
    EddyFlowProj%hf_meth = 'none'

    !> Raw statistical tests
    Test%sr = .false.
    Test%ar = .false.
    Test%do = .false.
    Test%al = .false.
    Test%sk = .false.
    Test%ds = .false.
    Test%tl = .false.
    Test%aa = .false.
    Test%ns = .false.
    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0

    !> Output files and other settings
    EddyFlowProj%out_md         = .true.
    EddyFlowProj%out_fluxnet    = .false.
    EddyFlowProj%out_full       = .false.
    EddyFlowProj%out_avrg_cosp  = .false.
    EddyFlowProj%out_biomet     = .false.
    RPsetup%out_st             = .false.
    RPsetup%filter_sr          = .false.
    RPsetup%filter_al          = .false.
    RPsetup%out_qc_details     = .false.
    RPsetup%out_raw            = .false.
    RPsetup%out_bin_sp         = .false.
    RPsetup%out_bin_og         = .false.
    RPsetup%out_full_sp        = .false.
    RPsetup%out_full_cosp      = .false.
    EddyFlowProj%fcc_follows    = .false.
    EddyFlowProj%make_dataset   = .true.
end subroutine ConfigureForMdRetrieval
