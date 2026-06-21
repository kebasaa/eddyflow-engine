!***************************************************************************
! configure_for_express.f90
! -------------------------
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
subroutine ConfigureForFluxnet
    use m_common_global_var
    implicit none

    !> Error code must be -9999
    EddyFlowProj%err_label = '-9999'
    

end subroutine ConfigureForFluxnet
