!***************************************************************************
! tag_run_mode.f90
! ----------------
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
! \brief       Add token to Timestamp_FilePadding to tag run mode
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TagRunMode()
    use m_common_global_var
    implicit none


    select case(EddyFlowProj%run_mode)
        case('express')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_exp'
        case('advanced')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_adv'
        case('md_retrieval')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_mdr'
        case default
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_nul'
    end select
end subroutine TagRunMode
