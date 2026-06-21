!***************************************************************************
! override_settings.f90
! ---------------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Forces some operations (regardless of user choice) based on instrument
!              models (e.g. CSAT3 no cross-wind correction) and logic
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OverrideSettings()
    use m_rp_global_var
    implicit none

    !> If biomet measurements are not to be used, they are also not to be output
    if (EddyFlowProj%biomet_data == 'none') EddyFlowProj%out_biomet = .false.

    !> if there is no LI-7500 among the instruments, Burba terms should not be calculated
    if (index(E2Col(co2)%Instr%model, 'li7500') == 0 &
        .and. index(E2Col(h2o)%Instr%model,'li7500') == 0) &
        RPsetup%bu_corr = 'none'
end subroutine OverrideSettings
