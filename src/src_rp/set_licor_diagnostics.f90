!***************************************************************************
! set_licor_diagnostics.f90
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
! \brief       Set AGC and RSSI as available
! \author      Gerardo Fratini
! \note        AGC72 and AGC75 each hold either an AGC or an RSSI, because no
!              analyser provides both. Which of the two it is is what
!              CecSignalIsRssi decides, and it matters only where the number is
!              compared against a threshold.
! \note        Every branch of this loop used to `exit`, not only the one that
!              matched its own analyser, so the first signal-strength column of
!              any kind ended the search: a site with an LI-7500 and an LI-7700
!              filled one slot and left the other at its error value. Each
!              analyser now takes the first column that is its own, and the
!              loop runs to the end.
! \sa          CecSignalColumnFor, which resolves the same columns per gas
!***************************************************************************
subroutine SetLicorDiagnostics(M)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: M
    !> local variables
    integer :: i
    character(32) :: model

    Essentials%AGC72 = error
    Essentials%AGC75 = error
    Essentials%RSSI77 = error

    !> First, the columns the project declares as signal strength.
    do i = 1, M
        if (.not. IsSignalColumn(i)) cycle
        model = UserCol(i)%instr%model
        if (index(model, 'li7200') /= 0 .and. Essentials%AGC72 == error) &
            Essentials%AGC72 = UserStats%Mean(i)
        if (index(model, 'li7500') /= 0 .and. Essentials%AGC75 == error) &
            Essentials%AGC75 = UserStats%Mean(i)
        if (index(model, 'li7700') /= 0 .and. Essentials%RSSI77 == error) &
            Essentials%RSSI77 = UserStats%Mean(i)
    end do
    !> Then, if AGC from LI-COR's flags are available, override any previous value with it
    !> This is because the user may call AGC a column that is not the actual AGC, while the
    !> value from LI-COR's flags are surely correct
    if (Diag7200%AGC  /= 0d0) Essentials%AGC72 = Diag7200%AGC
    if (Diag7500%AGC /= 0d0) Essentials%AGC75 = Diag7500%AGC

contains

!> Whether user column \a i is a signal strength.
!>
!> The project's records answer first: they name the column, so the answer does
!> not turn on how the metadata spelled the variable. Falling back to the name
!> is what a project written before the records existed has, and it is what the
!> whole of this used to do - case-sensitively, so `agc` from another tool was
!> not a signal strength and nothing said so.
logical function IsSignalColumn(i)
    integer, intent(in) :: i
    integer :: k
    character(32) :: name

    IsSignalColumn = .false.
    if (i < 1 .or. i > NumUserVar) return

    do k = 1, min(EddyFlowProj%agc_num, MaxNumAgcCols)
        if (EddyFlowProj%agc(k)%col /= UserCol(i)%orig_col) cycle
        IsSignalColumn = .true.
        return
    end do

    name = UserCol(i)%var
    call lowercase(name)
    IsSignalColumn = trim(adjustl(name)) == 'agc' &
        .or. trim(adjustl(name)) == 'rssi'
end function IsSignalColumn
end subroutine SetLicorDiagnostics
