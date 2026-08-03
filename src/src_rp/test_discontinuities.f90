!***************************************************************************
! test_discontinuities.f90
! ------------------------
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
! \brief       Checks for unphysical discontinuities in the time series \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestDiscontinuities(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: win_len
    integer :: nn = 0
    integer :: wdw_num = 0
    integer :: wdw = 0
    integer :: npoints_par
    integer :: gas
    !> How many variables this routine can flag: the anemometric block plus
    !> the gas block. Not GHGNumVar, which is the *width* of the flag arrays.
    integer :: nvars_tested
    integer :: hflags(GHGNumVar)
    integer :: sflags(GHGNumVar)
    real(kind = dbl) :: Mean(GHGNumVar)
    real(kind = dbl) :: Mean_up(GHGNumVar)
    real(kind = dbl) :: Mean_dw(GHGNumVar)
    real(kind = dbl) :: Var(GHGNumVar)
    real(kind = dbl) :: Var_up(GHGNumVar)
    real(kind = dbl) :: Var_dw(GHGNumVar)
    real(kind = dbl) :: HaarAvr(GHGNumVar)
    real(kind = dbl) :: HaarVar(GHGNumVar)
    real(kind = dbl), allocatable :: XX(:, :)
    real(kind = dbl), allocatable :: XX_dw(:, :)
    real(kind = dbl), allocatable :: XX_up(:, :)
    include '../src_common/interfaces_1.inc'


    write(*, '(a)', advance = 'no') '   Discontinuities test..'

    !> Additional control parameters
    win_len = RPsetup%avrg_len / 6
    if (win_len == 0) win_len = 1

    !> Initializations
    nn = idint((dble(win_len)) * Metadata%ac_freq * 6d1)
    wdw_num = idint(dble(N - nn) / 1d2) + 1

    hflags = 0
    sflags = 0
    !> u, v, w, ts and every gas slot. The early exit below counts against
    !> this rather than against the array width.
    nvars_tested = (ts - u + 1) + (lastGas - firstGas + 1)
    allocate(XX(nn, GHGNumVar))
    allocate(XX_dw(nn/2, GHGNumVar))
    allocate(XX_up(nn/2, GHGNumVar))
    
    do wdw = 1, wdw_num
        npoints_par = 0
        !> Full window
        do i = 1, nn
            XX(i, u:GHGNumVar) = Set(i + 100 * (wdw - 1), u:GHGNumVar)
        end do
        !> Half windows
        do i = 1, nn / 2
            XX_dw(i, :) = XX(i, :)
            XX_up(i, :) = XX(nn/2 + i, :)
        end do
        !> Convert instantaneous molar densities into mole fractions using
        !> standard air molar volume. The 1d3 is the umol basis every gas but
        !> water is held on, so it is asked of the species rather than spelled
        !> out per slot - four arms named four positions and a fifth gas was
        !> left on its raw molar density, which then failed the threshold
        !> comparisons against a mole-fraction limit.
        do gas = firstGas, lastGas
            if (E2Col(gas)%measure_type /= 'molar_density') cycle
            if (GasSlotIsWater(gas)) then
                XX(1:nn, gas) = XX(1:nn, gas) * StdVair
            else
                XX(1:nn, gas) = XX(1:nn, gas) * StdVair * 1d3
            end if
        end do

        !> Whole window mean values
        call AverageNoError(XX, size(XX, 1), size(XX, 2), Mean, error)
        !> Half windows mean values
        call AverageNoError(XX_dw, size(XX_dw, 1), size(XX_dw, 2), Mean_dw, error)
        call AverageNoError(XX_up, size(XX_up, 1), size(XX_up, 2), Mean_up, error)

        !> Whole window variance
        call StDevNoError(XX, size(XX, 1), size(XX, 2), Var, error)
        Var = Var**2
        !> Half windows variances
        call StDevNoError(XX_dw, size(XX_dw, 1), size(XX_dw, 2), Var_dw, error)
        call StDevNoError(XX_up, size(XX_up, 1), size(XX_up, 2), Var_up, error)

        !> Haar functions
        HaarAvr(:) = Mean_dw(:) - Mean_up(:)
        HaarVar(:) = (Var_dw(:) - Var_up(:)) / Var(:)

        !> Hard/soft flags for discontinuities beyond prescribed thresholds
        do j = u, v
            if (HaarAvr(j) > ds%hf_uv)  hflags(j) = 1
            if (HaarAvr(j) > ds%sf_uv)  sflags(j) = 1
            if (HaarVar(j) > ds%hf_var) hflags(j) = 1
            if (HaarVar(j) > ds%sf_var) sflags(j) = 1
        end do
        if (HaarAvr(w) > ds%hf_w)       hflags(w) = 1
        if (HaarAvr(w) > ds%sf_w)       sflags(w) = 1
        if (HaarVar(w) > ds%hf_var)     hflags(w) = 1
        if (HaarVar(w) > ds%sf_var)     sflags(w) = 1
        if (HaarAvr(ts) > ds%hf_t)      hflags(ts) = 1
        if (HaarAvr(ts) > ds%sf_t)      sflags(ts) = 1
        if (HaarVar(ts) > ds%hf_var)    hflags(ts) = 1
        if (HaarVar(ts) > ds%sf_var)    sflags(ts) = 1
        !> Sixteen comparisons over four named slots, all with the same shape
        !> and all reading a per-gas threshold that has been GHGNumVar-wide
        !> since the absolute-limit settings were made per record.
        do gas = firstGas, lastGas
            if (HaarAvr(gas) > ds%hf_gas(gas)) hflags(gas) = 1
            if (HaarAvr(gas) > ds%sf_gas(gas)) sflags(gas) = 1
            if (HaarVar(gas) > ds%hf_var)      hflags(gas) = 1
            if (HaarVar(gas) > ds%sf_var)      sflags(gas) = 1
        end do

        !> Stop early once every variable under test is flagged both ways -
        !> more windows cannot change the outcome.
        !>
        !> The count was against GHGNumVar, which is 68: the width of the flag
        !> arrays, not the number of variables set. Only u, v, w, ts and the
        !> gas slots are ever written, so the sum could not reach it and the
        !> exit never fired - every call ran every window. Counted over the
        !> slots this routine actually flags, it fires again.
        if (count(hflags(u:ts) == 1) + count(hflags(firstGas:lastGas) == 1) &
                == nvars_tested &
            .and. count(sflags(u:ts) == 1) &
                + count(sflags(firstGas:lastGas) == 1) == nvars_tested) exit
    end do
    if(allocated(XX)) deallocate(XX)
    if(allocated(XX_dw)) deallocate(XX_dw)
    if(allocated(XX_up)) deallocate(XX_up)

    ! Pack one digit per variable into the flag strings; absent variables
    ! are marked 9.
    do j = 1, GHGNumVar
        if (.not. E2Col(j)%present) then
            hflags(j) = 9
            sflags(j) = 9
        end if
    end do
    call PackFlagString(hflags, GHGNumVar, CharHF%ds)
    call PackFlagString(sflags, GHGNumVar, CharSF%ds)
    write(*,'(a)') ' Done.'
end subroutine TestDiscontinuities
