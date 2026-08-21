!***************************************************************************
! test_timelag.f90
! ----------------
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
! \brief       Check for scalar time-lags out of suggested ranges    \n
!              hard-flags and soft-flags file accordingly.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestTimeLag(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: rlag(GHGNumVar)
    !> One entry per variable, indexed by slot, as every other test's flag
    !> arrays are. They were four wide because the packing below was a
    !> four-digit integer, and that width is what bounded the test itself.
    integer :: hflags(GHGNumVar)
    integer :: sflags(GHGNumVar)
    integer :: min_rl(GHGNumVar)
    integer :: max_rl(GHGNumVar)
    integer :: def_rl(GHGNumVar)
    real(kind = dbl) :: FirstCol(N)
    real(kind = dbl) :: SecondCol(N)
    real(kind = dbl) :: tlag(GHGNumVar) = 0.d0
    real(kind = dbl) :: DefCov(GHGNumVar)
    real(kind = dbl) :: MaxCov(GHGNumVar)

    call LogSayNoAdv('   Time lag test..')

    !> Initializations
    hflags = 0
    sflags = 0

    !> Define min e max "row-lags" for scalars, using timelags retrieved from metadata file
    do j = firstGas, GHGNumVar
        min_rl(j) = nint(E2Col(j)%min_tl * Metadata%ac_freq)
        max_rl(j) = nint(E2Col(j)%max_tl * Metadata%ac_freq)
    end do
    !> Default values are taken from EddyFlow settings. tl%def_gas has always
    !> been GHGNumVar wide; only these four assignments were not.
    do j = firstGas, lastGas
        def_rl(j) = nint(tl%def_gas(j) * Metadata%ac_freq)
    end do

    !> Actual time-lags (tlag), maximum of the cov. (Rmax) \n
    !>  and cov. for default timelag (R0)
    !> Flags if the difference is too high
    !>
    !> Every configured gas. This ran co2..gas4 because hflags/sflags were four
    !> wide, which was in turn because the packing below was a four-digit
    !> integer - so the whole test, not just its reporting, stopped at the
    !> fourth gas. A fifth gas's lag was never compared against its default and
    !> its flag column read as "not performed" while the test was enabled.
    do i = firstGas, lastGas
        if (i - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (E2Col(i)%present) then
            FirstCol(:)  = Set(:, w)
            SecondCol(:) = Set(:, i)
            call CovMaxRS(def_rl(i), min_rl(i), max_rl(i), FirstCol, SecondCol &
                , MaxCov(i), DefCov(i), tlag(i), rlag(i), N)
            if((MaxCov(i) - DefCov(i)) * 1d2 / DefCov(i) >= tl%hf_lim) hflags(i) = 1
            if((MaxCov(i) - DefCov(i)) * 1d2 / DefCov(i) >= tl%sf_lim) sflags(i) = 1
        end if
    end do

    !> One digit per variable, as the other eight tests pack theirs.
    !>
    !> This was 90000 + sum(hflags(j) * 10**(4 - j)) rendered through int2char,
    !> the last base-10 packing in the engine. It caps the variable count at
    !> about nine before a 32-bit integer overflows, which is why the flag
    !> arrays - and therefore the test - could not grow.
    !>
    !> The emitted cell is unchanged at four gas records. int2char right-aligns
    !> the five-digit number in a zero-padded FlagStrLen field, so slots five
    !> to eight landed at positions 66 to 69; PackFlagString writes slot j at
    !> position j + 1, so they land at 6 to 9. Same four digits, same order -
    !> the writers slice from firstGas + 1 rather than from the end.
    call PackFlagString(hflags, GHGNumVar, CharHF%tl)
    call PackFlagString(sflags, GHGNumVar, CharSF%tl)
    call LogSay(' Done.')
end subroutine TestTimeLag

!***************************************************************************
!
! \brief       Performs covariance analysis for determinig the "optimal" time-lag
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovMaxRS(lagctr, lagmin, lagmax, Col1, Col2, MaxCov, DefCov, TLag, RLag, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N !< Number of raw data
    integer, intent(in) :: lagmin
    integer, intent(in) :: lagmax
    integer, intent(in) :: lagctr
    real(kind = dbl), intent(in) :: Col1(N) !< Fluctuations
    real(kind = dbl), intent(in) :: Col2(N) !< Fluctuations
    integer, intent(out) :: RLag
    real(kind = dbl), intent(out) :: TLag
    real(kind = dbl), intent(out) :: DefCov
    real(kind = dbl), intent(out) :: MaxCov
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: ii = 0
    integer :: N2 !< Number of rows of the raw dataset after time-lag compensation
    real(kind = dbl), allocatable :: ShLocSet(:, :)
    real(kind = dbl) :: Mean(2)
    real(kind = dbl) :: Cov

    Cov = 0.d0
    MaxCov = 0.d0
    DefCov = 0.d0
    TLag = 0.d0
    do i = lagmin, lagmax
        N2 = N - abs(i)
        allocate(ShLocSet(N2, 2))

        do ii = 1, N2
            if (i < 0) then
                ShLocSet(ii, 1) = Col1(ii - i)
                ShLocSet(ii, 2) = Col2(ii)
            else
                ShLocSet(ii, 1) = Col1(ii)
                ShLocSet(ii, 2) = Col2(ii + i)
            end if
        end do
        Mean = sum(ShLocSet, dim = 1)
        Mean = Mean / dble(N2)
        ShLocSet(:, 1) = ShLocSet(:, 1) - Mean(1)
        ShLocSet(:, 2) = ShLocSet(:, 2) - Mean(2)
        !> covariance
        do j = 1, N2
            Cov = Cov + ShLocSet(j, 1) * ShLocSet(j, 2)
        end do
        Cov = Cov / dble(N2 - 1)
        if (i == lagctr) DefCov = Cov
        !> max cov and actual time-lag
        if (abs(Cov) > abs(MaxCov)) then
            MaxCov = Cov
            TLag = dble(i) / Metadata%ac_freq
            RLag = i
        end if
        deallocate(ShLocSet)
    end do
end subroutine CovMaxRS
