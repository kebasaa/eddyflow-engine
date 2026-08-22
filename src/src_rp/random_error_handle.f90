!***************************************************************************
! random_error_handle.f90
! -----------------------
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
! \brief       Estimate flux random uncertainty according to the selected method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RandomUncertaintyHandle(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)


    call LogSay('  Estimating random uncertainty..')

    !> Calculate random uncertainty
    select case (RUsetup%meth)
        case('finkelstein_sims_01')
            call IntegralTurbulenceScale(Set, size(Set, 1), size(Set, 2))
            call RU_Finkelstein_Sims_01(Set, nrow, ncol)
        case('mann_lenschow_94')
            call IntegralTurbulenceScale(Set, size(Set, 1), size(Set, 2))
            call RU_Mann_Lenschow_04(nrow)
        case('lenschow_00')
            !> No integral turbulence scale either: the fit is over a fixed
            !> five lags, not over a scale this would measure.
            call RU_Lenschow_00(Set, nrow, ncol)
        case('billesbach_11')
            !> No integral turbulence scale: the shuffle destroys the
            !> autocorrelation this would measure, and the estimate does not
            !> use it. The other two need it because they integrate over lags.
            call RU_Billesbach_11(Set, nrow, ncol)
        case('none')
            Essentials%rand_uncer(u:lastGas) = error
            Essentials%rand_uncer_LE = error
            Essentials%rand_uncer_ET = error
        case('mahrt_98')
            !> Mahrt has been calculated already, so don't need to do anything
            continue
        case default
            call ExceptionHandler(42)
            Essentials%rand_uncer(u:lastGas) = error
            Essentials%rand_uncer_LE = error
            Essentials%rand_uncer_ET = error
            return
    end select
    call LogSay('  Done.')
end subroutine RandomUncertaintyHandle

!***************************************************************************
!
! \brief       Estimate random error according to \n
!              Finkelstein and Sims (2001), Eq. 8- 10
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RU_Finkelstein_Sims_01(Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    !> local variables
    integer :: var
    integer :: lag
    integer :: LagMax(M)
    integer :: errcnt
    real(kind = dbl), allocatable :: gam(:, :, :)
    real(kind = dbl) :: varcov
    real(kind = dbl), external :: LaggedCovarianceNoError

    !> Define max lag based on ITS
    LagMax(u:lastGas) = nint(ITS(u:lastGas) * Metadata%ac_freq)
    where (LagMax < 0) LagMax = nint(error)
    do var = u, lastGas
        if (var == v .or. var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error .and. LagMax(var) /= nint(error)) then
            allocate (gam(0:LagMax(var), 2, 2))
            gam = 0d0
            do lag = 0, LagMax(var)
                gam(lag, 1, 1) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, w), &
                                            size(Set, 1), lag, error)
                gam(lag, 2, 2) = &
                    LaggedCovarianceNoError(Set(:, var), Set(:, var), &
                                            size(Set, 1), lag, error)
                gam(lag, 1, 2) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, var), &
                                            size(Set, 1), lag, error)
                gam(lag, 2, 1) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, var), &
                                            size(Set, 1), -lag, error)
            end do

            !> variance of covariances, Eq. 8  in Finkelstein & Sims (2001, JGR)
            !> Initialize the value for lag = 0
            varcov = 0d0
            if (gam(0, 1, 1) /= error .and. gam(0, 2, 2) /= error) &
                varcov = gam(0, 1, 1) * gam(0, 2, 2) + gam(0, 1, 2) * gam(0, 2, 1)

            !> Now cycle on lag. Do it one sided and multiply by 2 (Eq. 9 and 10)
            errcnt = 0
            do lag = 1, LagMax(var)
                if (gam(lag, 1, 1) /= error .and. gam(0, 2, 2) /= error &
                    .and. gam(0, 1, 2) /= error .and. gam(0, 2, 1) /= error) then
                    varcov = varcov + 2d0 * gam(lag, 1, 1) * gam(lag, 2, 2) &
                        + 2d0 * gam(lag, 1, 2) * gam(lag, 2, 1)
                else
                    errcnt = errcnt + 1
                end if
            enddo
            deallocate (gam)
            !> Normalization (see Eq. 8 for why using N here)
            varcov = varcov / dfloat(N - errcnt)

            !> Random error is the square root of this variance
            if (varcov /= 0) then
                Essentials%rand_uncer(var) = dsqrt(abs(varcov))
            else
                Essentials%rand_uncer(var) = error
            end if
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RU_Finkelstein_Sims_01

!***************************************************************************
!
! \brief       Estimate random error according to Mann and Lenschow (1994)
!              See e.g. Eq. 5 in Finkelstein and Sims (2001)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RU_Mann_Lenschow_04(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: var
    real(kind = dbl) :: corr_coeff(E2NumVar)

    do var = u, lastGas
        if (var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error) then

            !> Correlation coefficient
            corr_coeff(var) = dabs(Stats%cov(w, var)) &
                / (dsqrt(Stats%cov(w, w)) * dsqrt(Stats%cov(var, var)))

            !> Random uncertainty
            Essentials%rand_uncer(var) = abs(Stats%cov(w, var)) &
                * dsqrt((1d0 + corr_coeff(var)**2) / corr_coeff(var)**2) &
                * dsqrt (2d0 * ITS(var) / (N / Metadata%ac_freq))
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RU_Mann_Lenschow_04


!***************************************************************************
!
! \brief       Estimate random error according to \n
!              Mahrt (1998), Eqs. 8 - 9
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo

!***************************************************************************
!
! \brief       Random uncertainty by Mahrt (1998) 6x6 sub-sampling
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Mahrt, L. (1998). Flux sampling errors for aircraft and towers.
!            Boundary-Layer Meteorol. 88: 163-187. Eqs. (8)-(10).
subroutine RU_Mahrt_98(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer, parameter :: n_sub    = 6
    integer, parameter :: n_subsub = 6
    integer :: sub_idx, subsub_idx, gas_var
    integer :: sub_len, subsub_len
    real(kind = dbl) :: cov_mat(GHGNumVar, GHGNumVar)
    real(kind = dbl) :: subsub_cov(n_subsub, GHGNumVar)
    real(kind = dbl) :: all_cov(n_sub*n_subsub, GHGNumVar)
    real(kind = dbl) :: sub_mean(GHGNumVar), sub_means(n_sub, GHGNumVar)
    real(kind = dbl) :: grand_mean(GHGNumVar), between_ss(GHGNumVar)
    real(kind = dbl) :: sigma_wi(n_sub, GHGNumVar), sigma_btw(GHGNumVar)
    real(kind = dbl), allocatable :: sub_chunk(:,:), ss_chunk(:,:)

    sub_len    = nrow / n_sub
    subsub_len = sub_len / n_subsub

    allocate(sub_chunk(sub_len, GHGNumVar))
    do sub_idx = 1, n_sub
        sub_chunk = Set(sub_len*(sub_idx-1)+1 : sub_len*sub_idx, 1:GHGNumVar)
        allocate(ss_chunk(subsub_len, GHGNumVar))
        do subsub_idx = 1, n_subsub
            ss_chunk = sub_chunk( &
                subsub_len*(subsub_idx-1)+1 : subsub_len*subsub_idx, :)
            call CovarianceMatrixNoError(ss_chunk, subsub_len, GHGNumVar, cov_mat, error)
            subsub_cov(subsub_idx, :) = cov_mat(w, :)
            all_cov(n_subsub*(sub_idx-1)+subsub_idx, :) = cov_mat(w, :)
        end do
        deallocate(ss_chunk)
        ! Sub-period mean covariance (F_i_bar, Mahrt 1998 Eq. 8)
        call AverageNoError(subsub_cov, n_subsub, GHGNumVar, sub_mean, error)
        sub_means(sub_idx, :) = sub_mean
        sigma_wi(sub_idx, :) = 0d0
        do subsub_idx = 1, n_subsub
            where (subsub_cov(subsub_idx,:) /= error .and. sub_mean /= error) &
                sigma_wi(sub_idx,:) = sigma_wi(sub_idx,:) + &
                    (subsub_cov(subsub_idx,:) - sub_mean)**2
        end do
        where (sigma_wi(sub_idx,:) > 0d0) &
            sigma_wi(sub_idx,:) = dsqrt(sigma_wi(sub_idx,:) / dble(n_subsub-1))
    end do
    deallocate(sub_chunk)

    ! Grand mean (F_bar, Mahrt 1998 Eq. 10)
    call AverageNoError(all_cov, n_sub*n_subsub, GHGNumVar, grand_mean, error)

    ! RE = mean of within-period sigmas / sqrt(n_subsub) (Mahrt 1998 Eq. 8)
    do gas_var = u, lastGas
        if (E2Col(gas_var)%present) then
            Essentials%rand_uncer(gas_var) = &
                sum(sigma_wi(:, gas_var)) / n_sub / dsqrt(dble(n_subsub))
        else
            Essentials%rand_uncer(gas_var) = error
        end if
    end do

    ! Between-sub-period sigma (sigma_btw, Mahrt 1998 Eq. 10)
    between_ss = 0d0
    do sub_idx = 1, n_sub
        where (sub_means(sub_idx,:) /= error .and. grand_mean /= error) &
            between_ss = between_ss + (sub_means(sub_idx,:) - grand_mean)**2
    end do
    sigma_btw = 0d0
    where (between_ss > 0d0) sigma_btw = dsqrt(between_ss / dble(n_sub-1))

    do gas_var = u, lastGas
        if (E2Col(gas_var)%present .and. Essentials%rand_uncer(gas_var) /= error &
            .and. Essentials%rand_uncer(gas_var) > 0d0) then
            Essentials%mahrt98_NR(gas_var) = &
                sigma_btw(gas_var) / Essentials%rand_uncer(gas_var)
        else
            Essentials%mahrt98_NR(gas_var) = error
        end if
    end do
end subroutine RU_Mahrt_98
!***************************************************************************
!
! \brief       Instrumental noise from the autocovariance intercept, after
!              Lenschow et al. (2000) as applied by Mauder et al. (2013).
! \author      Jonathan Muller
! \note
!              WHAT IT MEASURES. White instrument noise is uncorrelated
!              between samples, so it lands entirely in the lag-zero
!              autocovariance and nowhere else. The atmospheric part is
!              continuous across lag zero. Fit a line to the first few lags,
!              extrapolate it back to zero, and the gap between that line and
!              the measured lag-zero value IS the noise variance.
!
!              NOT THE SAME AS THE OTHER FOUR. Finkelstein & Sims and Mann &
!              Lenschow estimate a SAMPLING error - how much this half hour's
!              covariance would differ from another draw of the same process.
!              Billesbach estimates a floor below which a flux cannot be told
!              from nothing. This estimates the INSTRUMENT's own noise, and is
!              the smallest of the four: it says what the analyser contributes
!              and nothing about whether the atmosphere was sampled long
!              enough. They are not interchangeable.
!
!              LAGS 1 TO 5, IN SAMPLES. EddyUH's own window
!              (EC_Software_Preproc/EddyUH_unc_Preproc.m:157-160), fixed
!              rather than derived from the sampling rate, so at 10 Hz it
!              spans 0.1 to 0.5 s and at 20 Hz half that. Reproduced as it
!              stands, because the number of lags is what the linear
!              extrapolation is conditioned on; deriving it from the rate
!              would be a different method wearing the same citation.
!
!              THE VERTICAL WIND GATES IT. The noise variance of w is
!              computed and checked but does not enter the result. EddyUH
!              rejects the period when either intercept comes out
!              non-positive, and so does this: a w autocovariance with no
!              noise-like step at lag zero means the assumption the whole
!              method rests on does not hold here, whatever the scalar did.
! \sa          RU_Finkelstein_Sims_01, RU_Billesbach_11
!***************************************************************************
subroutine RU_Lenschow_00(Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    !> local variables
    !> The fit window, in SAMPLES. EddyUH's, and fixed by the method rather
    !> than by the acquisition rate - see the note above.
    integer, parameter :: first_lag = 1
    integer, parameter :: last_lag = 5
    integer :: var
    integer :: lag
    real(kind = dbl) :: acov_w(0:last_lag)
    real(kind = dbl) :: acov_g(0:last_lag)
    real(kind = dbl) :: w_noise
    real(kind = dbl) :: g_noise
    logical :: usable
    real(kind = dbl), external :: LaggedCovarianceNoError

    !> The vertical wind's autocovariance is the same for every gas, so it is
    !> measured once rather than once per column.
    usable = .true.
    do lag = 0, last_lag
        acov_w(lag) = LaggedCovarianceNoError(Set(:, w), Set(:, w), N, lag, error)
        if (acov_w(lag) == error) usable = .false.
    end do
    if (usable) then
        w_noise = acov_w(0) - InterceptAtZero(acov_w)
    else
        w_noise = error
    end if

    do var = u, lastGas
        if (var == v .or. var == w) cycle
        Essentials%rand_uncer(var) = error
        if (.not. E2Col(var)%present) cycle
        if (.not. usable) cycle

        do lag = 0, last_lag
            acov_g(lag) = LaggedCovarianceNoError(Set(:, var), Set(:, var), &
                                                  N, lag, error)
        end do
        if (any(acov_g == error)) cycle

        g_noise = acov_g(0) - InterceptAtZero(acov_g)

        !> Both intercepts must show a noise-like step. A non-positive one
        !> means the fitted line already sits at or above the measured lag-zero
        !> value, which is not a small noise estimate - it is the method
        !> failing to apply, and reporting it as a small number would be worse
        !> than reporting nothing.
        if (w_noise <= 0d0 .or. g_noise <= 0d0) cycle

        !> The scalar's noise variance against the TOTAL variance of w, not
        !> against w's noise. EddyUH's form: what reaches the flux is the
        !> scalar's own noise beaten against the full vertical wind signal.
        Essentials%rand_uncer(var) = dsqrt(g_noise * acov_w(0) / dble(N))
    end do

contains

    !> A straight line through lags first_lag..last_lag, read back at lag
    !> zero. Closed form rather than a general fit: the abscissa is a fixed
    !> short run of integers, so its mean and spread are constants and a
    !> polyfit call would only hide that.
    real(kind = dbl) function InterceptAtZero(acov)
        real(kind = dbl), intent(in) :: acov(0:last_lag)
        integer :: k
        real(kind = dbl) :: xbar, ybar, sxy, sxx, slope

        xbar = dble(first_lag + last_lag) / 2d0
        ybar = 0d0
        do k = first_lag, last_lag
            ybar = ybar + acov(k)
        end do
        ybar = ybar / dble(last_lag - first_lag + 1)

        sxy = 0d0
        sxx = 0d0
        do k = first_lag, last_lag
            sxy = sxy + (dble(k) - xbar) * (acov(k) - ybar)
            sxx = sxx + (dble(k) - xbar)**2
        end do
        slope = sxy / sxx
        InterceptAtZero = ybar - slope * xbar
    end function InterceptAtZero

end subroutine RU_Lenschow_00

!***************************************************************************
!
! \brief       Estimate the flux noise floor after Billesbach (2011), Eq. 3 -
!              the "random shuffle" method.
! \author      Gerardo Fratini, repaired by Jonathan Muller
!
! \note
! Reorder a scalar at random and its covariance with w should vanish: the
! shuffle destroys every real correlation, so whatever covariance survives is
! produced by noise alone. Repeat, and the average of those magnitudes is a
! floor below which a flux cannot be told from zero.
!
! This is NOT a sampling error and is not interchangeable with Finkelstein &
! Sims (2001), which answers the different question of how uncertain a
! measured flux is. Billesbach's estimate is dominated by instrument noise and
! is systematically the smaller of the two.
!
! The statistic is the mean of the ABSOLUTE covariances over the realisations,
! which is what Billesbach specifies and what EddyUH computes
! (EC_Software_Preproc/EddyUH_unc_Preproc.m, lines 99-116). Note that this is
! not a standard deviation: for a zero-mean Gaussian, E|X| is sqrt(2/pi) times
! sigma, so the value runs about four fifths of the scatter it describes.
!
! \sa          RU_Finkelstein_Sims_01, which fills the same slot for the
!              sampling error.
!***************************************************************************
subroutine RU_Billesbach_11(Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    !> local variables
    !> Twenty, which is what EddyUH uses and enough to steady the mean. Fixed
    !> rather than settable: unlike a window width it has no site-dependent
    !> right answer - more realisations simply cost more and wobble less - and
    !> the flat project-tag table has few free slots left to spend on it.
    integer, parameter :: ntimes = 20
    integer :: var
    integer :: i
    real(kind = dbl) :: shuffled(N)
    real(kind = dbl) :: cov
    real(kind = dbl) :: acc
    real(kind = dbl), external :: LaggedCovarianceNoError

    call SeedShuffleOnce()

    !> Same set of variables Finkelstein & Sims fills, for the same reason:
    !> these are the covariances with w that become fluxes.
    do var = u, lastGas
        if (var == v .or. var == w) cycle
        if (.not. E2Col(var)%present) then
            Essentials%rand_uncer(var) = error
            cycle
        end if

        acc = 0d0
        do i = 1, ntimes
            !> The scalar is shuffled, not w. Shuffling w once and reusing it
            !> for every gas would be cheaper, and is what the unreachable
            !> version of this routine did, but it makes every gas's estimate
            !> a draw from the same realisation - they would rise and fall
            !> together for no physical reason. Billesbach and EddyUH both
            !> shuffle the scalar, once per gas per realisation.
            call RandomShuffle(Set(1:N, var), shuffled, N)
            cov = LaggedCovarianceNoError(Set(1:N, w), shuffled, N, 0, error)
            if (cov == error) then
                acc = error
                exit
            end if
            acc = acc + dabs(cov)
        end do

        if (acc /= error) then
            Essentials%rand_uncer(var) = acc / dble(ntimes)
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RU_Billesbach_11

!***************************************************************************
!
! \brief       Put the shuffle's generator into a known state, once per run.
! \author      Jonathan Muller
!
! \note
! Without this the sequence comes from whatever gfortran seeds itself with,
! which since GCC 7 is drawn from the operating system: the same project over
! the same data would then publish a different noise floor every time it ran.
! A fixed seed makes the run reproducible, and because it is set once rather
! than per period, successive periods still draw independent permutations.
!
! Deliberately not a setting. It is not a knob anyone should want to turn -
! two different seeds are two equally valid answers - and the point of fixing
! it is that the number in the output file can be checked by running again.
!
! random_number is used nowhere else in the engine, so seeding it here reaches
! nothing but the shuffle.
!***************************************************************************
subroutine SeedShuffleOnce()
    use m_rp_global_var
    implicit none
    !> Arbitrary, and only has to stay put. Shares the spirit of PWBSetup's
    !> default seed without borrowing its value, which belongs to a different
    !> method.
    integer, parameter :: base = 20110
    integer :: n
    integer :: i
    integer, allocatable :: seed(:)
    logical, save :: seeded = .false.

    if (seeded) return
    call random_seed(size = n)
    allocate(seed(n))
    seed = [(base + 37 * i, i = 0, n - 1)]
    call random_seed(put = seed)
    deallocate(seed)
    seeded = .true.
end subroutine SeedShuffleOnce
!***************************************************************************
!
! \brief       shuffle array elements randomly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
subroutine RandomShuffle(arr, arrout, N)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: arr(N)
    real(kind = dbl), intent(out) :: arrout(N)
    !> Local variables
    integer :: work
    integer :: ix(N)
    integer :: i
    integer :: j
    integer, external :: RandomBetween

    !> Create array of indexes from 1 to size(arr)
    do i = 1, size(arr)
        ix(i) = i
    end do

    !> Shuffle indexes
    do i = N, 2, -1
        j = RandomBetween(1, i)
        !> swap
        work = ix(j)
        ix(j) = ix(i)
        ix(i) = work
    end do
    !> Assign to shuffled array
    arrout = arr(ix)
end subroutine RandomShuffle

!***************************************************************************
!
! \brief       Generate a random integer in [min, max], both ends included.
! \author      Gerardo Fratini
!
! \note
! Both ends, which it did not use to give. The form was
! `int((max - min) * x + min)` with x in [0,1), which spans min to max-1 and
! can never return max. Its only caller is the Fisher-Yates loop in
! RandomShuffle, where j must be able to equal i or the element at i is
! guaranteed to move: that draws uniformly from the permutations in which no
! element stays put, not from all of them, and the shuffle it produces is
! biased. Harmless while the shuffle was unreachable, and not while it feeds
! a published noise floor.
!***************************************************************************
integer function RandomBetween(min, max)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: min
    integer, intent(in) :: max
    real(kind = dbl) :: x

    call random_number(x)
    RandomBetween = min + int(dble(max - min + 1) * x)
    !> x < 1 by contract, so this cannot trigger; it is here because the cost
    !> of being wrong is an out-of-bounds swap in the caller.
    if (RandomBetween > max) RandomBetween = max
end function RandomBetween

!***************************************************************************
!
! \brief       Initialize random number generation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
subroutine InitRandomSeed()
    use m_rp_global_var
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, dt(8), pid, t(2), s
    integer(8) :: count, tms


    call random_seed(size = n)
    allocate(seed(n))

    !> XOR:ing the current time and pid. The PID is
    !> useful in case one launches multiple instances of the same
    !> program in parallel.
    call system_clock(count)
    if (count /= 0) then
        t = transfer(count, t)
    else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
        t = transfer(tms, t)
    end if
    s = ieor(t(1), t(2))
    pid = getpid() + 1099279 ! Add a prime
    s = ieor(s, pid)
    if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
    else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    end if
    call random_seed(put=seed)
end subroutine InitRandomSeed
