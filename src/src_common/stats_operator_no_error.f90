!***************************************************************************
! stats_operator_no_error.f90
! ---------------------------
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
! \brief       Calculate column-wise sums on a 2d array \n
!              ignoring specified error values \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SumNoError(Set, nrow, ncol, Summ, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Summ(ncol)
    !> local variables
    integer :: i
    integer :: j
    logical :: data_exist

    Summ = 0d0
    do j = 1, ncol
        data_exist = .false.
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                data_exist = .true.
                Summ(j) = Summ(j) + Set(i, j)
            end if
        end do
        if (.not. data_exist) Summ(j) = err_float
    end do
end subroutine SumNoError

!***************************************************************************
!
! \brief       Calculate column-wise averages on a 2d array \n
!              ignoring specified error values \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AverageNoError(Set, nrow, ncol, Mean, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Mean(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: RawMean(ncol)


    RawMean = 0d0
    do j = 1, ncol
        Nact = 0
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                Nact = Nact + 1
                RawMean(j) = RawMean(j) + Set(i, j)
            end if
        end do
        if (Nact /= 0) then
            RawMean(j) = RawMean(j) / dble(Nact)
        else
            RawMean(j) = err_float
        end if
    end do
    Mean = 0.d0
    do j = 1, ncol
        if (RawMean(j) /= err_float) then
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= err_float) then
                    Nact = Nact + 1
                    Mean(j) = Mean(j) + Set(i, j) - RawMean(j)
                end if
            end do
            if (Nact /= 0) then
                Mean(j) = Mean(j) / dble(Nact)
            else
                Mean(j) = err_float
            end if
        else
            Mean(j) = err_float
        end if
    end do

    where (Mean(:) /= err_float)
        Mean(:) = Mean(:) + RawMean(:)
    elsewhere
        Mean(:) = err_float
    end where
end subroutine AverageNoError

!***************************************************************************
!
! \brief       Calculate column-wise angular averages on a 2d array \n
!              ignoring specified error values. In EddyFlow, mainly meant for \n
!              calculation of mean wind direction given a set of wind direction \n
!              measurements.
!
!              Implementation reference:
!              "Circular Statistics in R"
!              by A. Pewsey, M. Neuhaeuser and G. D. Ruxton.
!              Yamartino, 1984:
!              https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AngularAverageNoError(Set, nrow, ncol, Mean, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Mean(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: CosSum
    real(kind = dbl) :: SinSum


    do j = 1, ncol

        !> Calculate a (CosSum) and b (SinSum)
        CosSum = 0d0
        SinSum = 0d0
        Nact = 0
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                Nact = Nact + 1
                CosSum = CosSum - dcos(Set(i, j) / 180d0 * p)
                SinSum = SinSum - dsin(Set(i, j) / 180d0 * p)
            end if
        end do
        if (Nact /= 0) then
            CosSum = CosSum / dble(Nact)
            SinSum = SinSum / dble(Nact)
        else
            Mean(j) = err_float
            cycle
        end if

        !> Angular average is atan2 of b and a
        !> "+p" adjust quadrant, then express in degrees
        Mean(j) = (datan2(SinSum, CosSum) + p) * 180d0 / p
    end do
end subroutine AngularAverageNoError

!***************************************************************************
!
! \brief       Calculate column-wise angular stdev on a 2d array \n
!              ignoring specified error values. In EddyFlow, mainly meant for \n
!              calculation of wind direction standard deviation given a set of wind direction. \n
!
!              Implementation reference:
!              Yamartino, 1984:
!              https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!
! \brief       Angular standard deviation (Yamartino 1984 single-pass approx.)
!              Ignores error-coded values.
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Yamartino, R.J. (1984). A comparison of several "single-pass"
!            estimators of the standard deviation of wind direction.
!            J. Climate Appl. Meteor. 23: 1362-1366.
subroutine AngularStDevApproxNoError(Set, nrow, ncol, AngStDev, err_float)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: AngStDev(ncol)
    integer :: icol, irow, nvalid
    real(kind = dbl) :: cos_sum, sin_sum, circ_r2, eps_val

    do icol = 1, ncol
        cos_sum = 0d0;  sin_sum = 0d0;  nvalid = 0
        do irow = 1, nrow
            if (Set(irow, icol) /= err_float) then
                nvalid = nvalid + 1
                cos_sum = cos_sum + dcos(Set(irow, icol) / 180d0 * p)
                sin_sum = sin_sum + dsin(Set(irow, icol) / 180d0 * p)
            end if
        end do
        if (nvalid == 0) then
            AngStDev(icol) = err_float;  cycle
        end if
        cos_sum = cos_sum / dble(nvalid)
        sin_sum = sin_sum / dble(nvalid)
        ! Yamartino (1984) Eq. (3)-(4): sigma = arcsin(eps)*(1+C*eps^3)
        circ_r2 = cos_sum**2 + sin_sum**2
        eps_val = dsqrt(max(0d0, 1d0 - circ_r2))
        AngStDev(icol) = (dasin(eps_val) * &
            (1d0 + (2d0/dsqrt(3d0) - 1d0)*eps_val**3)) * 180d0 / p
    end do
end subroutine AngularStDevApproxNoError


!***************************************************************************
!
! \brief       Calculates standard deviation of array (column-wise) ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StDevNoError(Set, nrow, ncol, StDev, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: StDev(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: Mean(ncol)


    !> Initializations
    StDev = 0d0

    !> Calculate mean values
    call AverageNoError(Set, nrow, ncol, Mean, err_float)

    !> Sum of squared residuals
    do j = 1, ncol
        if (Mean(j) == err_float) then
            StDev(j) = err_float
        else
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= err_float) then
                    Nact = Nact + 1
                    StDev(j) = StDev(j) + (Set(i, j) - Mean(j)) **2
                end if
            end do
            if (Nact /= 0 .and. StDev(j) >= 0d0) then
                StDev(j) = dsqrt(StDev(j) / dble(Nact-1))
            else
                StDev = err_float
            end if
        end if
    end do

end subroutine StDevNoError

!***************************************************************************
!
! \brief       Calculates covariance matrix of given array, ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovarianceMatrixNoError(Set, nrow, ncol, Cov, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Cov(ncol, ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: Nact = 0
    real(kind = dbl) :: sumi
    real(kind = dbl) :: sumj

    do i = 1, ncol
        do j = 1, ncol
            sumi = 0d0
            sumj = 0d0
            Cov(i, j) = 0d0
            Nact = 0
            do k = 1, nrow
                if (Set(k, i) /= err_float .and. Set(k, j) /= err_float) then
                    Nact = Nact + 1
                    Cov(i, j) = Cov(i, j) + Set(k, i) * Set(k, j)
                    sumi = sumi + Set(k, i)
                    sumj = sumj + Set(k, j)
                end if
            end do
            if (Nact /= 0) then
                sumi = sumi / dble(Nact)
                sumj = sumj / dble(Nact)
                Cov(i, j) = Cov(i, j) / dble(Nact)
                Cov(i, j) = Cov(i, j) - sumi * sumj
            else
                Cov(i, j) = err_float
            end if
        end do
    end do
end subroutine CovarianceMatrixNoError


!***************************************************************************
!
! \brief       Full pairwise correlation matrix, ignoring error-coded values
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine CorrelationMatrixNoError(Set, nrow, ncol, Corr, err_float)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Corr(ncol, ncol)
    real(kind = dbl) :: cov_mat(ncol, ncol)
    real(kind = dbl) :: sigma(ncol)
    integer :: ci, cj

    call CovarianceMatrixNoError(Set, nrow, ncol, cov_mat, err_float)
    call StDevNoError(Set, nrow, ncol, sigma, err_float)

    do ci = 1, ncol
        do cj = 1, ncol
            if (cov_mat(ci,cj) /= err_float .and. &
                sigma(ci) > 0d0 .and. sigma(cj) > 0d0) then
                Corr(ci, cj) = cov_mat(ci, cj) / (sigma(ci) * sigma(cj))
            else
                Corr(ci, cj) = err_float
            end if
        end do
    end do
end subroutine CorrelationMatrixNoError

!***************************************************************************
!
! \brief       Calculates covariance matrix of given arrays applying \n
!              given lag, ignoring provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function LaggedCovarianceNoError(col1, col2, nrow, rlag, err_float)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: rlag
    real(kind = dbl), intent(in) :: col1(nrow)
    real(kind = dbl), intent(in) :: col2(nrow)
    real(kind = dbl), intent(in) :: err_float
    !> Local variables
    integer :: lag
    integer :: i
    integer :: n
    real(kind = dbl) :: cov
    real(kind = dbl) :: sumi
    real(kind = dbl) :: sumj


    sumi = 0d0
    sumj = 0d0
    cov = 0d0
    n = 0
    if (rlag >= 0) then
        !> Positive lags are interpreted as col2 being "late"
        do i = 1, nrow - rlag
            if (col1(i) /= err_float .and. col2(i+rlag) /= err_float) then
                n = n + 1
                cov = cov + col1(i) * col2(i+rlag)
                sumi = sumi + col1(i)
                sumj = sumj + col2(i+rlag)
            end if
        end do
    else
        !> Positive lags are interpreted as col1 being "late"
        lag = -rlag
        do i = 1, nrow - lag
            if (col2(i) /= err_float .and. col1(i+lag) /= err_float) then
                n = n + 1
                cov = cov + col2(i) * col1(i+lag)
                sumi = sumi + col2(i)
                sumj = sumj + col1(i+lag)
            end if
        end do
    end if

    !> Finish up
    if (n /= 0) then
        sumi = sumi / dble(n)
        sumj = sumj / dble(n)
        cov = cov / dble(n)
        LaggedCovarianceNoError = cov - sumi * sumj
    else
        LaggedCovarianceNoError = err_float
    end if
end function LaggedCovarianceNoError

!***************************************************************************
!
! \brief       Column-wise skewness, ignoring error-coded values
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine SkewnessNoError(Set, nrow, ncol, Skw, err_float)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol), err_float
    real(kind = dbl), intent(out) :: Skw(ncol)
    integer :: icol, irow, nvalid
    real(kind = dbl) :: sigma(ncol), col_mean(ncol), cube_sum

    call AverageNoError(Set, nrow, ncol, col_mean, err_float)
    call StDevNoError(Set, nrow, ncol, sigma, err_float)

    Skw = err_float
    do icol = 1, ncol
        if (sigma(icol) == err_float .or. sigma(icol) == 0d0) cycle
        cube_sum = 0d0;  nvalid = 0
        do irow = 1, nrow
            if (Set(irow, icol) /= err_float) then
                nvalid = nvalid + 1
                cube_sum = cube_sum + (Set(irow, icol) - col_mean(icol))**3
            end if
        end do
        if (nvalid > 1) Skw(icol) = cube_sum / (sigma(icol)**3 * dble(nvalid))
    end do
end subroutine SkewnessNoError


!***************************************************************************
!
! \brief       Column-wise kurtosis, ignoring error-coded values
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine KurtosisNoError(Set, nrow, ncol, Kur, err_float)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol), err_float
    real(kind = dbl), intent(out) :: Kur(ncol)
    integer :: icol, irow, nvalid
    real(kind = dbl) :: sigma(ncol), col_mean(ncol), quad_sum

    call AverageNoError(Set, nrow, ncol, col_mean, err_float)
    call StDevNoError(Set, nrow, ncol, sigma, err_float)

    Kur = err_float
    do icol = 1, ncol
        if (sigma(icol) == err_float .or. sigma(icol) == 0d0) cycle
        quad_sum = 0d0;  nvalid = 0
        do irow = 1, nrow
            if (Set(irow, icol) /= err_float) then
                nvalid = nvalid + 1
                quad_sum = quad_sum + (Set(irow, icol) - col_mean(icol))**4
            end if
        end do
        if (nvalid > 1) Kur(icol) = quad_sum / (sigma(icol)**4 * dble(nvalid))
    end do
end subroutine KurtosisNoError


!***************************************************************************
!
! \brief       Column-wise quantile, ignoring error-coded values
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine QuantileNoError(Set, nrow, ncol, Quantile, qin, err_float)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float, qin
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Quantile(ncol)
    integer :: icol, irow, nvalid
    real(kind = dbl), allocatable :: vals(:)
    real(kind = dbl), external :: quantile_sas5

    Quantile = err_float
    do icol = 1, ncol
        nvalid = count(Set(:, icol) /= err_float)
        if (nvalid > 1) then
            allocate(vals(nvalid))
            nvalid = 0
            do irow = 1, nrow
                if (Set(irow, icol) /= err_float) then
                    nvalid = nvalid + 1;  vals(nvalid) = Set(irow, icol)
                end if
            end do
            Quantile(icol) = quantile_sas5(vals, nvalid, qin)
            deallocate(vals)
        end if
    end do
end subroutine QuantileNoError


subroutine unbiased_correlation(arr1, arr2, n, err_float, lag, r, t, m)
    ! Unbiased Pearson correlation of arr1 vs arr2 at given lag,
    ! skipping error-coded values.
    implicit none
    integer, intent(in) :: n, lag
    real, dimension(n), intent(in) :: arr1, arr2
    real, intent(in) :: err_float
    real, intent(out) :: r, t
    integer, intent(out) :: m
    real :: xmn, ymn, covxy, sx, sy
    real, allocatable :: xa(:), ya(:)

    allocate(xa(3*n), ya(3*n))
    xa = err_float;  ya = err_float
    xa(n+1:2*n) = arr1;  ya(n+1:2*n) = arr2
    ya = eoshift(ya, shift = -lag, boundary = err_float)
    where (xa == err_float) ya = err_float
    where (ya == err_float) xa = err_float
    m = count(xa /= err_float)
    if (m < 2) then
        r = err_float;  t = err_float;  deallocate(xa, ya);  return
    end if
    xmn = sum(xa, mask = xa /= err_float) / float(m)
    ymn = sum(ya, mask = ya /= err_float) / float(m)
    covxy = sum((xa-xmn)*(ya-ymn), &
        mask = xa /= err_float .and. ya /= err_float) / float(m)
    sx = sqrt(sum((xa-xmn)**2, mask = xa /= err_float) / float(m))
    sy = sqrt(sum((ya-ymn)**2, mask = ya /= err_float) / float(m))
    if (sx > 0.0 .and. sy > 0.0) then
        r = covxy / (sx * sy)
    else
        r = err_float;  t = err_float;  deallocate(xa, ya);  return
    end if
    if (abs(r) < 1.0) then
        t = r * sqrt(float(m-2) / (1.0 - r*r))
    else
        t = err_float
    end if
    deallocate(xa, ya)
end subroutine unbiased_correlation


!***************************************************************************
!
! \brief       Cross-correlation function for a specified lag range
!
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine CrossCorrelation(arr1, arr2, nrow, lagmin, lagmax, CCF)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, lagmin, lagmax
    real(kind = dbl), intent(in) :: arr1(nrow), arr2(nrow)
    real(kind = dbl), intent(out) :: CCF(lagmin:lagmax)
    integer :: lag
    real(kind = dbl) :: sig1(1), sig2(1)
    real(kind = dbl), external :: LaggedCovarianceNoError

    do lag = lagmin, lagmax
        CCF(lag) = LaggedCovarianceNoError(arr1, arr2, nrow, lag, error)
    end do
    call StDevNoError(arr1, nrow, 1, sig1, error)
    call StDevNoError(arr2, nrow, 1, sig2, error)
    if (sig1(1) > 0d0 .and. sig2(1) > 0d0) then
        CCF = CCF / (sig1(1) * sig2(1))
    else
        CCF = error
    end if
end subroutine CrossCorrelation