!***************************************************************************
! cross_corr_test.f90
! -------------------
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
! \brief       CCF degradation test: R^2 of the raw w-vs-variable CCF
!              against its de-flatlined counterpart, per variable.
! \author      Jonathan Muller, ETH Zurich
! \note        Ported from RFlux's qcStat.R (the lrt_h/lrt_fc/lrt_le/lrt_tau
!              block), Vitale, D. et al. (2020), Biogeosciences Discussion.
!              RFlux computes this only as an exported diagnostic column -
!              it is never one of cleanFlux.R's SevEr/ModEr tests, so this
!              stays purely informational here too rather than feeding
!              VitaleFlag, which would be inventing a threshold RFlux
!              itself never defines.
!
!              Stored per variable paired with w (CCF(u), CCF(v), CCF(ts),
!              CCF(gas)...), unlike RFlux, which only exports one combined
!              lrt_tau = max(lrt_wu, lrt_wv) for momentum flux. Keeping
!              CCF(u) and CCF(v) separate is a strict superset of that
!              information, not a divergence from it - lrt_tau is
!              max(CCF(u), CCF(v)) whenever both are available.
!***************************************************************************
subroutine CrossCorrTest(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer, parameter :: lagmin = -25
    integer, parameter :: lagmax = 25
    integer, parameter :: nlag = lagmax - lagmin + 1
    integer :: icol, irec, d0
    real(kind = dbl) :: dedup_set(nrow, ncol)
    integer :: n_flat(ncol)
    real(kind = dbl) :: raw_ccf(lagmin:lagmax), dup_ccf(lagmin:lagmax)
    real(kind = dbl), external :: CCFDegradationR2

    write(*, '(a)', advance = 'no') &
        '  Evaluating R2 on CCFs with and without repeated values..'
    write(ulog, '(a)', advance = 'no') &
        '  Evaluating R2 on CCFs with and without repeated values..'

    !> Same flat-line marking RFlux's own ind_w <- which(diff(x)==0)+1
    !> does, plus (needed for the D0 threshold below, which RFlux's own
    !> version never had to compute separately) a running count of how
    !> many samples each column had marked.
    dedup_set = Set
    n_flat = 0
    do icol = u, lastGas
        if (OutVarPresent(icol)) then
            do irec = 2, nrow
                if (dabs(Set(irec, icol) - Set(irec-1, icol)) < 1d-8) then
                    dedup_set(irec, icol) = error
                    n_flat(icol) = n_flat(icol) + 1
                end if
            end do
        end if
    end do

    !> w has no self-pair.
    Essentials%CCF(w) = error

    !> H and every gas: CCF(w, icol).
    do icol = ts, lastGas
        if (OutVarPresent(icol)) then
            call CrossCorrelation(Set(:, w), Set(:, icol), nrow, lagmin, lagmax, raw_ccf)
            call CrossCorrelation(dedup_set(:, w), dedup_set(:, icol), &
                                  nrow, lagmin, lagmax, dup_ccf)
            d0 = max(n_flat(w), n_flat(icol))
            Essentials%CCF(icol) = CCFDegradationR2(raw_ccf, dup_ccf, nlag, d0, nrow)
        else
            Essentials%CCF(icol) = error
        end if
    end do

    !> TAU: CCF(w, u) and CCF(w, v), kept separate - see this file's own
    !> header for why that is a superset of RFlux's single lrt_tau, not a
    !> divergence from it.
    do icol = u, v
        if (OutVarPresent(icol)) then
            call CrossCorrelation(Set(:, w), Set(:, icol), nrow, lagmin, lagmax, raw_ccf)
            call CrossCorrelation(dedup_set(:, w), dedup_set(:, icol), &
                                  nrow, lagmin, lagmax, dup_ccf)
            d0 = max(n_flat(w), n_flat(icol))
            Essentials%CCF(icol) = CCFDegradationR2(raw_ccf, dup_ccf, nlag, d0, nrow)
        else
            Essentials%CCF(icol) = error
        end if
    end do

    call LogSayList(' Done.')
end subroutine CrossCorrTest


!***************************************************************************
!
! \brief       R^2 of a through-origin regression of one CCF against its
!              de-flatlined counterpart, with RFlux's own D0-based
!              branching (qcStat.R: ifelse(D0 < N*0.9, ifelse(D0 > 0,
!              summary(lm(CORori~CORsub-1))$r.squared, 1), -1)).
! \author      Jonathan Muller, ETH Zurich
! \note        R's summary.lm() for a no-intercept fit uses an uncentered
!              model sum of squares (mss = sum(fitted^2), not sum around
!              the mean), so R^2 = mss/(mss+rss) reduces algebraically to
!              (sum(x*y))^2 / (sum(x^2)*sum(y^2)) - the squared cosine
!              similarity between the two CCF vectors, not the ordinary
!              (mean-centered) squared correlation a naive port would use.
!***************************************************************************
real(kind = dbl) function CCFDegradationR2(y, x, nlag, d0, ntot) result(r2)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nlag, d0, ntot
    real(kind = dbl), intent(in) :: y(nlag), x(nlag)
    real(kind = dbl) :: sxy, sxx, syy

    !> Too much of the wider column was flat-lined for the de-flatlined
    !> CCF to mean anything.
    if (dble(d0) >= 0.9d0 * dble(ntot)) then
        r2 = -1d0
        return
    end if
    !> Neither column had anything flat-lined to begin with: the two CCFs
    !> are identical, so RFlux short-circuits to a perfect 1 rather than
    !> fitting a 0/0 regression.
    if (d0 == 0) then
        r2 = 1d0
        return
    end if
    if (any(y == error) .or. any(x == error)) then
        r2 = error
        return
    end if

    sxy = sum(x * y)
    sxx = sum(x * x)
    syy = sum(y * y)
    if (sxx <= 0d0 .or. syy <= 0d0) then
        r2 = error
        return
    end if
    r2 = (sxy * sxy) / (sxx * syy)
end function CCFDegradationR2
