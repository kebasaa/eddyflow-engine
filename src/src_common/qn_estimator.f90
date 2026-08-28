!***************************************************************************
! qn_estimator.f90
! ----------------
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
! \brief       Qn scale estimator and helper functions
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Rousseeuw, P.J. & Croux, C. (1993). Alternatives to the Median
!            Absolute Deviation. J. Amer. Statist. Assoc. 88(424):1273-1283.

! Qn(x,n): Qn = dn * 2.2219 * {|xi-xj|; i<j}_(k), k = C(h,2), h=floor(n/2)+1.
double precision function Qn(x, n)
    use m_numeric_kinds
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    integer :: ii, jj, half_n, k_rank, npairs, ip
    real(kind = dbl) :: dn
    real(kind = dbl), allocatable :: pair_diffs(:)
    double precision, external :: pull

    ! Reference: Rousseeuw & Croux (1993) JASA 88:1273-1283, Eq. (1)
    half_n = n / 2 + 1
    k_rank = half_n * (half_n - 1) / 2
    npairs = n * (n - 1) / 2
    allocate(pair_diffs(npairs))
    ip = 0
    do ii = 1, n - 1
        do jj = ii + 1, n
            ip = ip + 1
            pair_diffs(ip) = dabs(x(ii) - x(jj))
        end do
    end do
    Qn = pull(pair_diffs, npairs, k_rank)
    deallocate(pair_diffs)

    ! Small-sample correction dn (Rousseeuw & Croux 1993, Table 1)
    select case (n)
        case (2);    dn = 0.399d0
        case (3);    dn = 0.994d0
        case (4);    dn = 0.512d0
        case (5);    dn = 0.844d0
        case (6);    dn = 0.611d0
        case (7);    dn = 0.857d0
        case (8);    dn = 0.669d0
        case (9);    dn = 0.872d0
        case default
            if (mod(n, 2) == 1) then
                dn = dble(n) / (dble(n) + 1.4d0)
            else
                dn = dble(n) / (dble(n) + 3.8d0)
            end if
    end select
    ! Asymptotic constant 2.2219 = 1/Phi^{-1}(5/8), Rousseeuw & Croux 1993 Eq. (1)
    Qn = dn * 2.2219d0 * Qn
end function Qn


! Weighted high median: smallest a(j) such that cumulative weight > total/2.
! Uses insertion sort on a local copy; independent of the StatLib original.
double precision function whimed(a, iw, n)
    use m_numeric_kinds
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: a(n)
    integer, intent(inout) :: iw(n)
    real(kind = dbl) :: vbuf(n), vtmp
    integer :: wbuf(n), itmp, spos, ipos, wtot, wcum

    vbuf = a;  wbuf = iw
    do spos = 2, n
        vtmp = vbuf(spos);  itmp = wbuf(spos);  ipos = spos - 1
        do while (ipos >= 1 .and. vbuf(ipos) > vtmp)
            vbuf(ipos+1) = vbuf(ipos);  wbuf(ipos+1) = wbuf(ipos)
            ipos = ipos - 1
        end do
        vbuf(ipos+1) = vtmp;  wbuf(ipos+1) = itmp
    end do
    wtot = sum(wbuf);  wcum = 0;  whimed = vbuf(n)
    do spos = 1, n
        wcum = wcum + wbuf(spos)
        if (2*wcum > wtot) then
            whimed = vbuf(spos);  return
        end if
    end do
end function whimed


! Shell sort: sorts a(1:n) into b(1:n).
subroutine qn_sort(a, n, b)
    use m_numeric_kinds
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: a(n)
    real(kind = dbl), intent(inout) :: b(n)
    integer :: gap, pos, cmp
    real(kind = dbl) :: tmp

    b = a;  gap = n / 2
    do while (gap > 0)
        do pos = gap + 1, n
            tmp = b(pos);  cmp = pos
            do while (cmp > gap .and. b(cmp-gap) > tmp)
                b(cmp) = b(cmp-gap);  cmp = cmp - gap
            end do
            b(cmp) = tmp
        end do
        gap = gap / 2
    end do
end subroutine qn_sort


! Quickselect: returns k-th order statistic of a(1:n) in expected O(n) time.
double precision function pull(a, n, k)
    use m_numeric_kinds
    implicit none
    integer, intent(in) :: n, k
    real(kind = dbl), intent(in) :: a(n)
    real(kind = dbl) :: buf(n), piv, tmp
    integer :: lo, hi, sl, sh

    buf = a;  lo = 1;  hi = n
    do while (lo < hi)
        piv = buf((lo+hi)/2);  sl = lo;  sh = hi
        do while (sl <= sh)
            do while (buf(sl) < piv); sl = sl + 1; end do
            do while (buf(sh) > piv); sh = sh - 1; end do
            if (sl <= sh) then
                tmp = buf(sl);  buf(sl) = buf(sh);  buf(sh) = tmp
                sl = sl + 1;  sh = sh - 1
            end if
        end do
        if (sh < k) lo = sl
        if (k < sl) hi = sh
    end do
    pull = buf(k)
end function pull
