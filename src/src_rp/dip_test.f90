!***************************************************************************
! dip_test.f90
! ------------
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
! \brief       Hartigan's dip statistic and its p-value, RFlux's DIP
!              diagnostic (Vitale et al. 2020, Biogeosciences).
! \author      Jonathan Muller, ETH Zurich
! \note
! Reference for the statistic itself:
!   Hartigan, J. A. & Hartigan, P. M. (1985). The Dip Test of Unimodality.
!   Annals of Statistics, 13(1), 70-84.
!   Hartigan, P. M. (1985). Computation of the Dip Statistic to Test for
!   Unimodality. Applied Statistics (AS 217), 34(3), 320-325.
!
! DipStatistic below is a direct port of the greatest-convex-minorant /
! least-concave-majorant algorithm as given in RUrlus/diptest's
! src/diptest-core/include/diptest/diptest.hpp (function `diptst`, MIT/GPL,
! itself a C++ rewrite of R diptest's C code, itself an f2c translation of
! the original AS 217 Fortran) - ported here to Fortran rather than
! reimplemented from the paper's description, because the algorithm's
! bookkeeping (which index is "low", which change points survive a cycle)
! is exactly the kind of detail an independent re-derivation gets subtly
! wrong in ways a compiler cannot catch. Validated instead: this port's
! output matches the `diptest` PyPI package (itself validated against R's
! diptest) to machine precision (~1e-16) across 336 test vectors spanning
! n=1..36000, unimodal/bimodal/trimodal/constant/heavily-tied data - see
! the session that added this file for the validation script.
!
! DipPvalue reproduces the same package's default p-value method: linear
! interpolation, on sqrt(n)*dip, of the tabulated critical-value grid
! (diptest/consts.py's _SAMPLE_SIZE x _ALPHA table, the same qDiptab R's
! own dip.test() interpolates - which is what RFlux calls). One behaviour
! had to be reverse-engineered rather than read off the docs: numpy's
! np.interp resolves an x that exactly matches a run of duplicate grid
! nodes (common at the table's small-n end) to the LAST matching node, not
! the first - confirmed by testing against its actual output. LinInterp1D
! below reproduces that by scanning brackets from the high end.
!***************************************************************************
subroutine DipStatistic(x, n, dip)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)  !< must already be sorted ascending
    real(kind = dbl), intent(out) :: dip

    integer :: mn(n), mj(n), gcm(n), lcm(n)
    integer :: low, high, i
    integer :: gcm_rel, lcm_rel, gcm_x, gcm_y, lcm_x, lcm_y
    real(kind = dbl) :: running_dip, d, dip_l_val, dip_u_val, tmp_val
    logical :: flag_done

    low = 1
    high = n
    running_dip = 0d0

    if (n < 2) then
        dip = 0d0
        return
    end if
    if (x(n) == x(1)) then
        dip = 0d0
        return
    end if

    call DipIndicesMinorant(x, n, mn)
    call DipIndicesMajorant(x, n, mj)

    do
        gcm(1) = high
        i = 1
        do while (gcm(i) > low)
            gcm(i+1) = mn(gcm(i))
            i = i + 1
        end do
        gcm_rel = i
        gcm_x = gcm_rel
        gcm_y = gcm_rel - 1

        lcm(1) = low
        i = 1
        do while (lcm(i) < high)
            lcm(i+1) = mj(lcm(i))
            i = i + 1
        end do
        lcm_rel = i
        lcm_x = lcm_rel
        lcm_y = 2

        if (gcm_rel /= 2 .or. lcm_rel /= 2) then
            call DipMaxDistance(x, n, gcm, gcm_rel, gcm_y, gcm_x, lcm, lcm_rel, lcm_y, lcm_x, d)
        else
            d = 0d0
        end if

        if (d < running_dip) exit

        call DipEnvelopeDip(x, n, gcm, gcm_rel, gcm_x, 0, 1, dip_l_val)
        call DipEnvelopeDip(x, n, lcm, lcm_rel, lcm_x, 1, -1, dip_u_val)

        if (dip_l_val < dip_u_val) then
            tmp_val = dip_u_val
        else
            tmp_val = dip_l_val
        end if
        if (running_dip < tmp_val) running_dip = tmp_val

        flag_done = (low == gcm(gcm_x) .and. high == lcm(lcm_x))
        low = gcm(gcm_x)
        high = lcm(lcm_x)

        if (flag_done) exit
    end do

    dip = running_dip / dble(2 * n)
end subroutine DipStatistic


!> Change-point indices for the greatest convex minorant of the (sorted)
!> empirical CDF, walking forward from the left.
subroutine DipIndicesMinorant(arr, n, mn)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: arr(n)
    integer, intent(out) :: mn(n)
    integer :: i, ind_at_i, ind_at_i_iter
    logical :: rate_change_flag

    mn(1) = 1
    do i = 2, n
        mn(i) = i - 1
        do
            ind_at_i = mn(i)
            ind_at_i_iter = mn(ind_at_i)
            rate_change_flag = (arr(i)-arr(ind_at_i))*dble(ind_at_i-ind_at_i_iter) < &
                                (arr(ind_at_i)-arr(ind_at_i_iter))*dble(i-ind_at_i)
            if (ind_at_i == 1 .or. rate_change_flag) exit
            mn(i) = ind_at_i_iter
        end do
    end do
end subroutine DipIndicesMinorant


!> Change-point indices for the least concave majorant, walking backward
!> from the right.
subroutine DipIndicesMajorant(arr, n, mj)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: arr(n)
    integer, intent(out) :: mj(n)
    integer :: i, ind_at_i, ind_at_i_iter
    logical :: rate_change_flag

    mj(n) = n
    do i = n - 1, 1, -1
        mj(i) = i + 1
        do
            ind_at_i = mj(i)
            ind_at_i_iter = mj(ind_at_i)
            rate_change_flag = (arr(i)-arr(ind_at_i))*dble(ind_at_i-ind_at_i_iter) < &
                                (arr(ind_at_i)-arr(ind_at_i_iter))*dble(i-ind_at_i)
            if (ind_at_i == n .or. rate_change_flag) exit
            mj(i) = ind_at_i_iter
        end do
    end do
end subroutine DipIndicesMajorant


!> Greatest distance between the current GCM and LCM piecewise-linear fits,
!> walking the two envelopes inward from (gcm_y, lcm_y) and recording where
!> the maximum occurred in (gcm_x, lcm_x).
subroutine DipMaxDistance(arr, n, gcm, gcm_rel, gcm_y, gcm_x, lcm, lcm_rel, lcm_y, lcm_x, ret_d)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: arr(n)
    integer, intent(in) :: gcm(n), lcm(n)
    integer, intent(in) :: gcm_rel, lcm_rel
    integer, intent(inout) :: gcm_y, lcm_y, gcm_x, lcm_x
    real(kind = dbl), intent(out) :: ret_d
    integer :: gcm_yv, lcm_yv, is_maj, i, j, i1, sign_

    ret_d = 0d0
    do
        gcm_yv = gcm(gcm_y)
        lcm_yv = lcm(lcm_y)
        if (gcm_yv > lcm_yv) then
            is_maj = 1
        else
            is_maj = 0
        end if

        if (is_maj == 1) then
            i = gcm_yv
            j = lcm_yv
            i1 = gcm(gcm_y + 1)
            sign_ = 1
        else
            i = lcm_yv
            j = gcm_yv
            i1 = lcm(lcm_y - 1)
            sign_ = -1
        end if

        block
            real(kind = dbl) :: dx
            dx = dble(sign_) * ( dble(j - i1 + sign_) - (arr(j)-arr(i1))*dble(i-i1)/(arr(i)-arr(i1)) )

            if (is_maj == 0) gcm_y = gcm_y - 1
            if (is_maj == 1) lcm_y = lcm_y + 1

            if (dx >= ret_d) then
                ret_d = dx
                gcm_x = gcm_y + 1
                lcm_x = lcm_y - is_maj
            end if
        end block

        if (gcm_y < 1) gcm_y = 1
        if (lcm_y > lcm_rel) lcm_y = lcm_rel

        if (gcm(gcm_y) == lcm(lcm_y)) exit
    end do
end subroutine DipMaxDistance


!> The dip contribution of one envelope (GCM when offset=0/sign=+1, LCM
!> when offset=1/sign=-1) over its change-point segments from x_start
!> onward. Local max-trackers seeded at (0, then 1 per segment) rather than
!> 0 throughout - not a bug: this works in un-normalised 2n*dip units, where
!> even a perfectly unimodal sample has a discreteness floor near 1.
subroutine DipEnvelopeDip(arr, n, optimum, rel_length, x_start, offset, sign_, dip_val)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: arr(n)
    integer, intent(in) :: optimum(n)
    integer, intent(in) :: rel_length, x_start, offset, sign_
    real(kind = dbl), intent(out) :: dip_val
    integer :: j, j_start, j_end, jj
    real(kind = dbl) :: C, arr_j_start, d
    real(kind = dbl) :: tmp_val

    dip_val = 0d0
    tmp_val = 1d0

    do j = x_start, rel_length - 1
        j_start = optimum(j + 1 - offset)
        j_end   = optimum(j + offset)

        if (j_end - j_start > 1 .and. arr(j_end) /= arr(j_start)) then
            C = dble(j_end - j_start) / (arr(j_end) - arr(j_start))
            arr_j_start = arr(j_start)
            do jj = j_start, j_end
                d = dble(sign_) * ( dble(jj - j_start + sign_) - (arr(jj)-arr_j_start)*C )
                if (tmp_val < d) tmp_val = d
            end do
        end if

        if (dip_val < tmp_val) dip_val = tmp_val
        tmp_val = 1d0
    end do
end subroutine DipEnvelopeDip


!> p-value via interpolation of RFlux/R diptest's tabulated critical
!> values (diptest PyPI package's consts.py, itself matching R diptest's
!> qDiptab), the same table R's dip.test() interpolates by default.
subroutine DipPvalue(n, dip, pval)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: dip
    real(kind = dbl), intent(out) :: pval


    integer, parameter :: n_tab_sizes = 21, n_tab_alphas = 26
    integer :: sample_size(n_tab_sizes)
    real(kind = dbl) :: alpha(n_tab_alphas)
    real(kind = dbl) :: crit_vals(n_tab_alphas, n_tab_sizes)
    integer :: i, k, i0, i1, n0, n1
    real(kind = dbl) :: y0(n_tab_alphas), y1(n_tab_alphas), xp(n_tab_alphas), frac_n, sD
    real(kind = dbl), external :: DipLinInterp1D

    data sample_size / 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200, 500, &
                        1000, 2000, 5000, 10000, 20000, 40000, 72000 /

    data alpha / 0.0D0, 0.01D0, 0.02D0, 0.05D0, 0.1D0, 0.2D0, 0.3D0, 0.4D0, &
                 0.5D0, 0.6D0, 0.7D0, 0.8D0, 0.9D0, 0.95D0, 0.98D0, 0.99D0, &
                 0.995D0, 0.998D0, 0.999D0, 0.9995D0, 0.9998D0, 0.9999D0, &
                 0.99995D0, 0.99998D0, 0.99999D0, 1.0D0 /

    !> BEGIN GENERATED dip_crit_vals - from diptest.consts.Consts._CRIT_VALS
    !> (RUrlus/diptest, stable branch). 21 sample sizes x 26 alpha
    !> quantiles; row i is sqrt(sample_size(i))-scale-ready as-is (the
    !> sqrt(n) scaling happens below). Not hand-edited: regenerate by
    !> reading the same table out of the diptest PyPI package.
    data (crit_vals(i,1), i=1,26) / 1.2500000000000000D-01, 1.2500000000000000D-01, 1.2500000000000000D-01, &
        1.2500000000000000D-01, 1.2500000000000000D-01, 1.2500000000000000D-01, &
        1.2500000000000000D-01, 1.2500000000000000D-01, 1.2500000000000000D-01, &
        1.2500000000000000D-01, 1.3255954878268900D-01, 1.5749736904023501D-01, &
        1.8740187880755901D-01, 2.0726978858736000D-01, 2.2375580462922201D-01, &
        2.3179625886419200D-01, 2.3726374382677901D-01, 2.4199289268859300D-01, &
        2.4436983904963200D-01, 2.4596662550469101D-01, 2.4743959723326200D-01, &
        2.4823065965663799D-01, 2.4875426914641599D-01, 2.4930203997425901D-01, &
        2.4945965232322501D-01, 2.4974836247845000D-01 /
    data (crit_vals(i,2), i=1,26) / 1.0000000000000001D-01, 1.0000000000000001D-01, 1.0000000000000001D-01, &
        1.0000000000000001D-01, 1.0000000000000001D-01, 1.0000000000000001D-01, &
        1.0000000000000001D-01, 1.0872059357632900D-01, 1.2156379802641400D-01, &
        1.3431891869705301D-01, 1.4729879897625200D-01, 1.6108502570260400D-01, &
        1.7681199847607601D-01, 1.8639179602794401D-01, 1.9361253363045000D-01, &
        1.9650913979884499D-01, 1.9815996728757601D-01, 1.9924427936243300D-01, &
        1.9961752740616601D-01, 1.9980094149902800D-01, 1.9991708183427101D-01, &
        1.9995902909307500D-01, 1.9997839537608200D-01, 1.9999315140581500D-01, &
        1.9999552502567300D-01, 1.9999983563921100D-01 /
    data (crit_vals(i,3), i=1,26) / 8.3333333333333301D-02, 8.3333333333333301D-02, 8.3333333333333301D-02, &
        8.3333333333333301D-02, 8.3333333333333301D-02, 9.2451447094193298D-02, &
        1.0391343105994900D-01, 1.1388522064021200D-01, 1.2307118713778099D-01, &
        1.3186973390253001D-01, 1.4056479649794101D-01, 1.4941924112912999D-01, &
        1.5913706457262700D-01, 1.6476960851330200D-01, 1.7917654739278199D-01, &
        1.9186282799556301D-01, 2.0210197104296801D-01, 2.1301578111118599D-01, &
        2.1951862728241500D-01, 2.2433904739444599D-01, 2.2944933215424099D-01, &
        2.3271453044960200D-01, 2.3654812835896899D-01, 2.3908879119949999D-01, &
        2.4010356643629499D-01, 2.4467288361776801D-01 /
    data (crit_vals(i,4), i=1,26) / 7.1428571428571397D-02, 7.1428571428571397D-02, 7.1428571428571397D-02, &
        7.2571781625074203D-02, 8.1731547853948899D-02, 9.4059018192252694D-02, &
        1.0324449080087100D-01, 1.1096459999569699D-01, 1.1780784650433500D-01, &
        1.2421608683353100D-01, 1.3040901396831700D-01, 1.3663964212306801D-01, &
        1.4424066903512400D-01, 1.5990339567833600D-01, 1.7519655327122299D-01, &
        1.8411865912150099D-01, 1.9101439617430599D-01, 1.9821679523218200D-01, &
        2.0234101074826100D-01, 2.0537756634683199D-01, 2.0830656252687399D-01, &
        2.0986604785237900D-01, 2.1096757693345100D-01, 2.1223334855870199D-01, &
        2.1266103831250599D-01, 2.1353618608816999D-01 /
    data (crit_vals(i,5), i=1,26) / 6.2500000000000000D-02, 6.2500000000000000D-02, 6.5691199450328294D-02, &
        7.3865113607176194D-02, 8.2004591776251204D-02, 9.2270060113189195D-02, &
        9.9673718959936305D-02, 1.0573353180273699D-01, 1.1103512984770500D-01, &
        1.1592005574998800D-01, 1.2056147926246499D-01, 1.2555875903484501D-01, &
        1.4184106703389901D-01, 1.5397830399856099D-01, 1.6597856724751001D-01, &
        1.7298852827675901D-01, 1.7901041349637400D-01, 1.8650438871117800D-01, &
        1.9448404115793999D-01, 2.0086429700502600D-01, 2.0884999705022900D-01, &
        2.1255604040621900D-01, 2.1714917413729901D-01, 2.2170007640450301D-01, &
        2.2500083535753199D-01, 2.3377291968768299D-01 /
    data (crit_vals(i,6), i=1,26) / 5.5555555555555601D-02, 6.1301809029892400D-02, 6.5861585817931501D-02, &
        7.3265114253531702D-02, 8.0394162959347495D-02, 8.9043242091384797D-02, &
        9.5081142029792801D-02, 9.9938089781104605D-02, 1.0415356007586800D-01, &
        1.0800780236193200D-01, 1.1251261712495100D-01, 1.2291503348081700D-01, &
        1.3641263938708401D-01, 1.4660378495401899D-01, 1.5708406565316599D-01, &
        1.6416464365721700D-01, 1.7282167458233799D-01, 1.8255528356781800D-01, &
        1.8865883312190601D-01, 1.9408912076824600D-01, 1.9915700809389000D-01, &
        2.0288159843655801D-01, 2.0597979573512901D-01, 2.1054115498897999D-01, &
        2.1180033095039000D-01, 2.1537991431762499D-01 /
    data (crit_vals(i,7), i=1,26) / 5.0000000000000003D-02, 6.1013255562326903D-02, 6.5162733321401600D-02, &
        7.1832161965616495D-02, 7.7966212182458999D-02, 8.5283535983456393D-02, &
        9.0320417370709893D-02, 9.4333498374511701D-02, 9.7781763038472497D-02, &
        1.0218086669662800D-01, 1.0996094814295100D-01, 1.1884476721158700D-01, &
        1.3046214964481900D-01, 1.3961139513709900D-01, 1.5096172888248099D-01, &
        1.5968415885823500D-01, 1.6719524735673999D-01, 1.7541954085608200D-01, &
        1.8061119579735099D-01, 1.8528641605039600D-01, 1.9120308390504401D-01, &
        1.9580515933918399D-01, 2.0029398089673001D-01, 2.0565108964621900D-01, &
        2.0968204878585300D-01, 2.2153028218296300D-01 /
    data (crit_vals(i,8), i=1,26) / 3.4137817227791897D-02, 5.4628420804897500D-02, 5.7219126023181500D-02, &
        6.1008736768969202D-02, 6.4265713733044405D-02, 6.9223410798959106D-02, &
        7.4546211436516699D-02, 7.9203087898176205D-02, 8.3621033469191003D-02, &
        8.8119848220290495D-02, 9.3124666680253002D-02, 9.9669439339068897D-02, &
        1.1008749690090600D-01, 1.1876076920366400D-01, 1.2889047521005501D-01, &
        1.3598356863635999D-01, 1.4245248368127700D-01, 1.5017281653074199D-01, &
        1.5545613369632799D-01, 1.6089649910695800D-01, 1.6697940794624799D-01, &
        1.7111793515550999D-01, 1.7590050570443200D-01, 1.8185667601316599D-01, &
        1.8574345415100399D-01, 1.9224056333056200D-01 /
    data (crit_vals(i,9), i=1,26) / 3.3718563622064997D-02, 4.7433374069840099D-02, 4.9089138762709199D-02, &
        5.2719998201553001D-02, 5.6779550905674200D-02, 6.2013467446818099D-02, &
        6.6016387206904795D-02, 6.9650607506640094D-02, 7.3343774059271394D-02, &
        7.7646066288025395D-02, 8.2455840711837203D-02, 8.8344627001736994D-02, &
        9.7234601812290294D-02, 1.0513021827063600D-01, 1.1430970428125301D-01, &
        1.2062404333582100D-01, 1.2655237803673899D-01, 1.3360135382395000D-01, &
        1.3856990379176701D-01, 1.4336916123967999D-01, 1.4894011639488300D-01, &
        1.5283253818362200D-01, 1.5601016361897099D-01, 1.6131922583934499D-01, &
        1.6556825591674901D-01, 1.7583445952278901D-01 /
    data (crit_vals(i,10), i=1,26) / 2.6267448507564201D-02, 3.9587189040574899D-02, 4.1457460674167300D-02, &
        4.4446261406995598D-02, 4.7399852504268598D-02, 5.1667737037434901D-02, &
        5.5103751900162201D-02, 5.8265005347493001D-02, 6.1451085730434299D-02, &
        6.4916440805397796D-02, 6.8917876242544196D-02, 7.3924907407829102D-02, &
        8.1479137939012694D-02, 8.8168914312666602D-02, 9.6056438301364400D-02, &
        1.0147855889383700D-01, 1.0650487144103001D-01, 1.1272463652426200D-01, &
        1.1716414018441700D-01, 1.2142585990898699D-01, 1.2673305188940101D-01, &
        1.3119857889754200D-01, 1.3369173948344401D-01, 1.3783163795069400D-01, &
        1.4155750962435101D-01, 1.6383304605981699D-01 /
    data (crit_vals(i,11), i=1,26) / 2.1854478136454501D-02, 3.1440050199991597D-02, 3.2900816047083399D-02, &
        3.5302381904001597D-02, 3.7727997310248201D-02, 4.1069998439958198D-02, &
        4.3770459862266499D-02, 4.6292564267129903D-02, 4.8851155289607998D-02, &
        5.1614589786575703D-02, 5.4812193206601897D-02, 5.8823048285136598D-02, &
        6.4913632404676694D-02, 7.0273787719126901D-02, 7.6709588607917906D-02, &
        8.1199841535591802D-02, 8.5285464666213395D-02, 9.0484782749029394D-02, &
        9.4093010666624399D-02, 9.7490434491674299D-02, 1.0228420428399700D-01, &
        1.0468062433461101D-01, 1.0749669423503901D-01, 1.1140887547015001D-01, &
        1.1353660771741100D-01, 1.1788671686531201D-01 /
    data (crit_vals(i,12), i=1,26) / 1.6485259743840300D-02, 2.2831985803042999D-02, 2.3891748644284901D-02, &
        2.5655953797757900D-02, 2.7398741457094800D-02, 2.9810937083015299D-02, &
        3.1777149653025298D-02, 3.3607382159038697D-02, 3.5462176059211301D-02, &
        3.7480584455027201D-02, 3.9804617911659901D-02, 4.2728384679916603D-02, &
        4.7152783315717997D-02, 5.1127944286882700D-02, 5.5802205219520798D-02, &
        5.9024132304225999D-02, 6.2042506516514599D-02, 6.5801601146609906D-02, &
        6.8447973111802798D-02, 7.0916944399419299D-02, 7.4118348608126300D-02, &
        7.6257940290383797D-02, 7.8573596793497902D-02, 8.1345835688913307D-02, &
        8.3296301375552204D-02, 9.2678042309673705D-02 /
    data (crit_vals(i,13), i=1,26) / 1.1123638884968800D-02, 1.6501773542982500D-02, 1.7259415799248900D-02, &
        1.8525942603292600D-02, 1.9791761263752101D-02, 2.1523374577845401D-02, &
        2.2925976987042799D-02, 2.4243848341112002D-02, 2.5584358256487000D-02, &
        2.7025212981628799D-02, 2.8692026215051701D-02, 3.0800676634140600D-02, &
        3.3996781429350399D-02, 3.6841841387830698D-02, 4.0272985031639702D-02, &
        4.2686479977744801D-02, 4.4958959158760997D-02, 4.7764387374944900D-02, &
        4.9719800186743698D-02, 5.1611461180145098D-02, 5.4054397886465197D-02, &
        5.5870452618263802D-02, 5.7387705633022798D-02, 5.9336590165387802D-02, &
        6.0764631047391097D-02, 7.0530910788239504D-02 /
    data (crit_vals(i,14), i=1,26) / 7.5548859757619598D-03, 1.0640346112751499D-02, 1.1125557320829401D-02, &
        1.1935365532893099D-02, 1.2741130641180799D-02, 1.3852454275181400D-02, &
        1.4753600428847600D-02, 1.5596318575104800D-02, 1.6451923802528599D-02, &
        1.7383057902553001D-02, 1.8450394988773499D-02, 1.9816267978207101D-02, &
        2.1878131318220299D-02, 2.3729474263341099D-02, 2.5919578977656999D-02, &
        2.7451802276199699D-02, 2.8898636956430100D-02, 3.0681350505016299D-02, &
        3.2017099682318903D-02, 3.3245274733295901D-02, 3.4833569857616799D-02, &
        3.5983238931746098D-02, 3.6905199584064498D-02, 3.8722115925642397D-02, &
        3.9930259057649999D-02, 4.3144816361717797D-02 /
    data (crit_vals(i,15), i=1,26) / 5.4165812787212199D-03, 7.6028674530018697D-03, 7.9498783464479906D-03, &
        8.5216518343593992D-03, 9.0977560553325305D-03, 9.8892452101407794D-03, &
        1.0530929709048200D-02, 1.1132272679738400D-02, 1.1743900905255200D-02, &
        1.2405033293814000D-02, 1.3168417932080300D-02, 1.4137794260304700D-02, &
        1.5614805502305800D-02, 1.6934397006756401D-02, 1.8513067368104000D-02, &
        1.9608026048323401D-02, 2.0648956858736401D-02, 2.1928517676508202D-02, &
        2.2868916897266899D-02, 2.3738710122235000D-02, 2.4833415889143201D-02, &
        2.5612657343359602D-02, 2.6549133693682898D-02, 2.7578430100536001D-02, &
        2.8443073310800000D-02, 3.1364094198210797D-02 /
    data (crit_vals(i,16), i=1,26) / 3.9043999745055698D-03, 5.4166418179658303D-03, 5.6617138625232304D-03, &
        6.0712097113522897D-03, 6.4762535755248001D-03, 7.0357309859002898D-03, &
        7.4942125458929898D-03, 7.9208788960173308D-03, 8.3557372476800607D-03, &
        8.8243933381235099D-03, 9.3678582071706103D-03, 1.0056046038840000D-02, &
        1.1101911683759100D-02, 1.2038099032834100D-02, 1.3172101055257599D-02, &
        1.3965512228196900D-02, 1.4688912220448800D-02, 1.5607677964745400D-02, &
        1.6268561599624799D-02, 1.6887493778941502D-02, 1.7650509338815300D-02, &
        1.8194426540050400D-02, 1.8622603781852300D-02, 1.9300179656543300D-02, &
        1.9624151804061699D-02, 2.1308125407458401D-02 /
    data (crit_vals(i,17), i=1,26) / 2.4565778544043300D-03, 3.4480928223332599D-03, 3.6047394371303602D-03, &
        3.8632654801084901D-03, 4.1208950675269201D-03, 4.4764005013747899D-03, &
        4.7655569310227604D-03, 5.0370402975007198D-03, 5.3123924740821303D-03, &
        5.6092991935995902D-03, 5.9535272837794896D-03, 6.3909228056351700D-03, &
        7.0556612623462502D-03, 7.6506368153934998D-03, 8.3682168704721505D-03, &
        8.8635789285491408D-03, 9.3416278718615898D-03, 9.9321863632402894D-03, &
        1.0349879529162900D-02, 1.0778090707686200D-02, 1.1318431686828299D-02, &
        1.1732944646857099D-02, 1.1999594896837501D-02, 1.2441005202788600D-02, &
        1.2946739673312800D-02, 1.4396063834027001D-02 /
    data (crit_vals(i,18), i=1,26) / 1.7495426919956600D-03, 2.4459513388530199D-03, 2.5571080227561201D-03, &
        2.7399095522726499D-03, 2.9225480567908000D-03, 3.1737463842246498D-03, &
        3.3807225853352699D-03, 3.5724387653598201D-03, 3.7673471575220899D-03, &
        3.9788500724913202D-03, 4.2243001317623296D-03, 4.5343750814854202D-03, &
        5.0017880840236796D-03, 5.4237224283639500D-03, 5.9265668102285902D-03, &
        6.2803473288037398D-03, 6.6103064155087297D-03, 7.0225469996764798D-03, &
        7.3182262815645804D-03, 7.6065423418208000D-03, 7.9564036720748202D-03, &
        8.2270524584353993D-03, 8.5224098978625099D-03, 8.9286390554030298D-03, &
        9.1385393300021309D-03, 9.5223457956677294D-03 /
    data (crit_vals(i,19), i=1,26) / 1.1945881410609100D-03, 1.7343534689628699D-03, 1.8119443458468100D-03, &
        1.9425947048589301D-03, 2.0717371962386800D-03, 2.2499320208695501D-03, &
        2.3952083147341899D-03, 2.5303679282466501D-03, 2.6686316871811400D-03, &
        2.8181999035215999D-03, 2.9913754814207701D-03, 3.2102489992013499D-03, &
        3.5436222031415502D-03, 3.8433019024467900D-03, 4.2025879937825300D-03, &
        4.4577490215571098D-03, 4.6946151321274297D-03, 4.9941606912916802D-03, &
        5.2091775774321799D-03, 5.4039623592437198D-03, 5.6454020170459401D-03, &
        5.8046079229921400D-03, 5.9977473959315101D-03, 6.3309925437811396D-03, &
        6.5698710938676200D-03, 6.8582944804622698D-03 /
    data (crit_vals(i,20), i=1,26) / 8.5241564801177697D-04, 1.2288347931066501D-03, 1.2846930445701800D-03, &
        1.3761765052555299D-03, 1.4675150200632299D-03, 1.5937645367246601D-03, &
        1.6966844550615099D-03, 1.7925341833790599D-03, 1.8906126163597699D-03, &
        1.9964547188617899D-03, 2.1192974838170398D-03, 2.2745769870358098D-03, &
        2.5099908089039700D-03, 2.7237507348622301D-03, 2.9807295856838700D-03, &
        3.1594219404038801D-03, 3.3273652798147999D-03, 3.5398896569857900D-03, &
        3.6940004548662499D-03, 3.8334571537218202D-03, 4.0079346963469596D-03, &
        4.1489273722288503D-03, 4.2839159079760998D-03, 4.4187010443287903D-03, &
        4.5081860456917897D-03, 5.1347746756558298D-03 /
    data (crit_vals(i,21), i=1,26) / 6.4440005325699695D-04, 9.1687220448428302D-04, 9.5793294676553196D-04, &
        1.0264186387234700D-03, 1.0949515421800201D-03, 1.1890409036941500D-03, &
        1.2657519769987400D-03, 1.3375096636150600D-03, 1.4104970922847199D-03, &
        1.4893670929880200D-03, 1.5802754194562600D-03, 1.6965164386007401D-03, &
        1.8730618472582599D-03, 2.0317840161055501D-03, 2.2235609750605400D-03, &
        2.3578281477762701D-03, 2.4834358012706700D-03, 2.6421082633949801D-03, &
        2.7524322157581002D-03, 2.8608570740143000D-03, 2.9869504450800301D-03, &
        3.0934009203805899D-03, 3.1993276719880100D-03, 3.3268823461118698D-03, &
        3.3931609447735499D-03, 3.7633169700585899D-03 /
    !> END GENERATED dip_crit_vals

    i1 = n_tab_sizes + 1
    do k = 1, n_tab_sizes
        if (sample_size(k) >= n) then
            i1 = k
            exit
        end if
    end do
    i0 = i1 - 1
    if (i0 < 1) i0 = 1
    if (i1 > n_tab_sizes) i1 = n_tab_sizes

    n0 = sample_size(i0)
    n1 = sample_size(i1)

    do i = 1, n_tab_alphas
        y0(i) = dsqrt(dble(n0)) * crit_vals(i, i0)
    end do
    sD = dsqrt(dble(n)) * dip

    if (i0 == i1) then
        xp = y0
    else
        frac_n = dble(n - n0) / dble(n1 - n0)
        do i = 1, n_tab_alphas
            y1(i) = dsqrt(dble(n1)) * crit_vals(i, i1)
        end do
        xp = y0 + frac_n * (y1 - y0)
    end if

    pval = 1d0 - DipLinInterp1D(sD, xp, alpha, n_tab_alphas)
end subroutine DipPvalue
double precision function DipLinInterp1D(x0, xp, fp, m) result(y)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: m
    real(kind = dbl), intent(in) :: x0, xp(m), fp(m)
    integer :: k

    if (x0 < xp(1)) then
        y = fp(1)
    else if (x0 > xp(m)) then
        y = fp(m)
    else
        !> Search from the right: an x0 that exactly matches a run of
        !> duplicate xp nodes (common near the table's small-n end, where
        !> several alpha quantiles share one critical value) resolves to
        !> the LAST matching node - matching np.interp's actual behaviour,
        !> confirmed by testing rather than assumed from its docs, which
        !> do not specify tie-breaking. Strict "<"/">"" above is what lets
        !> an x0 sitting exactly on xp(1) or xp(m) still reach this search
        !> when that boundary value is itself part of a duplicate run.
        y = fp(1)
        do k = m - 1, 1, -1
            if (x0 >= xp(k) .and. x0 <= xp(k+1)) then
                y = fp(k) + (x0 - xp(k)) / (xp(k+1) - xp(k)) * (fp(k+1) - fp(k))
                exit
            end if
        end do
    end if
end function DipLinInterp1D
