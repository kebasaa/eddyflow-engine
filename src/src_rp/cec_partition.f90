!***************************************************************************
! cec_partition.f90
! -----------------
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
! \brief   Conditional Eddy Covariance partitioning (Zahn et al. 2022)
!          Partitions H2O and CO2 fluxes into stomatal/non-stomatal
!          components using octant-based conditional statistics.
!
! Method:
!   Octant O1 (non-stomatal, I_E): w'>0, q'>0, c'>0
!   Octant O2 (stomatal, I_T):     w'>0, q'>0, c'<0
!
!   Sample fluxes (N = total finite points):
!     f_E = sum(w'*q' | I_E) / N    f_T = sum(w'*q' | I_T) / N
!     f_R = sum(w'*c' | I_E) / N    f_P = sum(w'*c' | I_T) / N
!     r_ET = f_E / f_T              r_Fc = f_R / f_P
!
!   Partitioned fluxes (WPL-corrected totals):
!     E = ET / (1 + 1/r_ET)         T = ET / (1 + r_ET)
!     R = Fc / (1 + 1/r_Fc)         P = Fc / (1 + r_Fc)
!
! Reference: Zahn et al. (2022), doi:10.1111/gcb.16122
!***************************************************************************
subroutine CecFluxes(primes, ET_total, Fc_total, do_cec)
    use m_rp_global_var
    implicit none
    !> Detrended fluctuations array (rows=samples, cols=variables).
    real(kind = dbl), intent(in) :: primes(:, :)
    real(kind = dbl), intent(in) :: ET_total  !< WPL-corrected ET [mmol m-2 s-1]
    real(kind = dbl), intent(in) :: Fc_total  !< WPL-corrected NEE [umol m-2 s-1]
    integer, intent(in) :: do_cec             !< 1=H2O+CO2, 2=H2O only, 3=CO2 only
    !> local variables
    integer :: i, nr
    integer :: N, n_O1, n_O2
    real(kind = dbl) :: wpr, qpr, cpr
    real(kind = dbl) :: sum_fE, sum_fT, sum_fR, sum_fP
    real(kind = dbl) :: f_E, f_T, f_R, f_P
    real(kind = dbl) :: r_ET, r_Fc
    real(kind = dbl) :: frac_O1, frac_O2

    nr = size(primes, 1)

    !> Initialise output to error
    CECFlux%E_cec    = error
    CECFlux%Tr_cec   = error
    CECFlux%Reco_cec = error
    CECFlux%GPP_cec  = error
    CECFlux%NEE_cec  = error
    CECFlux%r_ET_cec = error
    CECFlux%r_Fc_cec = error
    CECFlux%ok       = .false.

    !> ===== Phase 1: octant conditional sample fluxes =====
    N = 0; n_O1 = 0; n_O2 = 0
    sum_fE = 0d0; sum_fT = 0d0
    sum_fR = 0d0; sum_fP = 0d0

    do i = 1, nr
        wpr = primes(i, gW)
        qpr = primes(i, h2o)
        cpr = primes(i, co2)
        if (wpr == error .or. qpr == error .or. cpr == error) cycle
        N = N + 1
        !> Octant O1: non-stomatal (evaporation/respiration)
        if (wpr > 0d0 .and. qpr > 0d0 .and. cpr > 0d0) then
            n_O1   = n_O1 + 1
            sum_fE = sum_fE + wpr * qpr
            sum_fR = sum_fR + wpr * cpr
        end if
        !> Octant O2: stomatal (transpiration/photosynthesis)
        if (wpr > 0d0 .and. qpr > 0d0 .and. cpr < 0d0) then
            n_O2   = n_O2 + 1
            sum_fT = sum_fT + wpr * qpr
            sum_fP = sum_fP + wpr * cpr
        end if
    end do

    if (N < 2) return

    f_E = sum_fE / dble(N)
    f_T = sum_fT / dble(N)
    f_R = sum_fR / dble(N)
    f_P = sum_fP / dble(N)

    !> Flux ratios (undefined if denominator is zero)
    r_ET = error
    r_Fc = error
    if (f_T /= 0d0) r_ET = f_E / f_T
    if (f_P /= 0d0) r_Fc = f_R / f_P

    CECFlux%r_ET_cec = r_ET
    CECFlux%r_Fc_cec = r_Fc

    !> ===== Phase 2: octant validity and flux partitioning =====
    frac_O1 = dble(n_O1) / dble(N)
    frac_O2 = dble(n_O2) / dble(N)

    !> Reject period if combined octant fraction is too small
    if ((frac_O1 + frac_O2) < 0.20d0) return

    !> H2O partitioning (do_cec = 1 or 2)
    if ((do_cec == 1 .or. do_cec == 2) .and. ET_total /= error) then
        if (frac_O1 < 0.05d0) then
            !> Almost no non-stomatal signal -> all transpiration
            CECFlux%E_cec  = 0d0
            CECFlux%Tr_cec = ET_total
        else if (frac_O2 < 0.05d0) then
            !> Almost no stomatal signal -> all evaporation
            CECFlux%E_cec  = ET_total
            CECFlux%Tr_cec = 0d0
        else if (r_ET /= error .and. r_ET /= 0d0) then
            CECFlux%E_cec  = ET_total / (1d0 + 1d0 / r_ET)
            CECFlux%Tr_cec = ET_total / (1d0 + r_ET)
        else
            CECFlux%Tr_cec = ET_total
            CECFlux%E_cec  = 0d0
        end if
    end if

    !> CO2 partitioning (do_cec = 1 or 3)
    if ((do_cec == 1 .or. do_cec == 3) .and. Fc_total /= error) then
        !> Near-singularity: r_Fc ~ -1 means Reco and GPP nearly cancel
        if (r_Fc /= error .and. abs(r_Fc + 1d0) < 0.05d0) then
            !> Reco_cec and GPP_cec remain error (undefined)
        else if (frac_O1 < 0.05d0) then
            CECFlux%Reco_cec = 0d0
            CECFlux%GPP_cec  = Fc_total
        else if (frac_O2 < 0.05d0) then
            CECFlux%Reco_cec = Fc_total
            CECFlux%GPP_cec  = 0d0
        else if (r_Fc /= error .and. r_Fc /= 0d0) then
            CECFlux%Reco_cec = Fc_total / (1d0 + 1d0 / r_Fc)
            CECFlux%GPP_cec  = Fc_total / (1d0 + r_Fc)
        else
            CECFlux%Reco_cec = Fc_total
            CECFlux%GPP_cec  = 0d0
        end if
        !> NEE_cec = Reco + GPP (GPP is negative by sign convention -> recovers Fc_total)
        if (CECFlux%Reco_cec /= error .and. CECFlux%GPP_cec /= error) then
            CECFlux%NEE_cec = CECFlux%Reco_cec + CECFlux%GPP_cec
        end if
    end if

    CECFlux%ok = .true.

end subroutine CecFluxes

!***************************************************************************
!
! \brief   Apply pre-computed CEC ratios (r_ET, r_Fc) to total fluxes.
!          Called by FCC after reading the CEC ratios file written by RP.
!          Populates the module-level CECFlux fields using the same
!          validity guards and formulas as phase 2 of CecFluxes.
!
! \param   ET_total  WPL-corrected ET [mmol m-2 s-1]
! \param   Fc_total  WPL-corrected NEE [umol m-2 s-1]
! \param   r_ET      Sample flux ratio fE/fT from RP (or error)
! \param   r_Fc      Sample flux ratio fR/fP from RP (or error)
! \param   do_cec    1=H2O+CO2, 2=H2O only, 3=CO2 only
!***************************************************************************
subroutine ApplyCecRatios(ET_total, Fc_total, r_ET, r_Fc, do_cec)
    use m_rp_global_var
    implicit none
    real(kind = dbl), intent(in) :: ET_total
    real(kind = dbl), intent(in) :: Fc_total
    real(kind = dbl), intent(in) :: r_ET
    real(kind = dbl), intent(in) :: r_Fc
    integer, intent(in) :: do_cec

    !> Initialise output to error
    CECFlux%E_cec    = error
    CECFlux%Tr_cec   = error
    CECFlux%Reco_cec = error
    CECFlux%GPP_cec  = error
    CECFlux%NEE_cec  = error
    CECFlux%r_ET_cec = r_ET
    CECFlux%r_Fc_cec = r_Fc
    CECFlux%ok       = .false.

    if (r_ET == error .and. r_Fc == error) return

    !> H2O partitioning (do_cec = 1 or 2)
    if ((do_cec == 1 .or. do_cec == 2) .and. ET_total /= error) then
        if (r_ET /= error .and. r_ET /= 0d0) then
            CECFlux%E_cec  = ET_total / (1d0 + 1d0 / r_ET)
            CECFlux%Tr_cec = ET_total / (1d0 + r_ET)
        end if
    end if

    !> CO2 partitioning (do_cec = 1 or 3)
    if ((do_cec == 1 .or. do_cec == 3) .and. Fc_total /= error) then
        if (r_Fc /= error .and. abs(r_Fc + 1d0) < 0.05d0) then
            !> Near-singularity: leave Reco/GPP as error
        else if (r_Fc /= error .and. r_Fc /= 0d0) then
            CECFlux%Reco_cec = Fc_total / (1d0 + 1d0 / r_Fc)
            CECFlux%GPP_cec  = Fc_total / (1d0 + r_Fc)
        end if
        if (CECFlux%Reco_cec /= error .and. CECFlux%GPP_cec /= error) then
            CECFlux%NEE_cec = CECFlux%Reco_cec + CECFlux%GPP_cec
        end if
    end if

    CECFlux%ok = .true.

end subroutine ApplyCecRatios
