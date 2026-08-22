!***************************************************************************
! bpcf_cospectral_models.f90
! --------------------------
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
! \brief       Calculates theoretical co-spectra, after \n
!              Moncrieff et al. (1997), Eq 13-16
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CospectraMoncrieff97(nf, kf, Cospectrum, zL, N)
    use m_common_global_var
    implicit none
    ! in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: kf(N)
    real(kind = dbl), intent(in) :: nf(N)
    real(kind = dbl), intent(in) :: zL
    type(SpectralType), intent(out)  :: Cospectrum(N)
    ! local variables
    integer :: i = 0
    integer :: gas
    real(kind = dbl) :: Ac
    real(kind = dbl) :: Au
    real(kind = dbl) :: Bc
    real(kind = dbl) :: Bu


    if(zL > 0.d0) then
        do i = 1, N
            !> 1.1 stable conditions
            !> Scalar fluxes (T, CO2, H2O)
            Ac = 0.284d0 * ((1.d0 + 6.4d0 * zL)**0.75d0)
            Bc = 2.34d0 * (Ac**(-1.1d0))
            Cospectrum(i)%of(firstGas) = kf(i) / (nf(i) * (Ac + Bc * kf(i)**2.1d0))
            ! Reynolds stress, however not corrected so far
            Au = 0.124d0 * ((1.d0 + 7.9d0 * zL)**0.75d0)
            Bu = 2.34d0 * (Au**(-1.1d0))
            Cospectrum(i)%of(w_u) = kf(i) / (nf(i) * (Au + Bu * kf(i)**2.1d0))
        end do
    else
        do i = 1, N
            !> 1.2 unstable conditions
            !> Scalar fluxes (T, CO2, H2O)
            if(kf(i) < 0.54d0) then
                Cospectrum(i)%of(firstGas) = 12.92d0 * kf(i) / (nf(i) * (1.d0 + 26.7d0 * kf(i))**1.375d0)
            else
                Cospectrum(i)%of(firstGas) = 4.378d0 * kf(i) / (nf(i) * (1.d0 + 3.8d0 * kf(i))**2.4d0)
            end if
            ! reynolds stress, however not corrected so far
            if(kf(i) < 0.24d0) then
                Cospectrum(i)%of(w_u) = 20.78d0 * kf(i) / (nf(i) * (1.d0 + 31.0d0 * kf(i))**1.575d0)
            else
                Cospectrum(i)%of(w_u) = 12.66d0 * kf(i) / (nf(i) * (1.d0 + 9.6d0 * kf(i))**2.4d0)
            end if
        end do
    end if
    !> Set the cospectral model for temperature and every gas equal to the
    !> first scalar slot, which is where the formulae above deposit it.
    !>
    !> Moncrieff's scalar cospectrum is the same curve for temperature and
    !> every scalar - the formulae carry no species term - so which slot it
    !> is computed in is arbitrary and the rest are copies.
    !>
    !> This used to name h2o, ch4 and gas4 one by one, so a gas past the fourth
    !> kept the error value the array is initialised to. That made its
    !> cospectrum "all error", which is how SpectralCorrectionFactors decides a
    !> correction factor is unavailable - so every analytic method silently
    !> declined to correct those gases, without the cospectra file being
    !> involved at all.
    Cospectrum(:)%of(w_ts) = Cospectrum(:)%of(firstGas)
    do gas = firstGas, lastGas
        Cospectrum(:)%of(gas) = Cospectrum(:)%of(firstGas)
    end do
end subroutine CospectraMoncrieff97
!***************************************************************************
!
! \brief       Fills the model cospectrum with whichever analytic shape the
!              project asked for, and always leaves the Reynolds stress on
!              Moncrieff's.
! \author      Jonathan Muller
! \note
!              WHY THIS EXISTS. Every analytic spectral correction in this
!              program is a ratio of two integrals of the SAME model
!              cospectrum - one filtered by the transfer function, one not
!              (SpectralCorrectionFactors, bpcf_aux_subs.f90). The shape
!              therefore decides how much of the flux sits at the frequencies
!              the instrument attenuates, and so the size of the correction.
!              Until now that shape was one curve, Moncrieff et al. (1997),
!              with no way to ask what a different one would give.
!
!              WHAT CANCELS. The correction is IntCO / IntTFCO, so any
!              constant multiplying the whole cospectrum divides out exactly.
!              Only the SHAPE matters. That is why Kristensen's 2*pi/z
!              prefactor is dropped here, and why his numerically integrated
!              normalisation - which in EddyUH depends on hard-coded indices
!              into the frequency vector it is handed - is not reproduced. It
!              could not change a correction factor if it were.
!
!              REYNOLDS STRESS. Sakai, Su, Moraes and Kristensen are scalar
!              cospectra; they have no momentum form. Rather than pair a
!              chosen scalar model with an unrelated momentum one for some
!              options and not others, Cospectrum%of(w_u) stays Moncrieff's
!              for every option, uniformly. It is what BPCF_AnemometricFluxes
!              integrates for tau, and nothing else reads it.
!
!              STABILITY. Only Moncrieff and Kaimal have stable and unstable
!              branches. The other four are single neutral-to-unstable forms
!              that EddyUH applies at any zL, and this does the same. Under
!              stable stratification the cospectral peak moves to higher
!              frequency, so a neutral-form model puts too little flux there
!              and will UNDERSTATE the high-frequency loss. Said in the
!              interface tooltip too.
!
!              PROVENANCE. The formulae are EddyUH's
!              EC_Software_Spectral_Analysis/modelcospectra.m, which uses them
!              for diagnostic plots and not for its own correction - EddyUH
!              corrects against measured and fitted cospectra instead. So this
!              is EddyUH's model library, not EddyUH's method.
!
!              Worth knowing: EddyUH's CMoore, after Moore (1986), is
!              term-for-term the curve this program calls Moncrieff et al.
!              (1997), scalar and momentum both. The two names are the same
!              cospectrum. Kaimal et al. (1972) below is a genuinely different
!              curve - and a different one again from the column the spectral
!              assessment output labels kaimal_cosp, which is computed by
!              kaimal() in kaimal_models.f90 and is also Moore's.
! \sa          CospectraMoncrieff97, SpectralCorrectionFactors
!***************************************************************************
subroutine CospectralModel(nf, kf, Cospectrum, zL, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: kf(N)
    real(kind = dbl), intent(in) :: nf(N)
    real(kind = dbl), intent(in) :: zL
    type(SpectralType), intent(out)  :: Cospectrum(N)
    !> local variables
    integer :: i
    integer :: gas
    real(kind = dbl) :: co
    real(kind = dbl) :: hs, n0, x

    !> Always. It fills the momentum slot, which nothing below touches, and
    !> leaves the default path bit-for-bit what it was before this routine
    !> existed: the select below has no branch for moncrieff_97.
    call CospectraMoncrieff97(nf, kf, Cospectrum, zL, N)

    select case (trim(adjustl(EddyFlowProj%cosp_model)))
        case ('kaimal_72')
            !> Kaimal et al. (1972), Kansas. Stable form peak-normalised about
            !> n0, the intercept of the extended -4/3 slope with fCo/cov = 1.
            if (zL > 0d0) then
                hs = 1d0 + 6.4d0 * zL
                n0 = 0.23d0 * hs**0.75d0
                do i = 1, N
                    x = kf(i) / n0
                    co = 0.88d0 * x / (1d0 + 1.5d0 * x**2.1d0)
                    call PutScalar(i, co / nf(i))
                end do
            else
                do i = 1, N
                    if (kf(i) < 1d0) then
                        co = 11d0 * kf(i) / ((1d0 + 13.3d0 * kf(i))**1.75d0)
                    else
                        co = 4d0 * kf(i) / ((1d0 + 3.8d0 * kf(i))**(7d0/3d0))
                    end if
                    call PutScalar(i, co / nf(i))
                end do
            end if

        case ('sakai_01')
            !> Sakai et al. (2001), over a rough surface - more flux at low
            !> frequency than Kaimal, which is the point of the paper.
            do i = 1, N
                co = 8d0 * kf(i) / (1d0 + (11d0 * kf(i))**(7d0/3d0))
                call PutScalar(i, co / nf(i))
            end do

        case ('su_03')
            !> Su et al. (2003), two mixed hardwood forests in non-flat terrain.
            do i = 1, N
                co = 10.4d0 * kf(i) / ((1d0 + 10.6d0 * kf(i))**2d0)
                call PutScalar(i, co / nf(i))
            end do

        case ('moraes_08')
            !> Moraes et al. (2008).
            do i = 1, N
                co = 4.1d0 * kf(i) / (1d0 + (6.6d0 * kf(i))**(7d0/3d0))
                call PutScalar(i, co / nf(i))
            end do

        case ('kristensen_97')
            !> Kristensen et al. (1997). Their 2*pi/z prefactor is a constant
            !> and cancels in the correction ratio, so it is not carried.
            !> The exponents are written as 2*mu and 7/(6*mu) rather than
            !> folded to 0.46 and 5.0725, so that mu = 0.23 stays visible.
            do i = 1, N
                co = kf(i) / ((1d0 + (1.2d0 * 2d0 * p * kf(i)) &
                    **(2d0 * 0.23d0))**(7d0 / (6d0 * 0.23d0)))
                call PutScalar(i, co / nf(i))
            end do
    end select

contains

    !> Temperature and every gas take the same curve - none of these models
    !> carries a species term. Writes the whole scalar range rather than
    !> naming slots, so a fifth gas is not left holding the error value the
    !> array is initialised to, which is how SpectralCorrectionFactors decides
    !> a correction is unavailable.
    subroutine PutScalar(k, value)
        integer, intent(in) :: k
        real(kind = dbl), intent(in) :: value

        Cospectrum(k)%of(w_ts) = value
        do gas = firstGas, lastGas
            Cospectrum(k)%of(gas) = value
        end do
    end subroutine PutScalar

end subroutine CospectralModel
