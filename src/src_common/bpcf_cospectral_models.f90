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
