!***************************************************************************
! fit_rh_to_cutoff.f90
! --------------------
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
! \brief       Fit exponential curve to actual Fco vs. RH \n
!              after Ibrom et al. (2007, AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FitRh2Fco()
    use m_fx_global_var
    use m_levenberg_marquardt
    implicit none

    !> Local variables
    integer :: RH
    integer :: m
    integer :: cls
    integer :: cls2
    integer :: cnt
    integer :: cnt2
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    integer, parameter :: npar_EXP = 3
    integer :: info, ipvt_EXP(npar_EXP)
    real(kind = dbl) :: EXPPar(npar_EXP)
    real(kind = dbl) :: tol = 1d-04
    real(kind = dbl) :: mean_fc
    integer :: wsl
    integer :: primary
    real(kind = dbl) :: nyquist
    include '../src_common/interfaces.inc'



    !> One fit per hygrometer, each into its own RegPar(wsl, dum)%e1..e3.
    !>
    !> This used to fit the primary alone and leave one set of coefficients for
    !> the whole project, so a second hygrometer's RH dependence was the
    !> primary's however far apart the two read - on CH-LAE, some fourteen
    !> points, more than an RH class. Its own RH-class cut-offs were fitted by
    !> FitTFModels, written to the assessment file and read back, and then
    !> overridden by the primary's exponential; the round trip existed and did
    !> nothing.
    !>
    !> RegPar(wsl, dum) is safe to write: e1..e3 are components of the element,
    !> distinct from the %fc and %Fn that RH class dum holds.
    primary = PrimaryWaterSlot()
    if (primary < firstGas) return

    !> Allocated after the guard, not before it. xFit, yFit and ddum are one
    !> global set shared with FitTFModels, which sizes them to its own
    !> maxval(nlong) under the same "only if not already allocated" test. This
    !> routine claimed them at 10 and then returned here on a project with no
    !> primary water, without freeing them - so the next FitTFModels found them
    !> allocated, skipped its own sizing, and ran off the end of a
    !> ten-element array. Only projects with no water ever reached it, and
    !> those never got this far while a missing binned-spectra directory was
    !> quietly demoting the run to Moncrieff.
    if (.not. allocated(xFit)) allocate(xFit(10))
    if (.not. allocated(yFit)) allocate(yFit(10))
    if (.not. allocated(ddum)) allocate(ddum(10))

    do wsl = firstGas, lastGas
    if (.not. GasSlotIsWater(wsl)) cycle

    !> Preliminary validation of calculated cut-offs for water vapour, against
    !> THIS hygrometer's Nyquist rather than the station's.
    !>
    !> A hygrometer slower than the station resolves nothing above its own
    !> half-rate, so the station's Nyquist accepts fitted cut-offs it cannot
    !> have measured - by exactly the ratio of the two rates. An ex record
    !> written before the analyser block carried a rate leaves GasAcFreq at the
    !> error code, and then the station's rate is the best answer there is.
    nyquist = FCCMetadata%ac_freq / 2d0
    if (FCCMetadata%GasAcFreq(wsl) > 0d0) &
        nyquist = min(FCCMetadata%GasAcFreq(wsl), FCCMetadata%ac_freq) / 2d0
    where (RegPar(wsl, RH10:RH90)%fc > nyquist .or. &
        RegPar(wsl, RH10:RH90)%fc < 0d0) &
        RegPar(wsl, RH10:RH90)%fc = error

    !> This hygrometer's own path type, not the primary's. Which arm a fit
    !> takes is a property of the analyser being fitted; on a site with one
    !> open path and one closed, the second used to be fitted by whichever arm
    !> the first happened to need.
    if (FCCMetadata%GasPathType(wsl) == 'open') then
        !> If the instrument associated to this H2O reading is an open path
        !> fit exponential model analytically, such that the exponential function
        !> provides a constant value, equal to the mean value of f_cutoff among all
        !> RH classes
        write(*, '(a)', advance = 'no') ' Open-path H2O analyser: analytic &
            &fitting of cut-off frequencies vs. RH.. '
        write(ulog, '(a)', advance = 'no') ' Open-path H2O analyser: analytic &
            &fitting of cut-off frequencies vs. RH.. '
        mean_fc = 0d0
        cnt2 = 0
        do RH = RH10, RH90
            if (RegPar(wsl, RH)%fc /= error) then
                mean_fc = mean_fc + RegPar(wsl, RH)%fc
                cnt2 = cnt2 + 1
            end if
        end do
        mean_fc = mean_fc / cnt2
        RegPar(wsl, dum)%e1 = 1d-15
        RegPar(wsl, dum)%e2 = 1d-15
        RegPar(wsl, dum)%e3 = dlog(mean_fc)
    else
        !> Fit exponential model by least squares minimization
        write(*, '(a)', advance = 'no') ' Fitting in-situ assessment of &
            &cut-off frequencies vs. RH.. '
        write(ulog, '(a)', advance = 'no') ' Fitting in-situ assessment of &
            &cut-off frequencies vs. RH.. '

        !> Initialization of function parameters (see Ibrom et al. 2007, AFM)
        !> EXP: Par(1) = A, Par(2)= B, Par(3)=C
        !> in function f(x) = exp**(A * x**2 + B * x + C)
        EXPPar(1) = -2.d0
        EXPPar(2) = -1.d0
        EXPPar(3) = -2.d0
        m = 0
        do cls = RH10, RH90
            if (RegPar(wsl, cls)%fc /= error) then
                m = m + 1
                xFit(m) = dfloat(cls) * 1d-1
                yFit(m) = RegPar(wsl, cls)%fc
            end if
        end do

        !> Extrapolate cut-off frequency with different policies, depending on how
        !> many have been correctly calculated
        if (m >= 4) then
            !> Perform regression if there are at least 3 RH/fc pairs
            TFShape = 'exponential'
            allocate(fvec(m), fjac(m, npar_EXP))
            call lmder1(fcn, m, npar_EXP, EXPPar, fvec, fjac, tol, info, ipvt_EXP)
            if ((EXPPar(1) == -2d0 .and. EXPPar(2) == -1.d0 .and. EXPPar(3) == -2.d0) &
                .or. info < 1 .or. info > 3) then
                EXPPar(1:3) = error
            end if
            deallocate(fvec, fjac)
            RegPar(wsl, dum)%e1 = EXPPar(1)
            RegPar(wsl, dum)%e2 = EXPPar(2)
            RegPar(wsl, dum)%e3 = EXPPar(3)
        elseif(m == 3 .or. m == 2) then
            cnt = 0
            mean_fc = 0d0
            do cls = RH10, RH90
                if (RegPar(wsl, cls)%fc /= error) then
                    cnt = cnt + 1
                    mean_fc = mean_fc + RegPar(wsl, cls)%fc
                    if (cnt == m) then
                        mean_fc = mean_fc / cnt
                        RegPar(wsl, RH10:RH90)%fc = mean_fc
                        RegPar(wsl, dum)%e1 = 1d-15
                        RegPar(wsl, dum)%e2 = 1d-15
                        RegPar(wsl, dum)%e3 = dlog(mean_fc)
                        exit
                    end if
                end if
            end do
        elseif(m == 1) then
            do cls = RH10, RH90
                if (RegPar(wsl, cls)%fc /= error) then
                    RegPar(wsl, RH10:RH90)%fc = RegPar(wsl, cls)%fc
                    RegPar(wsl, dum)%e1 = 1d-15
                    RegPar(wsl, dum)%e2 = 1d-15
                    RegPar(wsl, dum)%e3 = dlog(RegPar(wsl, cls)%fc)
                    exit
                end if
            end do
        elseif(m == 0) then
            RegPar(wsl, dum)%e1 = error
            RegPar(wsl, dum)%e2 = error
            RegPar(wsl, dum)%e3 = error
        endif

        !> Now that the regression is done, sets Fc which are at -9999, to the closest valid value.
        if (m >= 3) then
            !> For low RH classes, uses higher ones
            do cls = RH50, RH10, - 1
                if (RegPar(wsl, cls)%fc == error) then
                    do cls2 = cls + 1, RH90
                        if (RegPar(wsl, cls2)%fc /= error) then
                            RegPar(wsl, cls)%fc = RegPar(wsl, cls2)%fc
                            RegPar(wsl, cls)%Fn = RegPar(wsl, cls2)%Fn
                            exit
                        end if
                    enddo
                end if
            end do
            !> For high RH classes, uses lower ones
            do cls = RH60, RH90
                if (RegPar(wsl, cls)%fc == error) then
                    do cls2 = cls - 1, RH10, - 1
                        if (RegPar(wsl, cls2)%fc /= error) then
                            RegPar(wsl, cls)%fc = RegPar(wsl, cls2)%fc
                            RegPar(wsl, cls)%Fn = RegPar(wsl, cls2)%Fn
                            exit
                        end if
                    enddo
                end if
            end do
        end if
    end if

    end do

    !> The primary's coefficients are also the project's.
    !>
    !> RegPar(dum, dum) is what the assessment file's standalone exponential
    !> section carries and what a hygrometer with no fit of its own falls back
    !> to, so it keeps holding the primary's - which is exactly what it held
    !> when this routine fitted nothing else.
    RegPar(dum, dum)%e1 = RegPar(primary, dum)%e1
    RegPar(dum, dum)%e2 = RegPar(primary, dum)%e2
    RegPar(dum, dum)%e3 = RegPar(primary, dum)%e3

    if (allocated(xFit)) deallocate(xFit)
    if (allocated(yFit)) deallocate(yFit)
    if (allocated(ddum)) deallocate(ddum)

    call LogSay('Done.')
end subroutine FitRh2Fco
