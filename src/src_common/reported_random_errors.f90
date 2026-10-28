!***************************************************************************
! reported_random_errors.f90
! --------------------------
! Copyright (C) 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyFlow (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyFlow (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief       The random errors as they are reported: spectrally corrected
!              like the flux they belong to, when ru_spec_corr asks for it.
! \author      Jonathan Muller
!
! \note
! Every random uncertainty method estimates the error of the covariance as
! measured, before any correction. The flux beside it has been through the
! band-pass spectral correction, which multiplies it by BPCF%of. Langford et
! al. (2015): "any corrections applied to the calculated fluxes
! (e.g. low or high frequency loss corrections) should also be applied to the
! random error". So, with ru_spec_corr on, each random error is
! multiplied by the spectral correction factor of its own flux - Tau by
! w_u's, H by w_ts's, LE and ET by the designated hygrometer's, each gas by
! its own - and the Langford limit of detection, alpha times the random
! error, follows.
!
! Off by default, and for a project without the key: EddyPro reports the
! random error of the uncorrected covariance, and the defaults reproduce
! EddyPro's outputs.
!
! Only the spectral factor: it is the one correction that scales the flux as
! a whole. The density (WPL) and sonic-temperature terms are added to it,
! from other covariances whose errors no method here estimates.
!
! Applied where the factors exist and the final values are written: by RP
! when it writes the final outputs itself (FCC does not follow), by FCC when
! it does. RP's FLUXNET file is FCC's input, and stays uncorrected whenever
! FCC follows. FCC can also be run on the FLUXNET file of a project RP
! finished on its own - the regression harness always does - and that file is
! corrected already; FCC tells the two apart with FccFollowsProject, RP's own
! test, evaluated before FCC changes any of the settings it reads. So nothing
! is corrected twice.
!
! Stored values are never changed: Essentials%rand_uncer and lEx%rand_uncer
! stay as measured, and the callers write the copies this returns.
!***************************************************************************
subroutine ReportedRandomErrors(apply, re_in, le_in, et_in, re_out, le_out, et_out)
    use m_common_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: apply
    real(kind = dbl), intent(in) :: re_in(E2NumVar)
    real(kind = dbl), intent(in) :: le_in
    real(kind = dbl), intent(in) :: et_in
    real(kind = dbl), intent(out) :: re_out(E2NumVar)
    real(kind = dbl), intent(out) :: le_out
    real(kind = dbl), intent(out) :: et_out
    !> local variables
    integer :: var
    integer :: wsl
    include 'interfaces_1.inc'

    re_out = re_in
    le_out = le_in
    et_out = et_in
    if (.not. (apply .and. RUsetup%spec_corr)) return

    re_out(u) = Scaled(re_in(u), BPCF%of(w_u))
    re_out(ts) = Scaled(re_in(ts), BPCF%of(w_ts))
    do var = firstGas, min(lastGas, GHGNumVar)
        re_out(var) = Scaled(re_in(var), BPCF%of(var))
    end do
    wsl = PrimaryWaterSlot()
    if (wsl >= firstGas .and. wsl <= min(lastGas, GHGNumVar)) then
        le_out = Scaled(le_in, BPCF%of(wsl))
        et_out = Scaled(et_in, BPCF%of(wsl))
    end if

contains

    !> A missing error stays missing; a missing or unusable factor leaves the
    !> error as measured rather than inventing one.
    pure function Scaled(re, factor) result(res)
        real(kind = dbl), intent(in) :: re
        real(kind = dbl), intent(in) :: factor
        real(kind = dbl) :: res

        if (re == error .or. factor == error .or. factor <= 0d0) then
            res = re
        else
            res = re * factor
        end if
    end function Scaled
end subroutine ReportedRandomErrors

!***************************************************************************
!
! rief       Whether this project's final outputs are written by FCC: RP's
!              test, shared so FCC can tell whether the ex record it reads
!              carries corrected random errors. Must be evaluated straight
!              after the project is read - FCC later falls back to Moncrieff
!              and drops the ensemble outputs when it finds no spectra.
! uthor      Jonathan Muller
!***************************************************************************
logical function FccFollowsProject()
    use m_common_global_var
    implicit none

    FccFollowsProject = EddyFlowProj%out_avrg_cosp &
        .or. EddyFlowProj%out_avrg_spec &
        .or. (EddyFlowProj%hf_meth /= 'none' &
        .and. EddyFlowProj%hf_meth /= 'moncrieff_97' &
        .and. EddyFlowProj%hf_meth /= 'massman_00' &
        .and. EddyFlowProj%hf_meth /= 'eugster_senn_95')
end function FccFollowsProject
