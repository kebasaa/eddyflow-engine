!***************************************************************************
! metek_head_correction.f90
! -------------------------
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
! \brief       Three-dimensional flow distortion of a Metek USA-1 sonic.
! \author      Jonathan Muller
! \note
!              WHAT IT IS. The transducers and their supports deflect the
!              flow before it is measured, by an amount that depends on where
!              the wind comes from. Metek measured that in a wind tunnel and
!              published it as three tables of Fourier coefficients over
!              elevation angle - one each for the wind speed, the azimuth and
!              the elevation - which are evaluated at 3, 6 and 9 times the
!              azimuth. EddyUH applies them in Functions_Library/METEK_HC.m,
!              and this is that routine.
!
!              THE TABLES ARE NOT SHIPPED. They are Metek GmbH's measurements
!              (6th October 2003), redistributed by EddyUH under the
!              University of Helsinki's own agreement. This program does not
!              carry them: head_corr_dir names a directory holding
!              phicorr.dat, ucorr.dat and alphacorr.dat, and without those
!              files the correction declines and says so. Anyone entitled to
!              the tables can point at their own copy.
!
!              ONE INNER-BAR MODELS ONLY, as METEK_HC.m's own header says.
!              Nothing here checks that, because nothing in the metadata
!              distinguishes the variants.
!
!              THE TWO MODES. 'raw' is for a logger that applied no
!              correction. 'undo_2d' is for one that already applied Metek's
!              online two-dimensional correction: that is removed first, with
!              the closed form Metek publishes for it, and the full
!              three-dimensional correction applied to what is left. Applying
!              the 3-D correction on top of the 2-D one without removing it
!              would count the horizontal part twice.
!
!              EDGES OF THE TABLE. Below -50 degrees of elevation the first
!              row is held, above +45 the last. Metek measured that range and
!              extrapolating a ninth-harmonic fit past it would produce
!              confident nonsense.
! \sa          InclinometerTilt
!***************************************************************************
subroutine MetekHeadCorrection(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    !> The tables: 20 rows of elevation from -50 to +45 degrees in steps of
    !> five, each with a mean and three cosine/sine pairs.
    integer, parameter :: nrows = 20
    real(kind = dbl), parameter :: grd_first = -50d0
    real(kind = dbl), parameter :: grd_step = 5d0
    real(kind = dbl), save :: phi_tab(nrows, 7)
    real(kind = dbl), save :: u_tab(nrows, 7)
    real(kind = dbl), save :: alpha_tab(nrows, 7)
    logical, save :: loaded = .false.
    logical, save :: refused = .false.
    integer :: i
    real(kind = dbl) :: uu, vv, ww
    real(kind = dbl) :: alfa, vc, delta, corr_z
    real(kind = dbl) :: u_meas, phi_rad, phi_deg
    real(kind = dbl) :: cf(7), cu(7), ca(7)
    real(kind = dbl) :: phi_c, u_c, alfa_c
    real(kind = dbl) :: u_true, alfa_true, phi_true

    if (RPSetup%head_corr_meth == 'none') return
    if (refused) return

    if (.not. loaded) then
        call LoadTables()
        if (refused) return
        loaded = .true.
    end if

    do i = 1, nrow
        if (Set(i, u) == error .or. Set(i, v) == error &
            .or. Set(i, w) == error) cycle
        uu = Set(i, u)
        vv = Set(i, v)
        ww = Set(i, w)

        !> Azimuth of the incoming flow, in Metek's convention.
        alfa = -datan2(vv, -uu)
        if (alfa < 0d0) alfa = alfa + 2d0 * p

        if (RPSetup%head_corr_meth == 'undo_2d') then
            !> Back out the online two-dimensional correction before applying
            !> the three-dimensional one.
            vc = dsqrt(uu**2 + vv**2)
            delta = 1d0 + 0.015d0 * dsin(3d0 * alfa + p / 6d0)
            corr_z = 0.031d0 * vc * (dsin(3d0 * alfa) - 1d0)
            uu = uu / delta
            vv = vv / delta
            ww = ww - corr_z
        end if

        u_meas = dsqrt(uu**2 + vv**2 + ww**2)
        phi_rad = -datan2(ww, dsqrt(uu**2 + vv**2))
        phi_deg = phi_rad * 180d0 / p

        call Interpolate(phi_tab, phi_deg, cf)
        call Interpolate(u_tab, phi_deg, cu)
        call Interpolate(alpha_tab, phi_deg, ca)

        phi_c = Harmonics(cf, alfa)
        u_c = Harmonics(cu, alfa)
        alfa_c = Harmonics(ca, alfa)

        u_true = u_meas * u_c
        alfa_true = alfa + alfa_c * p / 180d0
        phi_true = phi_rad + phi_c * p / 180d0

        Set(i, u) = -u_true * dcos(alfa_true) * dcos(phi_true)
        Set(i, v) = -u_true * dsin(alfa_true) * dcos(phi_true)
        Set(i, w) = -u_true * dsin(phi_true)
    end do

contains

    !> The three tables, from wherever the project says they live. Read once
    !> per run: they do not change between periods, and re-reading them for
    !> every one of a month's worth would be the dominant cost of this
    !> routine.
    subroutine LoadTables()
        logical :: ok

        call ReadOneTable('phicorr.dat', phi_tab, ok)
        if (ok) call ReadOneTable('ucorr.dat', u_tab, ok)
        if (ok) call ReadOneTable('alphacorr.dat', alpha_tab, ok)
        if (.not. ok) then
            refused = .true.
            call LogSay('  Metek head correction asked for, but the look-up &
                &tables were not found or not readable under "' &
                // trim(RPSetup%head_corr_dir) // '" - skipped for this run. &
                &phicorr.dat, ucorr.dat and alphacorr.dat are Metek GmbH data &
                &and are not shipped with this program.')
        end if
    end subroutine LoadTables

    subroutine ReadOneTable(name, tab, ok)
        character(*), intent(in) :: name
        real(kind = dbl), intent(out) :: tab(nrows, 7)
        logical, intent(out) :: ok
        integer :: unt, k, ios
        real(kind = dbl) :: grd
        character(PathLen) :: path

        ok = .false.
        if (len_trim(RPSetup%head_corr_dir) == 0) return
        path = trim(adjustl(RPSetup%head_corr_dir)) // slash // name
        open(newunit = unt, file = trim(path), status = 'old', iostat = ios)
        if (ios /= 0) return
        do k = 1, nrows
            !> Comma separated, and the first field is the elevation the row
            !> stands for. It is read and discarded: the grid is a fixed five
            !> degrees from -50, and a file that disagrees would silently
            !> shift every correction rather than fail, so the row count is
            !> what this checks.
            read(unt, *, iostat = ios) grd, tab(k, 1:7)
            if (ios /= 0) then
                close(unt)
                return
            end if
        end do
        close(unt)
        ok = .true.
    end subroutine ReadOneTable

    !> Linear interpolation between two rows five degrees apart, holding the
    !> end rows outside the measured range.
    subroutine Interpolate(tab, deg, c)
        real(kind = dbl), intent(in) :: tab(nrows, 7)
        real(kind = dbl), intent(in) :: deg
        real(kind = dbl), intent(out) :: c(7)
        integer :: idx
        real(kind = dbl) :: frac

        if (deg < grd_first) then
            c = tab(1, :)
            return
        end if
        if (deg > grd_first + grd_step * dble(nrows - 1)) then
            c = tab(nrows, :)
            return
        end if
        idx = int((deg - grd_first) / grd_step) + 1
        if (idx < 1) idx = 1
        if (idx >= nrows) then
            c = tab(nrows, :)
            return
        end if
        frac = (deg - (grd_first + grd_step * dble(idx - 1))) / grd_step
        c = tab(idx, :) + frac * (tab(idx + 1, :) - tab(idx, :))
    end subroutine Interpolate

    !> Mean plus the third, sixth and ninth harmonics of the azimuth. The
    !> column order is the one the files carry: C0, C3, S3, C6, S6, C9, S9.
    real(kind = dbl) function Harmonics(c, az)
        real(kind = dbl), intent(in) :: c(7)
        real(kind = dbl), intent(in) :: az

        Harmonics = c(1) &
            + c(2) * dcos(3d0 * az) + c(3) * dsin(3d0 * az) &
            + c(4) * dcos(6d0 * az) + c(5) * dsin(6d0 * az) &
            + c(6) * dcos(9d0 * az) + c(7) * dsin(9d0 * az)
    end function Harmonics

end subroutine MetekHeadCorrection
