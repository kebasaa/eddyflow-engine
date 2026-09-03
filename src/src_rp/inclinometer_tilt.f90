!***************************************************************************
! inclinometer_tilt.f90
! ---------------------
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
! \brief       Rotate the wind by inclination angles measured alongside it.
! \author      Jonathan Muller
! \note
!              WHAT IT IS FOR. A mast that leans, or sways, tilts the sonic
!              with it. A planar fit or a double rotation removes the MEAN
!              tilt over an averaging period; neither can remove a tilt that
!              changes within one. An inclinometer logged at the same rate as
!              the wind can, sample by sample, and that is what this does -
!              EddyUH's EC_Software_Common/EddyUH_tiltangle.m.
!
!              WHERE THE ANGLES COME FROM. Three optional raw columns named
!              theta, phi and psi, matched by their metadata variable name
!              and read through the custom-column machinery that already
!              carries signal-strength channels. No new record family: there
!              is one sonic, so a name is enough to say which instrument an
!              angle belongs to, unlike a signal strength which has to name
!              its analyser.
!
!              The channels hold the inclinometer's OUTPUT VOLTAGE, not an
!              angle. angle = -asin(V / sensitivity), with the sensitivity in
!              volts per g. A channel that is absent contributes a zero angle,
!              which is what EddyUH does and what leaves that axis alone.
!
!              PSI IS ALWAYS ZERO, even when a psi channel exists. EddyUH
!              reads it, then overwrites it with zeros and comments "not
!              measured" (EddyUH_tiltangle.m:60). Reproduced, because the
!              rotation matrix and the swinging term are both built on that
!              assumption and a real psi would need the whole thing rederived.
!              So this is a two-angle correction with three channels declared.
!
!              THE TWO MODES. 'position' rotates the measured vector by the
!              inclination. 'position_swing' adds a term for the motion of the
!              sonic head as the mast swings, built from the lever arm L from
!              the pivot to the head and an omega assembled from the time
!              derivatives of the angles. The second needs the arm to be
!              right; a wrong arm adds a velocity that is not there.
!
!              THE SWINGING TERM IS A SCALAR, AND THAT IS ODD. EddyUH writes
!              `omega*T*L` with omega 1x3, T 3x3 and L 3x1, so the product is
!              a single number, and it is added to u, v and w ALIKE
!              (EddyUH_tiltangle.m:104). The velocity of a point on a rotating
!              body is omega CROSS (T L), a vector with three different
!              components; a dot product is not that. The units survive - rad
!              per second times metres is metres per second - which is why it
!              is easy to miss.
!
!              Reproduced as EddyUH writes it, because the option exists to
!              reproduce EddyUH's numbers and a silently corrected version
!              would reproduce nothing. Anyone wanting the physical form
!              should treat this mode as unavailable rather than adjust the
!              arm to compensate. Recorded in docs/eddyuh-comparison.md.
!
!              BEFORE ANY ROTATION. This is a hardware correction on the raw
!              wind, like the angle-of-attack calibration beside it, and the
!              double rotation or planar fit that follows expects a series
!              already in the sonic's true frame.
!
!              THE FLUX LOOP ONLY. The two pre-passes - time-lag optimisation
!              and the planar fit assessment - run before DefineUserSet, so
!              the angle columns do not exist yet and this cannot be applied
!              there. MetekHeadCorrection, which needs nothing but the wind,
!              is applied in all three. So a planar fit is fitted to a wind
!              that has not had its within-period tilt removed, while the
!              fluxes are computed from one that has.
!
!              That is a real inconsistency and not a rounding one, but it is
!              second order: the plane is fitted to per-period MEAN winds, and
!              this correction is aimed at variation WITHIN a period. Moving
!              DefineUserSet earlier would close it, at the cost of a layout
!              change to the one part of the program where a layout fault is
!              silent. Left as it is, and recorded here, until a dataset
!              exists that this correction actually moves.
! \sa          MetekHeadCorrection, TiltCorrection
!***************************************************************************
subroutine InclinometerTilt(Set, nrow, ncol, UserSet, unrow, uncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Passed rather than reached for: UserSet is a local of the main
    !> program, and only UserCol - which says what each of its columns is -
    !> is a module global.
    integer, intent(in) :: unrow, uncol
    real(kind = dbl), intent(in) :: UserSet(unrow, uncol)
    !> local variables
    integer :: i
    integer :: n_found
    logical :: metek
    real(kind = dbl), allocatable :: theta(:), phi(:), psi(:)
    real(kind = dbl), allocatable :: thetadot(:), phidot(:), psidot(:)
    real(kind = dbl) :: t(3, 3)
    real(kind = dbl) :: v_obs(3), v_true(3), omega(3)
    real(kind = dbl) :: dt
    logical :: swing

    if (RPSetup%tilt_sensor_meth == 'none') return
    swing = RPSetup%tilt_sensor_meth == 'position_swing'

    allocate(theta(nrow), phi(nrow), psi(nrow))
    n_found = 0
    call AngleChannel('theta', theta, n_found)
    call AngleChannel('phi', phi, n_found)
    !> Read for the count and the message only. See the note above: the angle
    !> itself is forced to zero exactly as EddyUH does.
    call AngleChannel('psi', psi, n_found)
    psi = 0d0

    if (n_found == 0) then
        call LogSay('  Inclinometer tilt correction asked for, but no theta, &
            &phi or psi column was found - skipped.')
        deallocate(theta, phi, psi)
        return
    end if

    !> Metek reports a left-handed frame, so v is flipped into a right-handed
    !> one for the rotation and flipped back afterwards.
    metek = index(E2Col(u)%instr%model, 'usa1') /= 0

    allocate(thetadot(nrow), phidot(nrow), psidot(nrow))
    if (swing) then
        dt = 1d0 / Metadata%ac_freq
        call CentredGradient(theta, thetadot, nrow, dt)
        call CentredGradient(phi, phidot, nrow, dt)
        call CentredGradient(psi, psidot, nrow, dt)
    end if

    do i = 1, nrow
        if (Set(i, u) == error .or. Set(i, v) == error &
            .or. Set(i, w) == error) cycle

        v_obs(1) = Set(i, u)
        v_obs(2) = Set(i, v)
        v_obs(3) = Set(i, w)
        if (metek) v_obs(2) = -v_obs(2)

        call RotationFromAngles(theta(i), phi(i), psi(i), t)
        v_true = matmul(t, v_obs)

        if (swing) then
            omega(1) = -thetadot(i) * dsin(psi(i)) &
                + phidot(i) * dcos(theta(i)) * dcos(psi(i))
            omega(2) = thetadot(i) * dcos(psi(i)) &
                + phidot(i) * dcos(theta(i)) * dsin(psi(i))
            omega(3) = psidot(i) - phidot(i) * dsin(theta(i))
            !> See the note above: one number, added to all three components.
            v_true = v_true &
                + dot_product(matmul(omega, t), RPSetup%tilt_arm)
        end if

        if (metek) v_true(2) = -v_true(2)
        Set(i, u) = v_true(1)
        Set(i, v) = v_true(2)
        Set(i, w) = v_true(3)
    end do

    deallocate(theta, phi, psi, thetadot, phidot, psidot)

contains

    !> One angle channel, from the custom columns, as an angle in radians.
    !> Absent leaves it at zero, which rotates that axis not at all.
    subroutine AngleChannel(name, angle, found)
        character(*), intent(in) :: name
        real(kind = dbl), intent(out) :: angle(nrow)
        integer, intent(inout) :: found
        integer :: j, k
        real(kind = dbl) :: s
        character(32) :: label

        angle = 0d0
        do j = 1, min(NumUserVar, uncol)
            label = UserCol(j)%var
            call lowercase(label)
            if (trim(adjustl(label)) /= name) cycle

            do k = 1, min(nrow, unrow)
                if (UserSet(k, j) == error) cycle
                !> asin is only defined on [-1, 1]. A reading past full scale
                !> is a fault in the channel, not a steeper mast, so it is
                !> clamped rather than allowed to produce a NaN that would
                !> then propagate into every wind component.
                s = UserSet(k, j) / RPSetup%tilt_sensor_v_g
                if (s > 1d0) s = 1d0
                if (s < -1d0) s = -1d0
                angle(k) = -dasin(s)
            end do
            if (RPSetup%tilt_lpf_s > 0d0) &
                call MovingAverage(angle, nrow, &
                    max(1, nint(RPSetup%tilt_lpf_s * Metadata%ac_freq)))
            found = found + 1
            return
        end do
    end subroutine AngleChannel

end subroutine InclinometerTilt

!***************************************************************************
!
! \brief       The 3-2-1 rotation matrix EddyUH_tiltangle.m:83-85 builds.
! \author      Jonathan Muller
! \note        Separated so the static check can exercise it on its own, and
!              because the swinging term needs the same matrix.
!***************************************************************************
subroutine RotationFromAngles(theta, phi, psi, t)
    use m_rp_global_var
    implicit none
    real(kind = dbl), intent(in) :: theta, phi, psi
    real(kind = dbl), intent(out) :: t(3, 3)

    t(1, 1) = dcos(theta) * dcos(psi)
    t(1, 2) = dsin(theta) * dsin(phi) * dcos(psi) - dcos(phi) * dsin(psi)
    t(1, 3) = dcos(phi) * dsin(theta) * dcos(psi) + dsin(phi) * dsin(psi)
    t(2, 1) = dcos(theta) * dsin(psi)
    t(2, 2) = dsin(phi) * dsin(theta) * dsin(psi) + dcos(phi) * dcos(psi)
    t(2, 3) = dcos(phi) * dsin(theta) * dsin(psi) - dsin(phi) * dcos(psi)
    t(3, 1) = -dsin(theta)
    t(3, 2) = dsin(phi) * dcos(theta)
    t(3, 3) = dcos(phi) * dcos(theta)
end subroutine RotationFromAngles

!***************************************************************************
!
! \brief       Centred difference, one-sided at the two ends - MATLAB's
!              gradient(), which EddyUH uses.
! \author      Jonathan Muller
!***************************************************************************
subroutine CentredGradient(x, dx, n, h)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(out) :: dx(n)
    real(kind = dbl), intent(in) :: h
    integer :: i

    dx = 0d0
    if (n < 2 .or. h <= 0d0) return
    dx(1) = (x(2) - x(1)) / h
    dx(n) = (x(n) - x(n - 1)) / h
    do i = 2, n - 1
        dx(i) = (x(i + 1) - x(i - 1)) / (2d0 * h)
    end do
end subroutine CentredGradient

!***************************************************************************
!
! \brief       Centred running mean over a fixed window, in place.
! \author      Jonathan Muller
! \note        EddyUH's movav(). Shorter windows at the two ends rather than
!              padding, so the filtered series keeps its length and its ends
!              are not pulled toward zero.
!***************************************************************************
subroutine MovingAverage(x, n, width)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer, intent(in) :: width
    integer :: i, lo, hi, half
    real(kind = dbl) :: smoothed(n)

    if (width <= 1) return
    half = width / 2
    do i = 1, n
        lo = max(1, i - half)
        hi = min(n, i + half)
        smoothed(i) = sum(x(lo:hi)) / dble(hi - lo + 1)
    end do
    x = smoothed
end subroutine MovingAverage
