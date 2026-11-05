!***************************************************************************
! potential_radiation.f90
! -----------------------
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
! \brief       Computes potential radiation, needed for identification \n
!              of night and daytime
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function PotentialRadiation(latit, longit, utc_off) result(RP)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: latit
    !> Longitude [deg E] and the hours the timestamps are ahead of UTC
    real(kind = dbl), intent(in) :: longit
    real(kind = dbl), intent(in) :: utc_off
    !> local variables
    integer :: i
    real(kind = dbl) :: t_h
    real(kind = dbl) :: clock
    real(kind = dbl) :: solar
    integer :: lDOY
    real(kind = dbl) :: lat_rad
    real(kind = dbl) :: rpot
    real(kind = dbl) :: theta       !< Some sort of solar angular position
    real(kind = dbl) :: theta_rad       !< Some sort of solar angular position
    real(kind = dbl) :: delta       !< solar declination (in radiants)
    real(kind = dbl) :: ET          !< Equation of Time
    real(kind = dbl) :: omega       !< solar hour angle
    real(kind = dbl) :: LAS         !< boh
    real(kind = dbl), parameter :: solar_constant = 1376d0
    real(kind = dbl) :: RP(17568)

    !> Compute potential radiation, based on lat/long, on a 30 min basis.
    !>
    !> RP(i) is the radiation at the middle of the year's i-th half-hour -
    !> the one that ENDS at i half-hours, which is how DateTimeToHalfHourNumber
    !> numbers a period by its end - so a period reads its own RP. It used to
    !> sit at i + 1/2 half-hours, and the readers averaged i - 1 and i - 2:
    !> the period's start, a quarter of an hour early.
    !>
    !> And at apparent solar time, not clock time. Clock noon was taken for
    !> solar noon whatever the longitude and time zone, and the equation of
    !> time below was computed and never used: at a site on UTC far from
    !> Greenwich the whole curve sat hours off, and where biomet has neither
    !> Rg nor PPFD that curve decides day and night - and with it Burba's day
    !> or night regression for an LI-7500.
    do i = 1, 17568
        !> Hours since the start of the year, standard time of the timestamps
        t_h = (dble(i) - 0.5d0) * 0.5d0
        !> Compute day
        lDOY = int(t_h / 24d0) + 1
        clock = t_h - 24d0 * dble(lDOY - 1)
        !> Compute theta
        theta = 2d0 * p * (lDOY - 1d0) / 365d0;
        !> Compute equation of time
        ET = 0.000075d0 + 0.001868d0 * dcos(theta) - 0.032077d0 * dsin(theta) &
          - 0.014615d0 * dcos(2d0 * theta) - 0.040849d0 * dsin(2d0 * theta)
        !> Compute solar declination in radiants
        delta = 0.006918d0 - 0.399912d0 * dcos(theta) + 0.070257d0 * dsin(theta) &
             - 0.006758d0 * dcos(2d0 * theta) + 0.000907d0 * dsin(2d0 * theta) &
             -0.002697d0 * dcos(3d0 * theta) + 0.00148d0 * dsin(3d0 * theta)
        !> Apparent solar time: the longitude's offset from the time zone's
        !> meridian, and the equation of time (Spencer 1971, in radians, so
        !> 229.18 minutes per unit)
        solar = clock + (longit / 15d0 - utc_off) + 229.18d0 * ET / 60d0
        !> Compute LAS, hours from solar noon
        LAS = 12.d0 - solar
        LAS = dabs(LAS)
        !> Compute omega
        omega = -15d0 * LAS
        !> Latitude and longitude in radiants
        lat_rad = latit * p / 180d0;
        !> Compute theta_rad
        theta_rad = dacos(dsin(delta) * dsin(lat_rad) &
                 + dcos(delta) * dcos(lat_rad) * dcos(omega * p / 180d0))
        !> Compute extraterrestrial solar radiation
        rpot = solar_constant * (1.00011d0 + 0.034221d0 * dcos(theta) + 0.00128d0 * dsin(theta) &
             + 0.000719d0 * dcos(2d0 * theta) + 0.000077d0 * dsin(2d0 * theta))
        RP(i) = rpot*dcos(theta_rad)
        if (RP(i) < 0d0) RP(i) = 0d0
    end do
end function PotentialRadiation
