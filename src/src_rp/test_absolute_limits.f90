!***************************************************************************
! test_absolute_limits.f90
! ------------------------
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
! \brief       Checks for data outside realistic ranges, hard-flag \n
!              file accordingly and eliminate those values if requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestAbsoluteLimits(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    logical :: printout
    !> local variables
    integer :: i = 0
    integer :: cnt1 = 0
    integer :: cnt2 = 0
    integer :: hflags(GHGNumVar)
    real(kind = dbl) :: HorVel
    !> Molar-density scale and rough outlier ceiling of the gas under test.
    !> Water differs from every other gas in both, because its mole fraction
    !> is reported in mmol mol-1 rather than umol mol-1.
    real(kind = dbl) :: dens_scale
    real(kind = dbl) :: rough_max


    if (printout) write(*, '(a)', advance = 'no') '   Absolute limits test..'

    !> initializations
    hflags = 0

    !> Flag and filter wind components
    cnt1 = 0
    cnt2 = 0
    do i = 1, N
        if (all(Set(i, u:w) /= error)) then
            !> Horizontal wind
            HorVel = sqrt((Set(i, u)**2) + (Set(i, v)**2))
            if (HorVel > al%u_max) then
                cnt1 = cnt1 + 1
                hflags(u) = 1
                hflags(v) = 1
                if(RPsetup%filter_al) then
                    Set(i, u) = error
                    Set(i, v) = error
                end if
            end if
            !> Vertical wind
            if (abs(Set(i, w)) > al%w_max) then
                cnt2 = cnt2 + 1
                hflags(w) = 1
                if(RPsetup%filter_al) Set(i, w) = error
            end if
        end if
    end do
    Essentials%al_s(u) = cnt1
    Essentials%al_s(v) = cnt1
    Essentials%al_s(w) = cnt2

    !> Flag and filter sonic temperature
    Essentials%al_s(ts) = count(Set(:, ts) /= error .and. &
                                (Set(:, ts) < al%t_min + 273.15d0 .or. &
                                 Set(:, ts) > al%t_max + 273.15d0))
    if (Essentials%al_s(ts) > 0) hflags(ts) = 1
    if(RPsetup%filter_al) then
        where (Set(:, ts) < al%t_min + 273.15d0 .or. &
            Set(:, ts) > al%t_max + 273.15d0)
            Set(:, ts) = error
        end where
    end if

    !> Flag every gas slot, expressed as [mmol m-3] if molar_density and as a
    !> mole fraction otherwise.
    !>
    !> This was four copies of the same twenty lines, one per historical gas,
    !> which is why gases past the fourth were never tested at all: their flag
    !> stayed at the "not performed" filler and their out-of-range data was
    !> never counted or filtered. The copies differed in exactly two values,
    !> both of which single out water rather than a slot number:
    !>
    !>   - the molar-density scale, because water's mole fraction is reported
    !>     in mmol mol-1 where the other gases use umol mol-1;
    !>   - the ceiling of the rough outlier filter that runs before mean T and
    !>     P are known, in those same units.
    !>
    !> Which slot holds water is settled the same way the analyser block
    !> settles where the krypton coefficients go: by comparing against h2o.
    do i = firstGas, lastGas
        if (.not. E2Col(i)%present) then
            Essentials%al_s(i) = ierror
            hflags(i) = 9
            cycle
        end if
        !> A gas whose limits were never configured cannot be tested.
        !>
        !> Only the four historical slots get limits from the fixed project
        !> keys; past those they come from the per-gas records, and a project
        !> that names a gas without them leaves the pair at 0/0. Testing
        !> against that rejects *every* value, and because filtering runs on
        !> the same pass it replaces the whole series with the error code -
        !> so the gas reaches the flux code with no data at all and is then
        !> dropped by EliminateCorruptedVariables. Absent limits mean the test
        !> was not performed, exactly as for an absent gas; they do not mean
        !> the data is out of range.
        if (al%gas_max(i) <= al%gas_min(i)) then
            Essentials%al_s(i) = ierror
            hflags(i) = 9
            cycle
        end if
        if (i == h2o) then
            dens_scale = StdVair
            rough_max = 80d0
        else
            dens_scale = StdVair * 1d3
            rough_max = 2000d0
        end if
        if (E2Col(i)%measure_type == 'molar_density') then
            Essentials%al_s(i) = count(Set(:, i) /= error .and. &
                                       (Set(:, i) * dens_scale < al%gas_min(i) .or. &
                                        Set(:, i) * dens_scale > al%gas_max(i)))
            !> Actual filtering of molar density gas data is deferred to a later time
            !> when mean T and P are computed. However, some rough filtering is necesasry
            !> to eliminate strong outliers which mess up things if present.
            where (Set(:, i) /= error .and. &
                (Set(:, i) * dens_scale < 0 .or. &
                Set(:, i) * dens_scale > rough_max))
                Set(:, i) = error
            end where
        else
            Essentials%al_s(i) = count(Set(:, i) /= error .and. &
                                       (Set(:, i) < al%gas_min(i) .or. &
                                        Set(:, i) > al%gas_max(i)))
            !> If filtering is requested and data is mole fraction / mixing ratio
            !> eliminate OOR data
            if(RPsetup%filter_al) then
                where (Set(:, i) /= error .and. &
                       (Set(:, i) < al%gas_min(i) .or. &
                        Set(:, i) > al%gas_max(i)))
                    Set(:, i) = error
                end where
            end if
        end if
        if (Essentials%al_s(i) > 0) hflags(i) = 1
    end do

    !> Pack one digit per variable into the flag string
    call PackFlagString(hflags, GHGNumVar, CharHF%al)


    if (RPsetup%filter_al) then
        where (E2Col(u:lastGas)%present)
            Essentials%al_s(u:lastGas) = 0
        end where
    end if

    if (printout) write(*,'(a)') ' Done.'
end subroutine TestAbsoluteLimits
