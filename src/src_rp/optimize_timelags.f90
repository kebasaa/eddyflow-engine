!***************************************************************************
! optimize_timelags.f90
! ---------------------
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
! \brief       Calculate most likely time-lag and range of variation for
!              all gases
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OptimizeTimelags(toSet, nrow, actn, M, h2o_n, MM, cls_size)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: M
    integer, intent(in) :: MM
    real(kind = dbl), intent(in) :: cls_size
    type (TimeLagDatasetType), intent(in) :: toSet(nrow)
    integer, intent(in) :: actn(M)
    !> Per class and per slot: only classed water slots are filled.
    integer, intent(out) :: h2o_n(MM, M)
    !> local variables
    integer :: gas
    integer :: wsl
    integer :: cls
    integer :: i
    integer :: N
    integer :: nn
    integer :: first
    integer :: last
    integer :: nup, ndw
    !> The report prints this number, so it lives in m_typedef where both
    !> can reach it.
    integer, parameter :: min_numerosity = toMinH2OClassN
    real(kind = dbl) :: medx
    real(kind = dbl) :: medup, meddw
    real(kind = dbl), allocatable :: tmpx(:)
    real(kind = dbl), allocatable :: devx(:)
    real(kind = dbl), allocatable :: tmpup(:),tmpdw(:)
    real(kind = dbl), allocatable :: devup(:),devdw(:)
    real(kind = dbl) :: MAD
    real(kind = dbl) :: MADup, MADdw
    real(kind = dbl) :: mvec
    real(kind = dbl) :: sdvec
    real(kind = dbl) :: tmpvec(nrow)
    real(kind = dbl) ,parameter :: min_range = 0.3d0
    logical, external :: WaterSlotClassed
    include '../src_common/interfaces_1.inc'

!TO REFINE integer :: read_status
!TO REFINE integer :: h2on

    !> The site's water is forced present here so its RH-class treatment below
    !> runs; slot six is water only when record two holds it.
    wsl = PrimaryWaterOutSlot()
    E2Col(wsl)%present = .true.

    !> Every class starts at nothing counted, before the gas loop rather than
    !> inside the RH block below - h2o_n is intent(out), so this routine owes
    !> the caller a defined array on EVERY path, and the RH block is reached
    !> only when the water gas has determinations to classify. The guard at
    !> "if (N <= 0) cycle" skips it for a water gas with none, and the counts
    !> then stayed undefined all the way into the report, which printed the
    !> stack. Zero is also the true answer for such a class.
    h2o_n = 0

    !> And the classes themselves, for the same reason and one worse. toH2O is
    !> not a local: it is module state, so on a run whose water never reaches
    !> the RH block it holds whatever the LAST caller left there - or, on the
    !> first call, nothing defined at all. Only the block below ever cleared
    !> it, and that block is exactly the one such a run does not enter.
    !>
    !> Two things then read it. The test at the foot of this routine decides
    !> from toH2O(1)%def whether any class could be filled, so the alert it
    !> guards fired or did not according to stale memory; and SetTimelags
    !> takes the water detection window from these same classes, which is a
    !> flux consequence rather than a cosmetic one.
    !>
    !> Seen on base_pwb_prefilt, where no gas settles anywhere: the class
    !> table printed -9999 or 0.00 for the same input depending on what had
    !> been in memory beforehand.
    toH2O%def = error
    toH2O%min = error
    toH2O%max = error

    do gas = firstGas, lastGas
        if (E2Col(gas)%present) then
            !> All gases, including H2O, are treated here
            toPasGas(gas)%def = error
            toPasGas(gas)%min = error
            toPasGas(gas)%max = error
            N = actn(gas)
            !> A gas the optimiser accepted nothing for keeps the error window
            !> set above, which SetTimelags reads as "no optimised window" and
            !> falls back on the declared one.
            !>
            !> Without this the median of a zero-length array left medx
            !> undefined and the window was built from it. It still passed the
            !> `max > min` guard downstream, because MAD is floored at 0.1, so
            !> a gas with no determinations at all reached the production pass
            !> carrying whatever had been on the stack.
            if (N <= 0) cycle
            allocate (tmpx(N), devx(N))
            tmpx(1:N) = toSet(1:N)%tlag(gas)
            call median(tmpx, N, medx)
            devx(1:N) = dabs(toSet(1:N)%tlag(gas) - medx)
            call median(devx, N, MAD)
            if (MAD < 0.1 ) MAD = 0.1 !< Set a minimum value for MAD
            toPasGas(gas)%def = medx
            toPasGas(gas)%max = medx + (TOSetup%pg_range * MAD / 0.6745d0)
            toPasGas(gas)%min = medx - (TOSetup%pg_range * MAD / 0.6745d0)
            deallocate (tmpx, devx)

            !> Every classed hygrometer, each by its own humidity - see
            !> WaterSlotClassed. This was the designated hygrometer alone.
            if (WaterSlotClassed(gas)) then
                !> Water vapour, the same as above, but for RH classes
                toH2O(:, gas)%def=error
                toH2O(:, gas)%min=error
                toH2O(:, gas)%max=error
                do cls = 1, MM
                    h2o_n(cls, gas) = 0
                    tmpvec = 0d0
                    do i = 1, actn(gas)
                        if(toSet(i)%RH(gas) >= dfloat(cls - 1) * cls_size &
                            .and. toSet(i)%RH(gas) <= dfloat(cls) * cls_size) then
                            h2o_n(cls, gas) = h2o_n(cls, gas) + 1
                            tmpvec(h2o_n(cls, gas)) = toSet(i)%tlag(gas)
                        end if
                    end do
                    N = h2o_n(cls, gas)
                    if (N < min_numerosity) cycle
                    !> Eliminate outliers in each class
                    !> and redefine "short" time-lag set (without outliers)
                    mvec = sum(tmpvec(:)) / N
                    sdvec = dsqrt(sum((tmpvec(:) - mvec)**2))
                    nn = 0
                    do i = 1, N
                        if (tmpvec(i) > mvec + 5.0d0 * sdvec &
                            .or. tmpvec(i) < mvec - 5.0d0 * sdvec) cycle
                        nn = nn + 1
                        tmpvec(nn) = tmpvec(i)
                    end do

                    !> Calculate plausibility range
                    N = nn
                    allocate(tmpx(N), devx(N))
                    allocate(tmpup(N), tmpdw(N))
                    allocate(devup(N), devdw(N))
                    tmpx(1:N) = tmpvec(1:N)
                    call median(tmpx, N, medx)
                    devx(1:N) = tmpvec(1:N) - medx
                    nup = 0
                    ndw = 0
                    do i = 1, N
                        if (devx(i) > 0d0) then
                            nup = nup + 1
                            tmpup(nup) = devx(i)
                        elseif (devx(i) < 0d0) then
                            ndw = ndw + 1
                            tmpdw(ndw) = dabs(devx(i))
                        end if
                    end do
                    call median(tmpup(1:nup), nup, medup)
                    call median(tmpdw(1:ndw), ndw, meddw)
                    devup(1:nup) = dabs(tmpup(1:nup) - medup)
                    devdw(1:ndw) = dabs(tmpdw(1:ndw) - meddw)
                    call median(devup(1:nup), nup, MADup)
                    call median(devdw(1:ndw), ndw, MADdw)
                    toH2O(cls, gas)%def = medx
                    toH2O(cls, gas)%max = medx + (TOSetup%pg_range * MADup / 0.6745d0)
                    toH2O(cls, gas)%min = medx - (TOSetup%pg_range * MADdw / 0.6745d0)
                    if (TOSetup%pg_range * MADup / 0.6745d0 < min_range) &
                        toH2O(cls, gas)%max = medx + min_range
                    if (TOSetup%pg_range * MADdw / 0.6745d0 < min_range) &
                        toH2O(cls, gas)%min = medx - min_range
                    deallocate (tmpx, devx)
                    deallocate (tmpup, tmpdw)
                    deallocate (devup, devdw)
                end do

                !> Adjust time-lags for classes not filled
                !> Detects first good class
                first = 1
                do cls = 1, MM
                    if (h2o_n(cls, gas) > min_numerosity) then
                        first = cls
                        exit
                    end if
                end do
                !> Detects last good class
                last = MM
                do cls = MM, 1, -1
                    if (h2o_n(cls, gas) > min_numerosity) then
                        last = cls
                        exit
                    end if
                end do
                !> Set initial classes equal to first good class
                if (first > 1) toH2O(1:first - 1, gas) = toH2O(first, gas)
                !> For intermediate classes, averages before and after
                if (last > first + 1) then
                    do cls = first + 1, last - 1
                        if (h2o_n(cls, gas) == min_numerosity) then
                            toH2O(cls, gas)%def = (toH2O(cls-1, gas)%def + toH2O(cls+1, gas)%def) * 0.5d0
                            toH2O(cls, gas)%min = (toH2O(cls-1, gas)%min + toH2O(cls+1, gas)%min) * 0.5d0
                            toH2O(cls, gas)%max = (toH2O(cls-1, gas)%max + toH2O(cls+1, gas)%max) * 0.5d0
                        end if
                    end do
                end if
                !> Set traling classes to linear extrapolation of last 2 good ones
                if (last < MM .and. last > 2) then
                    do cls = last + 1, MM
                        toH2O(cls, gas)%def = 2d0 * toH2O(cls - 1, gas)%def - toH2O(cls - 2, gas)%def
                        toH2O(cls, gas)%min = 2d0 * toH2O(cls - 1, gas)%min - toH2O(cls - 2, gas)%min
                        toH2O(cls, gas)%max = 2d0 * toH2O(cls - 1, gas)%max - toH2O(cls - 2, gas)%max
                    end do
                end if
            end if
        end if
    end do

    !> If time-lag optimization failed, switch to covariance maximization.
    !>
    !> Not under PWB. There the aggregate table is a by-product - a summary of
    !> per-gas windows and H2O relative-humidity classes written for the next
    !> run to read - while the method's actual output is the half-hourly
    !> table, one settled lag per period per gas, already complete by the time
    !> this runs. Switching the method here threw all of that away and
    !> processed with covariance maximization instead, reporting every
    !> _TLAG_PWB_SOURCE column as missing.
    !>
    !> What made it fire is water. The classes below need reliable water
    !> detections, and a hygrometer that has few or none leaves them empty -
    !> which says nothing at all about the other gases' lags, and used to be
    !> masked by the summary borrowing a donor from another analyser.
    !>
    !> Asked of the designated hygrometer alone, and only when it is classed
    !> at all. A second hygrometer with no classes falls back to its own
    !> single window, which says nothing about the method. And with
    !> classing off (one class, or an open-path primary) there was never a
    !> table to fill, yet this test - not gated on that - read it as failed
    !> and switched every such run to covariance maximisation.
    !>
    !> Nested, not one .and.: MM can be zero, and Fortran does not promise to
    !> stop before reading toH2O(MM, wsl).
    if (WaterSlotClassed(wsl)) then
        if (toH2O(1, wsl)%def == error .and. toH2O(MM, wsl)%def == error) then
            if (Meth%tlag == 'pwb') then
                call LogSay(' Alert> No H2O relative-humidity classes could be filled for the')
                call LogSay('        aggregate time-lag summary. The half-hourly PWB table is')
                call LogSay('        unaffected and is what this run uses.')
            else
                call ExceptionHandler(43)
                Meth%tlag = 'maxcov'
                TimeLagOptSelected = .false.
            end if
        end if
    end if
end subroutine OptimizeTimelags

!***************************************************************************
!
! \brief       Whether a slot's time lag is classed by relative humidity.
! \author      Jonathan Muller
! \note        Every closed-path hygrometer, each by its own humidity: water
!              adsorbs on the wall of the tube it is drawn through, so its lag
!              drifts with humidity, and each analyser has its own tube. An
!              open-path hygrometer has no tube and is not classed - which is
!              what forcing the class count to one for an open-path primary
!              used to say, for the primary alone and for every other
!              hygrometer with it.
!
!              Only the designated hygrometer used to be classed. A second one
!              took a single window, and its determinations were gated on the
!              designated hygrometer's humidity rather than its own.
!***************************************************************************
logical function WaterSlotClassed(slot)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: slot
    logical, external :: GasSlotIsWater

    WaterSlotClassed = .false.
    if (slot < firstGas .or. slot > lastGas) return
    if (TOSetup%h2o_nclass <= 1) return
    if (.not. GasSlotIsWater(slot)) return
    if (.not. E2Col(slot)%present) return
    if (E2Col(slot)%instr%path_type == 'open') return
    WaterSlotClassed = .true.
end function WaterSlotClassed
