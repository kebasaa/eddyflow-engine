!***************************************************************************
! pfd_handle.f90
! ---------------
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
! \brief       Post-flux despiking (test_pfd): accumulate NEE/H/LE across
!              the whole FCC run, then despike each series once the run is
!              complete, writing a dedicated output CSV.
!
! \details     Mirrors the RP side's PwbTimelagCache/PostProcessPwbTimelagCache
!              shape (pwb_timelag_handle.f90): a growable per-period cache,
!              populated during the main loop, consumed once after it by a
!              routine that needs the whole run in hand. Off by default;
!              StorePfdCache is only called when EddyFlowProj%test_pfd is
!              set, so a run with the feature off never allocates the cache
!              or writes the new file.
!
!              The actual despiking arithmetic is m_flux_despike_core.f90 -
!              this file is only the FCC-side plumbing: accumulate, sort,
!              derive the averaging period's periods-per-day, call it once
!              per flux, write the result.
!
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          m_flux_despike_core.f90, pwb_timelag_handle.f90,
!              period_ordering.f90 (SortPeriodsByMinutes/ModalPeriodStep,
!              shared with stor_clean_handle.f90's RP-side cache)
!***************************************************************************

!> Append one period's NEE/H/LE to the whole-run cache. Called from
!> ex_loop, guarded by EddyFlowProj%test_pfd at the call site.
subroutine StorePfdCache(date, time, nee, h, le)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: date, time
    real(kind = dbl), intent(in) :: nee, h, le
    type(PfdCacheEntryType), allocatable :: tmp(:)

    allocate(tmp(PfdCacheN + 1))
    if (PfdCacheN > 0) tmp(1:PfdCacheN) = PfdCache(1:PfdCacheN)
    tmp(PfdCacheN + 1)%date = date
    tmp(PfdCacheN + 1)%time = time
    tmp(PfdCacheN + 1)%nee = nee
    tmp(PfdCacheN + 1)%h = h
    tmp(PfdCacheN + 1)%le = le
    call move_alloc(tmp, PfdCache)
    PfdCacheN = PfdCacheN + 1
end subroutine StorePfdCache


!***************************************************************************
!> Runs once, after the whole FCC run has populated PfdCache. Sorts
!> chronologically, derives periods-per-day from the modal spacing between
!> consecutive periods (robust to a gap or two, unlike trusting the first
!> pair), despikes NEE/H/LE, and writes one CSV with the result.
!***************************************************************************
subroutine PostProcessFluxDespiking()
    use m_fx_global_var
    use m_flux_despike_core, only: DespikeFluxSeries
    implicit none
    integer, allocatable :: ord(:)
    integer(8), allocatable :: tmin(:)
    real(kind = dbl), allocatable :: nee(:), h(:), le(:)
    integer, allocatable :: nee_spike(:), h_spike(:), le_spike(:)
    real(kind = dbl), allocatable :: nee_clean(:), h_clean(:), le_clean(:)
    integer :: i, mfreq, pfd_unit
    integer(8) :: modal_step
    character(PathLen) :: PfdPath, Test_Path
    integer :: dot, open_status
    integer(8), external :: PeriodMinutes, ModalPeriodStep

    if (PfdCacheN < 20) then
        !> Too short for a seasonal-trend decomposition to mean anything -
        !> stlplus.default itself refuses n.p < 4, and despiking()'s own
        !> minimum is ten days of data. Say why nothing was written rather
        !> than fail or silently skip.
        call LogSay('   Post-flux despiking: run has fewer than 20 periods, skipped.')
        return
    end if

    allocate(ord(PfdCacheN), tmin(PfdCacheN))
    do i = 1, PfdCacheN
        ord(i) = i
        tmin(i) = PeriodMinutes(PfdCache(i)%date, PfdCache(i)%time)
    end do
    call SortPeriodsByMinutes(ord, tmin, PfdCacheN)

    !> Periods per day, from the most common spacing between consecutive
    !> sorted periods rather than the first pair, so one irregular gap at
    !> the start of the run cannot mis-set it.
    modal_step = ModalPeriodStep(tmin, ord, PfdCacheN)
    if (modal_step <= 0_8) then
        call LogSay('   Post-flux despiking: could not determine the averaging period, skipped.')
        return
    end if
    mfreq = nint(1440d0 / dble(modal_step))
    if (mfreq < 4) then
        call LogSay('   Post-flux despiking: fewer than 4 periods/day, skipped.')
        return
    end if

    allocate(nee(PfdCacheN), h(PfdCacheN), le(PfdCacheN))
    allocate(nee_spike(PfdCacheN), h_spike(PfdCacheN), le_spike(PfdCacheN))
    allocate(nee_clean(PfdCacheN), h_clean(PfdCacheN), le_clean(PfdCacheN))
    do i = 1, PfdCacheN
        nee(i) = ErrToMissing(PfdCache(ord(i))%nee)
        h(i) = ErrToMissing(PfdCache(ord(i))%h)
        le(i) = ErrToMissing(PfdCache(ord(i))%le)
    end do

    call DespikeFluxSeries(nee, PfdCacheN, mfreq, error, 0.01d0, nee_spike, nee_clean)
    call DespikeFluxSeries(h, PfdCacheN, mfreq, error, 0.01d0, h_spike, h_clean)
    call DespikeFluxSeries(le, PfdCacheN, mfreq, error, 0.01d0, le_spike, le_clean)

    !> Same construction as FLUXNET_Path in init_out_files.f90, with its own
    !> padding, so the two live side by side and neither depends on the
    !> other's open file handle.
    Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
              // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
              // PFD_FilePadding // Timestamp_FilePadding // CsvExt
    dot = index(Test_Path, CsvExt, .true.) - 1
    PfdPath = Test_Path(1:dot) // CsvExt

    open(newunit = pfd_unit, file = trim(PfdPath), status = 'replace', iostat = open_status, encoding = 'utf-8')
    if (open_status /= 0) then
        call LogSay('   Post-flux despiking: could not open output file, skipped.')
        return
    end if

    write(pfd_unit, '(a)') 'TIMESTAMP,NEE,NEE_SPIKE,NEE_CLEANED,H,H_SPIKE,H_CLEANED,LE,LE_SPIKE,LE_CLEANED'
    do i = 1, PfdCacheN
        write(pfd_unit, '(a,a1,a,11(a1,a))') &
            trim(PfdCache(ord(i))%date) // ' ' // trim(PfdCache(ord(i))%time), ',', &
            trim(FormatOrErr(nee(i))), ',', trim(IntToStr(nee_spike(i))), &
            ',', trim(FormatOrErr(nee_clean(i))), &
            ',', trim(FormatOrErr(h(i))), ',', trim(IntToStr(h_spike(i))), &
            ',', trim(FormatOrErr(h_clean(i))), &
            ',', trim(FormatOrErr(le(i))), ',', trim(IntToStr(le_spike(i))), &
            ',', trim(FormatOrErr(le_clean(i)))
    end do
    close(pfd_unit)

    deallocate(ord, tmin, nee, h, le, nee_spike, h_spike, le_spike, nee_clean, h_clean, le_clean)

contains

    real(kind = dbl) function ErrToMissing(v)
        real(kind = dbl), intent(in) :: v
        if (v == error) then
            ErrToMissing = error
        else
            ErrToMissing = v
        end if
    end function ErrToMissing

    character(24) function FormatOrErr(v)
        real(kind = dbl), intent(in) :: v
        if (v == error) then
            FormatOrErr = trim(adjustl(EddyFlowProj%err_label))
        else
            write(FormatOrErr, '(es14.6e2)') v
            FormatOrErr = adjustl(FormatOrErr)
        end if
    end function FormatOrErr

    character(4) function IntToStr(v)
        integer, intent(in) :: v
        write(IntToStr, '(i0)') v
    end function IntToStr

end subroutine PostProcessFluxDespiking
