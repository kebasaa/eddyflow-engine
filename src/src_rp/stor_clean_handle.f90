!***************************************************************************
! stor_clean_handle.f90
! ----------------------
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
! \brief       Storage-flux cleaning (test_stor_clean): accumulate every
!              configured gas's storage term across the whole RP run, then
!              clean each series once the run is complete, writing a
!              dedicated output CSV.
!
! \details     Mirrors pfd_handle.f90's shape (a growable per-period cache,
!              populated during the main loop, consumed once after it), but
!              lives entirely in RP: the storage term is computed here and
!              never reaches FCC's ex record, unlike NEE/H/LE. Off by
!              default; StoreStorCache is only called when Test%stor_clean
!              is set, so a run with the feature off never allocates the
!              cache or writes the new file.
!
!              The actual cleaning arithmetic is m_storage_clean_core.f90 -
!              this file is only the RP-side plumbing: accumulate, sort,
!              derive the averaging period's periods-per-day, call it once
!              per configured gas, write the result.
!
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          m_storage_clean_core.f90, pfd_handle.f90, period_ordering.f90
!***************************************************************************

!> Append one period's per-gas storage terms to the whole-run cache. Called
!> from the RP flux loop, guarded by Test%stor_clean at the call site.
subroutine StoreStorCache(date, time)
    use m_rp_global_var
    implicit none
    character(*), intent(in) :: date, time
    type(StorCacheEntryType), allocatable :: tmp(:)

    allocate(tmp(StorCacheN + 1))
    if (StorCacheN > 0) tmp(1:StorCacheN) = StorCache(1:StorCacheN)
    tmp(StorCacheN + 1)%date = date
    tmp(StorCacheN + 1)%time = time
    tmp(StorCacheN + 1)%strg = Stor%of
    call move_alloc(tmp, StorCache)
    StorCacheN = StorCacheN + 1
end subroutine StoreStorCache


!***************************************************************************
!> Runs once, after the whole RP run has populated StorCache. Sorts
!> chronologically, derives periods-per-day from the modal spacing between
!> consecutive periods, cleans every configured gas's storage series, and
!> writes one CSV with the result.
!***************************************************************************
subroutine PostProcessStorClean()
    use m_rp_global_var
    use m_storage_clean_core, only: CleanStorageSeries
    implicit none
    integer, allocatable :: ord(:), hod(:)
    integer(8), allocatable :: tmin(:)
    real(kind = dbl), allocatable :: strg(:, :), strg_clean(:, :)
    integer, allocatable :: strg_spike(:, :)
    integer :: i, gas, mfreq, stor_unit, ngas
    integer(8) :: modal_step, phase
    character(PathLen) :: StorPath, Test_Path
    integer :: dot, open_status
    integer(8), external :: PeriodMinutes, ModalPeriodStep
    character(32), external :: GasOutputLabel
    character(LongOutstringLen) :: hdr
    character(32) :: species

    if (StorCacheN < 20) then
        !> Too short for a whole-run, time-of-day-binned outlier test to
        !> mean anything - RFlux's own cleanFlux() warns below ten days of
        !> half-hourly data for the same reason. Say why nothing was
        !> written rather than fail or silently skip.
        call LogSay('   Storage-flux cleaning: run has fewer than 20 periods, skipped.')
        return
    end if

    allocate(ord(StorCacheN), tmin(StorCacheN))
    do i = 1, StorCacheN
        ord(i) = i
        tmin(i) = PeriodMinutes(StorCache(i)%date, StorCache(i)%time)
    end do
    call SortPeriodsByMinutes(ord, tmin, StorCacheN)

    !> Periods per day, from the most common spacing between consecutive
    !> sorted periods rather than the first pair, so one irregular gap at
    !> the start of the run cannot mis-set it.
    modal_step = ModalPeriodStep(tmin, ord, StorCacheN)
    if (modal_step <= 0_8) then
        call LogSay('   Storage-flux cleaning: could not determine the averaging period, skipped.')
        return
    end if
    mfreq = nint(1440d0 / dble(modal_step))
    if (mfreq < 4) then
        call LogSay('   Storage-flux cleaning: fewer than 4 periods/day, skipped.')
        return
    end if

    !> Time-of-day class per period, from its own minutes-since-midnight
    !> divided by the period length - not from array position, so a gap in
    !> the run cannot shift the phase.
    allocate(hod(StorCacheN))
    do i = 1, StorCacheN
        phase = mod(tmin(ord(i)), 1440_8)
        hod(i) = int(phase / modal_step) + 1
        if (hod(i) < 1) hod(i) = 1
        if (hod(i) > mfreq) hod(i) = mfreq
    end do

    allocate(strg(StorCacheN, GHGNumVar), strg_clean(StorCacheN, GHGNumVar))
    allocate(strg_spike(StorCacheN, GHGNumVar))
    ngas = 0
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        ngas = ngas + 1
        do i = 1, StorCacheN
            strg(i, gas) = StorCache(ord(i))%strg(gas)
        end do
        call CleanStorageSeries(strg(:, gas), StorCacheN, hod, mfreq, error, &
            strg_spike(:, gas), strg_clean(:, gas))
    end do

    if (ngas == 0) then
        call LogSay('   Storage-flux cleaning: no gas configured, skipped.')
        return
    end if

    !> Same construction as FLUXNET_Path in init_fluxnet_file_rp.f90, with
    !> its own padding, so the two live side by side and neither depends on
    !> the other's open file handle.
    Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
              // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
              // STOR_FilePadding // Timestamp_FilePadding // CsvExt
    dot = index(Test_Path, CsvExt, .true.) - 1
    StorPath = Test_Path(1:dot) // CsvExt

    open(newunit = stor_unit, file = trim(StorPath), status = 'replace', iostat = open_status, encoding = 'utf-8')
    if (open_status /= 0) then
        call LogSay('   Storage-flux cleaning: could not open output file, skipped.')
        return
    end if

    hdr = 'TIMESTAMP'
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        species = GasOutputLabel(gas)
        call lowercase(species)
        hdr = trim(hdr) // ',' // trim(species) // '_strg' &
            // ',' // trim(species) // '_strg_spike' &
            // ',' // trim(species) // '_strg_cleaned'
    end do
    write(stor_unit, '(a)') trim(hdr)

    do i = 1, StorCacheN
        write(stor_unit, '(a)', advance = 'no') &
            trim(StorCache(ord(i))%date) // ' ' // trim(StorCache(ord(i))%time)
        do gas = firstGas, lastGas
            if (.not. E2Col(gas)%present) cycle
            write(stor_unit, '(a)', advance = 'no') ',' // trim(FormatOrErr(strg(i, gas))) &
                // ',' // trim(IntToStr(strg_spike(i, gas))) &
                // ',' // trim(FormatOrErr(strg_clean(i, gas)))
        end do
        write(stor_unit, '(a)') ''
    end do
    close(stor_unit)

    deallocate(ord, tmin, hod, strg, strg_clean, strg_spike)

contains

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

end subroutine PostProcessStorClean
