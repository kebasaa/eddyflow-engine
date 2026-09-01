!***************************************************************************
! biomet_retrieve_embedded_data.f90
! ---------------------------------
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
! \brief       Finalize retrieval of biomet data from already created bSet \n
!              in case of embedded biomet files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometRetrieveEmbeddedData(proceed, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: proceed
    logical, intent(in) :: printout
    !> local variables
    integer :: i
    !> Said once per run, not once per averaging period: the project states
    !> these column numbers once, so a wrong one is wrong every period.
    logical, save :: slot_warned = .false.

    !> Retrieve embedded biomet data if they exist (the option was
    !> selected and data was successfully read with at least one
    !> valid biomet record)
    if (printout) write(*,'(a)') '  Retrieving biomet data..'
    if (printout) write(ulog,'(a)') '  Retrieving biomet data..'

    !> Initialize biomet data to error
    if (allocated(bAggr)) bAggr = error
    if (allocated(bAggrFluxnet)) bAggrFluxnet = error
    if (allocated(bAggrEddyFlow)) bAggrEddyFlow = error

    if (proceed) then
        if (printout) write(LogInteger, '(i3)') nbRecs
        if (printout) write(*, '(a)') '   ' // trim(adjustl(LogInteger)) &
            // ' biomet records imported.'
        if (printout) write(ulog, '(a)') '   ' // trim(adjustl(LogInteger)) &
            // ' biomet records imported.'

        !> Aggregate biomet variables over the averaging interval
        call BiometAggregate(bSet, size(bSet, 1), size(bSet, 2), bAggr)

        !> Convert data to standard units
        call BiometStandardEddyFlowUnits()

        !> Aggregate biomet variables over the averaging interval
        call BiometAggregate(bSet, size(bSet, 1), size(bSet, 2), bAggrEddyFlow)

        !> Convert aggregated values to FLUXNET units
        call BiometStandardFluxnetUnits()
    else
        if (printout) call ExceptionHandler(72)
    end if

    !> Associate values to variables, as selected by user
    !> The - 2 is to account for the DATE and TIME columns in the file, which
    !> are not included in bAggr. The 2 shall eventually be replaced by nbTimestamp
    !> as per read_biomet_meta_file.f90
    !>
    !> The bounds test is not defensive padding: biom_ta and its siblings are
    !> column NUMBERS a project states by hand, and nothing has checked them
    !> against the file they are supposed to index. A project written for one
    !> site and pointed at another - biom_ta=2 against a biomet file whose
    !> second column is TIME - lands on bAggrEddyFlow(0), one element before
    !> the array. A release build reads whatever is there and reports it as a
    !> temperature; only -fcheck=all makes it stop. Out-of-range slots are
    !> left at their error value, which is what an unstated slot already gets.
    do i = bTa, bRg
        if (bSetup%sel(i) <= 0) cycle
        if (bSetup%sel(i) - 2 < 1 .or. bSetup%sel(i) - 2 > size(bAggrEddyFlow)) then
            if (.not. slot_warned) then
                call LogSay('   Warning> A biomet column selected in the project ' &
                    // '(biom_ta and its siblings) is outside the range the biomet ' &
                    // 'file describes, so that variable is not used. Check those ' &
                    // 'column numbers against this file''s own columns.')
                slot_warned = .true.
            end if
            cycle
        end if
        biomet%val(i) = bAggrEddyFlow(bSetup%sel(i) - 2)
    end do

    !> Deallocate variables no longer used
    if (allocated(bSet)) deallocate(bSet)
    if (allocated(bTs)) deallocate(bTs)
    if (printout) write(*,'(a)') '  Done.'
    if (printout) write(ulog,'(a)') '  Done.'

end subroutine BiometRetrieveEmbeddedData



