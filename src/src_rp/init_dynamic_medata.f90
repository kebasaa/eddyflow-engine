!***************************************************************************
! init_dynamic_medata.f90
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
!*******************************************************************************
!
! \brief       Read dynamic metadata file and figure out available parameters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!*******************************************************************************
subroutine InitDynamicMetadata(N)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: N
    !> Local variables
    integer :: open_status
    integer :: io_status


    write(*, '(a)', advance = 'no') ' Initializing dynamic metadata usage..'

    !> Open file
    open(udf, file = AuxFile%DynMD, status = 'old', iostat = open_status)

    !> Interpret dynamic metadata file header and control in case of error
    if (open_status == 0) then
        call ReadDynamicMetadataHeader(udf)
    else
        call ExceptionHandler(68)
        EddyFlowProj%use_dynmd_file = .false.
    end if

    !> Count number of rows in the file (all of them, no matter if well formed),
    !> to give a maximum number to calibration data arrays
    N = 0
    countloop: do
        read(udf, *, iostat = io_status)
        if (io_status < 0 .or. io_status == 5001 .or. io_status == 5008) exit
        N = N + 1
    end do countloop
    close(udf)

    write(*, '(a)') ' Done.'
end subroutine InitDynamicMetadata

!***************************************************************************
!
! \brief       Reads and interprets file header, searching for knwon variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine ReadDynamicMetadataHeader(unt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    !> local variables
    character(LongInstringLen) :: dataline
    character(64) :: Headerlabels(NumStdDynMDVars)
    integer :: read_status
    integer :: sepa
    integer :: cnt
    integer :: i
    integer :: j
    integer :: slot
    integer :: nfields
    character(64) :: field_suffix(nDynMDGasFields)
    integer, external :: GasSlotFromDynMDTag


    read(unt, '(a)', iostat = read_status) dataline
    cnt = 0
    do
        sepa = index(dataline, ',')
        if (sepa == 0) sepa = len_trim(dataline) + 1
        if (len_trim(dataline) == 0) exit
        cnt = cnt + 1
        Headerlabels(cnt) = dataline(1:sepa - 1)
        dataline = dataline(sepa + 1: len_trim(dataline))
    end do

    DynamicMetadataOrder = nint(error)
    do i = 1, cnt
        do j = 1, NumStdDynMDVars
        if(trim(adjustl(StdDynMDVars(j))) &
            == trim(adjustl(Headerlabels(i)))) then
            DynamicMetadataOrder(j) = i
            exit
        end if
        end do
    end do

    !> The same header, resolved per gas slot.
    !>
    !> StdDynMDVars is a fixed list ending at gas4_irga_tau, so a project with
    !> more than four gases had no name it could give the fifth analyser -
    !> nothing in the file could reach it. Here every configured gas is
    !> offered `<label>_irga_*` under its own record label, so a COS record
    !> answers to `cos_irga_model`, and the four historical spellings keep
    !> working because GasSlotFromDynMDTag accepts them as aliases.
    !>
    !> Both passes fill in; for the first four gases they agree by
    !> construction, and the reader takes this one.
    call DynMDGasFieldNames(field_suffix, nfields)
    DynMDGasOrder = nint(error)
    do i = 1, cnt
        do j = 1, nfields
            slot = GasSlotFromDynMDTag(Headerlabels(i), field_suffix(j))
            if (slot < firstGas .or. slot > lastGas) cycle
            DynMDGasOrder(slot, j) = i
            exit
        end do
    end do
end subroutine ReadDynamicMetadataHeader
