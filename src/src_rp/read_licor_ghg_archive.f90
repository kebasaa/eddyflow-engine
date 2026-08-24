!***************************************************************************
! read_licor_ghg_archive.f90
! --------------------------
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
! \brief       Read one LI-COR GHG raw file and store all data and metadata
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadLicorGhgArchive(ZipFile, FirstRecord, LastRecord, LocCol, &
    LocBypassCol, MetaIsNeeded, BiometIsNeeded, DataIsNeeded, ValidateMetadata, &
    fRaw, nrow, ncol, skip_file, passed, faulty_col, N, FileEndReached, printout, &
    NextZipFile)

    use m_rp_global_var
    use m_ghg_prefetch
    implicit none
    !> in/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow, ncol
    character(*), intent(in) :: ZipFile
    !> The archive most likely wanted next, so its extraction can be under
    !> way while this one's data is being turned into a flux. Empty when
    !> there is no next - the preamble reads one file and stops. Required
    !> rather than optional for the same reason as UnZipArchive's arguments:
    !> no explicit interface exists for an external procedure.
    character(*), intent(in) :: NextZipFile
    logical, intent(in) :: MetaIsNeeded
    logical, intent(in) :: BiometIsNeeded
    logical, intent(in) :: DataIsNeeded
    logical, intent(in) :: ValidateMetadata
    logical, intent(in) :: printout
    integer, intent(out) :: N
    integer, intent(out) :: faulty_col
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: skip_file
    logical, intent(out) :: passed(32)
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    type(ColType), intent(inout) :: LocBypassCol(MaxNumCol)
    logical, intent(out) :: FileEndReached
    !> local variables
    integer :: del_status
    character(PathLen) :: MetaFile
    character(PathLen) :: DataFile
    character(PathLen) :: BiometFile
    character(PathLen) :: BiometMetaFile
    character(CommLen) :: comm
    logical :: skip_biomet_file
    logical :: prefetched
    character(PathLen) :: SrcDir


    skip_file = .false.
    passed = .true.

    !> Unzip archive, unless something already did it for this one. A claim
    !> is granted only for this exact archive and only once the extraction
    !> has signalled that it finished, so failing to get one costs nothing
    !> but the ordinary extraction.
    call GhgPrefetchClaim(ZipFile, prefetched)
    if (prefetched) then
        SrcDir = GhgPrefetchDir()
    else
        SrcDir = trim(adjustl(TmpDir))
    end if
    call UnZipArchive(ZipFile, 'metadata','data', MetaFile, DataFile, &
        BiometFile, BiometMetaFile, skip_file, SrcDir, prefetched)
    if (skip_file) then
        call GhgPrefetchRelease()
        return
    end if

    if (MetaFile /= 'none') &
        MetaFile = trim(adjustl(SrcDir)) // trim(Metafile)
    if (DataFile /= 'none') &
        DataFile = trim(adjustl(SrcDir)) // trim(DataFile)
    if (BiometMetaFile /= 'none') &
        BiometMetaFile = trim(adjustl(SrcDir)) // trim(BiometMetaFile)
    if (BiometFile /= 'none') &
        BiometFile = trim(adjustl(SrcDir)) // trim(BiometFile)

    !> First, handle biomet data and metadata files
    if (BiometIsNeeded) then
        if (BiometFile == 'none' .or. BiometMetaFile == 'none') then
            fnbRecs = 0
            call ExceptionHandler(5)
        else
            call ReadBiometMetaFile(BiometMetaFile, skip_biomet_file)
            if (.not. skip_biomet_file) &
                call ReadBiometFile(BiometFile, skip_biomet_file)
            if (skip_biomet_file) then
                fnbRecs = 0
                call ExceptionHandler(44)
            end if
        end if
    end if

    !> Handle metadata file
    if (MetaIsNeeded) then
        if (MetaFile == 'none') then
            call ExceptionHandler(3)
            skip_file = .true.
            call GhgPrefetchRelease()
            return
        end if
        call ReadMetadataFile(LocCol, MetaFile, skip_file, printout)
        if (skip_file) then
            call GhgPrefetchRelease()
            return
        end if
        if (DataIsNeeded) then

            !> If it's in the raw file processing loop, define used variables
            !> based on variables already identified (LocBypassCol)
            call RetrieveVarsSelection(LocBypassCol, LocCol)
            !> Re-apply user column selection: ReadMetadataFile resets var names
            !> (e.g. col_ts reverts to 'fast_t'), and BypassCol is never populated
            !> for GHG files, so RetrieveVarsSelection alone cannot restore useit.
            call DefineUsedVariables(LocCol)
        else
            !> In the preamble phase
            !> Embedded mode: define variables to be used,
            !> based on availability in the metadata file
            if (EddyFlowProj%run_env == 'embedded' &
                .and. EddyFlowProj%run_mode == 'express') &
                call DefaultVarsSelection(LocCol)
            !> Desktop mode: define used variables based on user selection
            !> at processing-time from GUI
            call DefineUsedVariables(LocCol)
        end if
        if (ValidateMetadata) then
            call MetadataFileValidation(Col, passed, faulty_col)
            if (.not. passed(1)) then
                call GhgPrefetchRelease()
                return
            end if
        else
            passed(1) = .true.
        end if
    end if

    !> Handle raw data file
    if (DataIsNeeded) then
        if (DataFile == 'none') then
            call ExceptionHandler(4)
            skip_file = .true.
            call GhgPrefetchRelease()
            return
        end if
        call ImportNativeData(DataFile, FirstRecord, LastRecord, &
            LocCol, fRaw, size(fRaw, 1), size(fRaw, 2), &
            skip_file, N, FileEndReached)
        if (skip_file) then
            call GhgPrefetchRelease()
            return
        end if
    end if

    !> Delete data and metadata files
    comm = (trim(comm_del) // ' ' // DataFile(1:len_trim(DataFile)) // ' ' &
        // trim(adjustl(MetaFile)) // ' ' &
        // trim(adjustl(BiometFile)) // ' ' &
        // trim(adjustl(BiometMetaFile)) &
        // ' *.status ' // comm_err_redirect)

    del_status = system(trim(comm))

    !> Nothing is reading the prefetch directory any more, so the next
    !> archive may be extracted into it - and this is the moment to ask,
    !> because everything the caller does from here is computation this
    !> waits behind otherwise.
    call GhgPrefetchRelease()
    call GhgPrefetchStart(NextZipFile)
end subroutine ReadLicorGhgArchive
