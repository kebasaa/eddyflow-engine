!***************************************************************************
! read_binned_file.f90
! --------------------
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
! \brief       Reads file containing binned co-spectra and import relevant ones
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadBinnedFile(InFile, BinSpec, BinCosp, nrow, nbins, skip)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(FileListType), intent(in) :: InFile
    type(SpectraSetType), intent(out) :: BinSpec(nrow)
    type(SpectraSetType), intent(out) :: BinCosp(nrow)
    logical, intent(out) :: skip
    integer, intent(out) :: nbins
    !> local variables
    integer :: i
    integer :: j
    integer :: var
    integer :: nvar
    integer :: open_status
    integer :: read_status
    integer :: spec_ord(GHGNumVar)
    integer :: cosp_ord(GHGNumVar)
    real(kind = dbl), allocatable :: aux(:)
    character(ShortInstringLen) :: dataline
    character(72) :: field
    character(72) :: probe
    character(72) :: speclabs(GHGNumVar)
    character(72) :: cosplabs(GHGNumVar)
    character(72) :: oldspeclabs(GHGNumVar)
    character(72) :: oldcosplabs(GHGNumVar)
    character(64) :: vartags(GHGNumVar)
    character(32), external :: LegacySpectralVarTag
    include '../src_common/interfaces_1.inc'

    !> Open file
    open(udf, file = InFile%path, iostat = open_status)
    !> Control on error in file opening
    skip = .false.
    if (open_status /= 0) then
        skip = .true.
        call ExceptionHandler(62)
        return
    end if

    BinSpec = ErrSpec
    BinCosp = ErrSpec

    !> The column names this file should carry, from the same helper the
    !> writer uses. This replaced a single positional read of exactly 18
    !> items - three frequencies, eight spectra, seven cospectra - which is
    !> why the on-the-fly spectral assessment could never see a fifth gas
    !> however wide the loops downstream were: its input had four.
    call SpectralVarTags(vartags)
    speclabs = ''
    cosplabs = ''
    oldspeclabs = ''
    oldcosplabs = ''
    do j = 1, GHGNumVar
        if (len_trim(vartags(j)) > 0) then
            speclabs(j) = 'f_nat*spec(' // trim(vartags(j)) // ')/var(' &
                // trim(vartags(j)) // ')'
            cosplabs(j) = 'f_nat*cospec(w_' // trim(vartags(j)) &
                // ')/cov(w_' // trim(vartags(j)) // ')'
            call uppercase(speclabs(j))
            call uppercase(cosplabs(j))
        end if
        !> A file written before the tags became record-derived spells the
        !> first four gas slots co2/h2o/ch4/gas4. Accepted as well, so an
        !> existing binned directory still imports.
        if (len_trim(LegacySpectralVarTag(j)) > 0) then
            oldspeclabs(j) = 'f_nat*spec(' // trim(LegacySpectralVarTag(j)) &
                // ')/var(' // trim(LegacySpectralVarTag(j)) // ')'
            oldcosplabs(j) = 'f_nat*cospec(w_' &
                // trim(LegacySpectralVarTag(j)) // ')/cov(w_' &
                // trim(LegacySpectralVarTag(j)) // ')'
            call uppercase(oldspeclabs(j))
            call uppercase(oldcosplabs(j))
        end if
    end do

    !> Skip to the column-name line. Counting header lines would break the
    !> moment the preamble changed length, and the previous reader discarded
    !> it only by letting the numeric read fail - which cannot tell a header
    !> from a corrupt row.
    dataline = ''
    do
        read(udf, '(a)', iostat = read_status) dataline
        if (read_status /= 0) exit
        if (index(dataline, '#_freq') > 0) exit
    end do
    if (read_status /= 0) then
        close(udf)
        skip = .true.
        call ExceptionHandler(62)
        return
    end if

    !> Match each column by name. A slot with no column keeps ErrSpec and so
    !> declines to be assessed, rather than picking up its neighbour's data.
    spec_ord = 0
    cosp_ord = 0
    nvar = 0
    do
        if (index(dataline, ',') /= 0) then
            field = dataline(1:index(dataline, ',') - 1)
        else
            field = dataline(1:len_trim(dataline))
        end if
        nvar = nvar + 1
        !> Case-insensitively. The writer this replaces named the fourth gas
        !> from SpecCol%label, which carries whatever case the metadata used
        !> - 'COS' where every other file says 'cos'. A binned directory kept
        !> from an earlier run is the normal case rather than the exception:
        !> the project points sa_bin_spectra at one.
        probe = field
        call uppercase(probe)
        do j = 1, GHGNumVar
            if (len_trim(speclabs(j)) > 0 .and. probe == speclabs(j)) &
                spec_ord(j) = nvar
            if (len_trim(cosplabs(j)) > 0 .and. probe == cosplabs(j)) &
                cosp_ord(j) = nvar
        end do
        do j = 1, GHGNumVar
            if (spec_ord(j) == 0 .and. len_trim(oldspeclabs(j)) > 0 &
                .and. probe == oldspeclabs(j)) spec_ord(j) = nvar
            if (cosp_ord(j) == 0 .and. len_trim(oldcosplabs(j)) > 0 &
                .and. probe == oldcosplabs(j)) cosp_ord(j) = nvar
        end do
        if (index(dataline, ',') == 0) exit
        dataline = dataline(index(dataline, ',') + 1:len_trim(dataline))
    end do

    if (nvar < 4) then
        close(udf)
        skip = .true.
        call ExceptionHandler(62)
        return
    end if

    allocate(aux(nvar))
    i = 0
    do
        i = i + 1
        if (i > nrow) exit
        read(udf, *, iostat = read_status) aux(1:nvar)
        if (read_status > 0) then
            i = i - 1
            cycle
        end if
        if (read_status < 0) exit
        BinSpec(i)%fnum  = nint(aux(1))
        BinSpec(i)%fn    = aux(2)
        BinSpec(i)%fnorm = aux(3)
        do j = 1, GHGNumVar
            if (spec_ord(j) > 0) BinSpec(i)%of(j) = aux(spec_ord(j))
            if (cosp_ord(j) > 0) BinCosp(i)%of(j) = aux(cosp_ord(j))
        end do
        BinCosp(i)%fnum = BinSpec(i)%fnum
        BinCosp(i)%fn   = BinSpec(i)%fn
        BinCosp(i)%fnorm = BinSpec(i)%fnorm
    end do
    nbins = i - 1
    close(udf)
    if (allocated(aux)) deallocate(aux)

    !> Un-normalize binned spectra by dividing by the frequency
    do var = w, lastGas
        where (BinSpec(1:nbins)%fn /= error &
            .and. BinSpec(1:nbins)%fn /= 0d0 .and. BinSpec(1:nbins)%of(var) /= error)
            BinSpec(1:nbins)%of(var) = BinSpec(1:nbins)%of(var) / BinSpec(1:nbins)%fn
        else where
            BinSpec(1:nbins)%of(var) = error
        end where
    end do

    !> Un-normalize binned cospectra by dividing by the frequency
    do var = ts, lastGas
        where (BinCosp(1:nbins)%fn /= error &
            .and. BinCosp(1:nbins)%fn /= 0d0 .and. BinCosp(1:nbins)%of(var) /= error)
            BinCosp(1:nbins)%of(var) = BinCosp(1:nbins)%of(var) / BinCosp(1:nbins)%fn
        else where
            BinCosp(1:nbins)%of(var) = error
        end where
    end do

    !> Check that spectra values are within reasonable values (not too high).
    !> Individually discard spectra if this is not the case
    !> This is somewhat arbitrary, introduced to eliminate observed implausible spectra
    !> It's very strict: one only outranged value will eliminate the whole (co)spectrum
    ol: do var = w, lastGas
        il: do i = 1, nbins
            if (BinSpec(i)%of(var) > MaxNormSpecValue) then
                BinSpec(:)%of(var) = error
                exit il
            end if
        end do il
    end do ol

    ol2: do var = ts, lastGas
        il2: do i = 1, nbins
            if (dabs(BinCosp(i)%of(var)) > MaxNormSpecValue) then
                BinCosp(:)%of(var) = error
                exit il2
            end if
        end do il2
    end do ol2

    !> Similar filter as above, but now imposes that f*spectrum < 1 for each frequency
    ol3: do var = w, lastGas
        il3: do i = 1, nbins
            if (BinSpec(i)%fn /= error .and. BinSpec(i)%of(var) /= error &
                .and. BinSpec(i)%fn * BinSpec(i)%of(var) > 1d0) then
                BinSpec(:)%of(var) = error
                exit il3
            end if
        end do il3
    end do ol3

    ol4: do var = ts, lastGas
        il4: do i = 1, nbins
            if (BinCosp(i)%fn /= error .and. BinCosp(i)%of(var) /= error .and. &
                dabs(BinCosp(i)%fn * BinCosp(i)%of(var)) > 10d0) then
                BinCosp(:)%of(var) = error
                exit il4
            end if
        end do il4
    end do ol4
end subroutine ReadBinnedFile
