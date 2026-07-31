!*******************************************************************************
! report_imported_spectra.f90
! ---------------------------
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
! \brief       Output on stdout number of spectra and co-spectra imported
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine ReportImportedSpectra(nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nbins
    !> Local variables
    integer :: gas
    integer :: cnt
    character(16) :: name
    character(16), external :: GasName
    logical, external :: GasSlotIsWater

    write(*, '(a)') '  Imported gas binned spectra:'
    !> One line per configured gas, named from its record. This used to be
    !> four fixed lines labelled CO2/H2O/CH4/Gas 4, which on any project
    !> whose records are ordered differently named the wrong species, and
    !> which said nothing at all about gases past the fourth.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        !> Water is binned by relative humidity class, everything else by
        !> month - the same split SpectraSortingAndAveraging applies, so the
        !> count has to be summed the same way or a hygrometer reads as zero.
        if (GasSlotIsWater(gas)) then
            cnt = sum(MeanBinSpec(nbins/2, RH10:RH90)%cnt(gas))
        else
            cnt = MeanBinSpec(nbins/2, 1)%cnt(gas)
        end if
        name = GasName(gas)
        write(*, '(a, i5)') '   ' // trim(name) // ': ', cnt
    end do
    write(*, '(a)') '  Done.'
end subroutine ReportImportedSpectra
