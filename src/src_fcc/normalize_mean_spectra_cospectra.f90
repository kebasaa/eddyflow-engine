!***************************************************************************
! normalize_mean_spectra_cospectra.f90
! ------------------------------------
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
! \brief       Normalize sums to obtain average spectra and cospectra
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NormalizeMeanSpectraCospectra(nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nbins
    !> Local variables
    integer :: bin
    integer :: cls


    !> Normalize to complete the averaging process
    do cls = 1, MaxGasClasses
        do bin = 1, nbins
            !> Sorted spectra
            where (MeanBinSpec(bin, cls)%cnt(firstGas:lastGas) >= FCCsetup%SA%min_smpl)
                MeanBinSpec(bin, cls)%fnum(firstGas:lastGas) = &
                    MeanBinSpec(bin, cls)%fnum(firstGas:lastGas) / MeanBinSpec(bin, cls)%cnt(firstGas:lastGas)
                MeanBinSpec(bin, cls)%fn(firstGas:lastGas) = &
                    MeanBinSpec(bin, cls)%fn(firstGas:lastGas)   / dfloat(MeanBinSpec(bin, cls)%cnt(firstGas:lastGas))
                MeanBinSpec(bin, cls)%of(firstGas:lastGas) = &
                    MeanBinSpec(bin, cls)%of(firstGas:lastGas)   / dfloat(MeanBinSpec(bin, cls)%cnt(firstGas:lastGas))
                MeanBinSpec(bin, cls)%ts(firstGas:lastGas) = &
                    MeanBinSpec(bin, cls)%ts(firstGas:lastGas)   / dfloat(MeanBinSpec(bin, cls)%cnt(firstGas:lastGas))
            elsewhere
                MeanBinSpec(bin, cls)%fnum(firstGas:lastGas) = nint(error)
                MeanBinSpec(bin, cls)%fn(firstGas:lastGas) = error
                MeanBinSpec(bin, cls)%of(firstGas:lastGas) = error
                MeanBinSpec(bin, cls)%ts(firstGas:lastGas) = error
            end where
            !> Sorted cospectra
            where (MeanBinCosp(bin, cls)%cnt(ts:lastGas) > 0)
                MeanBinCosp(bin, cls)%fnum(ts:lastGas) = &
                    MeanBinCosp(bin, cls)%fnum(ts:lastGas) / MeanBinCosp(bin, cls)%cnt(ts:lastGas)
                MeanBinCosp(bin, cls)%fn(ts:lastGas) = &
                    MeanBinCosp(bin, cls)%fn(ts:lastGas)   / dfloat(MeanBinCosp(bin, cls)%cnt(ts:lastGas))
                MeanBinCosp(bin, cls)%of(ts:lastGas) = &
                    MeanBinCosp(bin, cls)%of(ts:lastGas)   / dfloat(MeanBinCosp(bin, cls)%cnt(ts:lastGas))
            elsewhere
                MeanBinCosp(bin, cls)%fnum(ts:lastGas) = nint(error)
                MeanBinCosp(bin, cls)%fn(ts:lastGas) = error
                MeanBinCosp(bin, cls)%of(ts:lastGas) = error
            end where
        end do
    end do
end subroutine NormalizeMeanSpectraCospectra
