!***************************************************************************
! read_timelag_opt_file.f90
! -------------------------
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
! \brief       Read time-lag optimization file and import relevant parameters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadTimelagOptFile(ncls)
    use m_rp_global_var
    use m_pwb_timelag, only: TimelagOptGasLabel
    implicit none
    integer, external :: PrimaryWaterSlot
    !> in/out variables
    integer :: ncls
    !> local variables
    integer :: open_status
    integer :: read_status
    integer :: gas
    logical :: matched
    character(32) :: gname
    logical, external :: GasSlotIsWater
    character(500) :: strg


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading time-lag optimization file: ' // AuxFile%to(1:len_trim(AuxFile%to))
    open(udf, file = AuxFile%to, status = 'old', iostat = open_status)

    if (open_status == 0) then
        write(*, '(a)') ' Time lag optimization file found, retrieving content..'
        do
            read(udf, '(a)', iostat = read_status) strg
            if (read_status /= 0) exit

            !> One lookup per configured gas, under the same names
            !> TimelagOptGasLabel gave the writer. These were four hand-written
            !> triples, so a gas past the fourth was written by the optimiser
            !> and then never read back - its window was lost between runs and
            !> it fell to its nominal one without saying so.
            !>
            !> List-directed, not '(f6.2)': the sentinel a gas with no
            !> determinations is written as needs nine characters, and a fixed
            !> f6.2 read of it returns garbage rather than failing. This also
            !> reads a file written by the previous, narrower format.
            matched = .false.
            do gas = firstGas, lastGas
                gname = TimelagOptGasLabel(gas)
                if (index(strg, 'Median_' // trim(gname) // '_timelag_[s]') /= 0) then
                    read(strg(index(strg, ':')+1:len_trim(strg)), *) toPasGas(gas)%def
                    !> A plain median row for the PRIMARY means it was not
                    !> classed by RH. A second hygrometer always has one and
                    !> says nothing about whether the primary was classed.
                    if (gas == PrimaryWaterSlot()) ncls = 0
                    matched = .true.
                    exit
                end if
                if (index(strg, 'Mimimum_' // trim(gname) // '_timelag_[s]') /= 0) then
                    read(strg(index(strg, ':')+1:len_trim(strg)), *) toPasGas(gas)%min
                    matched = .true.
                    exit
                end if
                if (index(strg, 'Maximum_' // trim(gname) // '_timelag_[s]') /= 0) then
                    read(strg(index(strg, ':')+1:len_trim(strg)), *) toPasGas(gas)%max
                    matched = .true.
                    exit
                end if
            end do
            if (matched) cycle

            !> h2o as a function of RH
            ncls = 0
            if (index(strg, 'H2O_timelag_determinations_as_a_function') /= 0) then
                !> Skip one line
                read(udf, *)
                read(udf, *)
                !> Read as many classes as available
                do
                    read(udf, '(a)', iostat = read_status) strg
                    if (read_status /= 0 .or. index(strg, '%') == 0) exit
                    ncls = ncls + 1
                    read(strg(20:len_trim(strg)), '(3(f14.2))') &
                        toH2O(ncls)%def,  toH2O(ncls)%min,  toH2O(ncls)%max
                end do
                exit
            end if
        end do
    else
       !> If the specified file is not found or is empty, switches to covariance maximization without default
        Meth%tlag = 'maxcov'
        call ExceptionHandler(39)
    end if
    write(*,'(a)')   ' Done.'
end subroutine ReadTimelagOptFile
