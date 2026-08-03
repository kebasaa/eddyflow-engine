!***************************************************************************
! read_spectral_assessment_file.f90
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
! \brief       Read file containing results of spectral assessment to \n
!              retrieve cutoff frequencies
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadSpectralAssessmentFile()
    use m_fx_global_var
    implicit none
    logical, external :: GasSlotIsWater
    !> local variables
    integer :: gas
    integer :: cls
    integer :: i
    integer :: open_status
    integer :: read_status
    integer :: slot
    logical :: short_file
    real(kind = dbl) :: skipFn, skipfc
    character(64) :: sa_tags(GHGNumVar)
    character(64) :: blockname
    character(ShortInstringLen) :: dataline


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading spectral assessment file: '
    write(*,'(a)') '  ' // AuxFile%sa(1:len_trim(AuxFile%sa))
    open(udf, file = AuxFile%sa, status = 'old', iostat = open_status)

    RegPar%Fn = 0d0
    RegPar%fc = 0d0
    RegPar%f2 = 0d0
    read_status = 0
    short_file = .false.
    if (open_status == 0) then
        write(*, '(a)') '  Spectral assessment file found, importing content..'
        !> skip 7 lines
        do i = 1, 7
            read(udf, *)
        end do
        !> Read H2O transfer functions for IIR filter
        do cls = RH10, RH90
            read(udf, '(a)') dataline
            dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
            read(dataline, *)  RegPar(h2o, cls)%Fn, RegPar(h2o, cls)%fc
        end do

        !> One block per configured gas but water, matching what
        !> OutputSpectralAssessmentResults writes over the same range.
        !>
        !> A file written before this range widened carries only three blocks,
        !> and a file written for a smaller project carries fewer than this one
        !> expects. Neither may be read as though the blocks were there: the
        !> next thing in the file is the exponential-fit section, and a blind
        !> `read` of its title line succeeds, so a count mismatch would be
        !> consumed as transfer-function parameters rather than noticed.
        !>
        !> The block header is therefore checked, not skipped. On a mismatch
        !> the two peeked lines are put back and the loop stops, leaving the
        !> file positioned exactly where a full-length read would have left it
        !> - so the sections below still parse and the remaining gases stay
        !> unfitted, falling back to an analytic transfer function.
        !> Start every gas block unfitted. RegPar is zeroed above, and a
        !> cut-off of zero is not a missing value - it is an infinitely
        !> aggressive correction, which is what a gas absent from the file
        !> would silently receive. The readiness check keys on `error`, so
        !> writing the sentinel here is what makes an absent gas fall back to
        !> the analytic method instead of inventing a correction for it.
        !> Water is left alone: its classes are the RH table read above.
        do gas = firstGas, lastGas
            if (GasSlotIsWater(gas)) cycle
            RegPar(gas, JAN:DEC)%Fn = error
            RegPar(gas, JAN:DEC)%fc = error
        end do

        !> Blocks are matched by the name in their own header, not by their
        !> position in the file.
        !>
        !> This loop used to walk the expected gas list and assign the Nth
        !> block to the Nth slot, testing only that the header contained the
        !> word TFP - the gas name it carries was never read. So a file whose
        !> block set differed from what this project expects, by even one
        !> gas, had every block after the difference silently assigned to the
        !> wrong species. That is not hypothetical: it is what stops the water
        !> carve-out here from being widened to a second hygrometer, because
        !> every file written so far contains a block for one.
        !>
        !> Driven by the file now: read blocks until the headers stop, resolve
        !> each name through SpectralGasNames - the same helper the writer
        !> names them with - and consume a block for a gas this project does
        !> not have rather than mis-assigning it.
        call SpectralGasNames(sa_tags)
        do gas = firstGas, lastGas
            call uppercase(sa_tags(gas))
        end do

        do
            !> Blank separator, then the block header naming the gas
            read(udf, *, iostat = read_status)
            if (read_status /= 0) exit
            read(udf, '(a)', iostat = read_status) dataline
            if (read_status /= 0) exit
            if (index(dataline, 'TFP') == 0) then
                !> End of the per-gas section. Put both lines back: what
                !> follows is counted from here.
                backspace(udf)
                backspace(udf)
                exit
            end if

            blockname = adjustl(dataline(1:index(dataline, 'TFP') - 1))
            call uppercase(blockname)
            slot = 0
            do gas = firstGas, lastGas
                if (GasSlotIsWater(gas)) cycle
                if (len_trim(sa_tags(gas)) == 0) cycle
                if (trim(adjustl(sa_tags(gas))) == trim(blockname)) then
                    slot = gas
                    exit
                end if
            end do

            do cls = JAN, DEC
                read(udf, '(a)', iostat = read_status) dataline
                if (read_status /= 0) exit
                dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
                if (slot > 0) then
                    read(dataline, *, iostat = read_status) &
                        RegPar(slot, cls)%Fn, RegPar(slot, cls)%fc
                else
                    !> A gas this project does not carry. Consume the block so
                    !> the file stays aligned rather than skipping it.
                    read(dataline, *, iostat = read_status) skipFn, skipfc
                end if
                if (read_status /= 0) exit
            end do
            if (read_status /= 0) exit
        end do

        !> Short means a gas this project wants got no block, which is now
        !> answered by what was found rather than by how far the loop got.
        do gas = firstGas, lastGas
            if (GasSlotIsWater(gas)) cycle
            if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
            if (all(RegPar(gas, JAN:DEC)%fc == error)) short_file = .true.
        end do

        !> A truncated block is a malformed file, not an older one: the
        !> position is no longer trustworthy, so stop rather than parse on.
        if (read_status /= 0) then
            close(udf)
            EddyFlowProj%hf_meth = 'moncrieff_97'
            call ExceptionHandler(65)
            return
        end if
        if (short_file) call ExceptionHandler(65)

        !> skip 4 lines
        do i = 1, 4
            read(udf, *)
        end do

        !> Read parameters of exponential fit fc vs. RH
        read(udf, *) RegPar(dum, dum)%e1, RegPar(dum, dum)%e2, RegPar(dum, dum)%e3

        !> skip 6 lines
        do i = 1, 6
            read(udf, *)
        end do
        !> Read parameters of Ibrom's model for spectral correction factor
        read(udf, '(a)') dataline
        dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
        read(dataline, *)  UnPar(1), UnPar(2)
        read(udf, '(a)') dataline
        dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
        read(dataline, *)  StPar(1), StPar(2)
        close(udf)
        write(*, '(a)') ' Done.'
    else
        !> If the specified file was not found or is empty,
        !> switches to an analytic method
        EddyFlowProj%hf_meth = 'moncrieff_97'
        call ExceptionHandler(65)
    end if
end subroutine ReadSpectralAssessmentFile
