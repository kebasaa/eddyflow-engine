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
    logical, external :: GasHasSpectralFit
    integer, external :: PrimaryWaterOutSlot
    integer, external :: SlotFromSpectralStamp
    !> local variables
    integer :: gas
    integer :: cls
    integer :: wsl
    integer :: i
    integer :: open_status
    integer :: read_status
    integer :: slot
    !> The hygrometer the unnamed primary table belongs to, and the header line
    !> that says so.
    integer :: water_slot
    character(ShortInstringLen) :: water_header
    !> A hygrometer block's `exp=` token, and whether it parsed.
    character(96) :: exp_text
    integer :: exp_status
    logical :: short_file
    !> Gases this file actually carries a fit for.
    integer :: n_fitted
    integer :: n_configured
    real(kind = dbl) :: skipFn, skipfc
    !> One block's rows as they sit in the file: indexed by MONTH.
    real(kind = dbl) :: monthFn(12), monthfc(12)
    character(64) :: sa_tags(GHGNumVar)
    character(64) :: blockname
    character(ShortInstringLen) :: dataline


    !> Open planar fit file and read rotation matrices
    call LogSay(' Reading spectral assessment file: ')
    write(*,'(a)') '  ' // AuxFile%sa(1:len_trim(AuxFile%sa))
    write(ulog,'(a)') '  ' // AuxFile%sa(1:len_trim(AuxFile%sa))
    open(udf, file = AuxFile%sa, status = 'old', iostat = open_status)

    RegPar%Fn = 0d0
    RegPar%fc = 0d0
    RegPar%f2 = 0d0
    read_status = 0
    short_file = .false.
    if (open_status == 0) then
        call LogSay('  Spectral assessment file found, importing content..')
        !> Six free-text lines, then the primary water table's own header,
        !> which is kept rather than skipped: it carries the stamp saying
        !> which hygrometer the nine rows below it were fitted from.
        do i = 1, 6
            read(udf, *)
        end do
        read(udf, '(a)', iostat = read_status) water_header
        if (read_status /= 0) water_header = ''
        !> Read H2O transfer functions for IIR filter.
        !>
        !> Into the slot the writer put them in. OutputSpectralAssessmentResults
        !> resolves PrimaryWaterOutSlot and writes RegPar(wsl, cls); this read
        !> them back from the second slot unconditionally, so a project whose
        !> water is record three wrote its RH regression into slot seven and
        !> read it into slot six - the water block landing on whatever gas
        !> record two happened to hold, and the real hygrometer left with no
        !> RH-dependent cutoff at all.
        wsl = PrimaryWaterOutSlot()
        !> Onto the hygrometer the file says it came from, when it says.
        !>
        !> This table is identified by its position, so it landed on whichever
        !> record this project calls primary - which need not be the one that
        !> was primary when the file was written. Two hygrometers and a
        !> re-ordered project is enough to put one analyser's RH-dependent
        !> cut-offs on the other.
        water_slot = SlotFromSpectralStamp(water_header)
        if (water_slot < firstGas) water_slot = wsl
        if (.not. GasSlotIsWater(water_slot)) water_slot = wsl
        do cls = RH10, RH90
            read(udf, '(a)') dataline
            dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
            read(dataline, *)  RegPar(water_slot, cls)%Fn, &
                RegPar(water_slot, cls)%fc
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
        !>
        !> This range is CLASSES, not months. The two are both 1..12 and were
        !> both spelled JAN:DEC, which is how a month index came to be stored
        !> where a class index belongs.
        do gas = firstGas, lastGas
            if (gas == wsl) cycle
            RegPar(gas, 1:MaxGasClasses)%Fn = error
            RegPar(gas, 1:MaxGasClasses)%fc = error
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
            !> Files written before the second hygrometer's block was named
            !> correctly say '<TAG> VAPOUR'. The word was never part of the
            !> tag, so it matched nothing and the block was discarded; drop it
            !> and those files resolve as they were always meant to.
            if (len_trim(blockname) > 7) then
                if (blockname(len_trim(blockname) - 6:len_trim(blockname)) &
                    == ' VAPOUR') &
                    blockname = blockname(1:len_trim(blockname) - 7)
            end if

            !> What the block says it is for wins over where it sits among
            !> repeats of its species. A block naming a species and an analyser
            !> this project measures goes to that record whatever the project's
            !> record order is now; the name match below is what reads every
            !> file written before the stamp existed.
            slot = SlotFromSpectralStamp(dataline)
            if (slot == water_slot) slot = 0
            if (slot == 0) then
                do gas = firstGas, lastGas
                    !> Only the slot the unnamed table above filled is
                    !> excluded - not "the primary", which is the same thing
                    !> only while the file agrees with this project about which
                    !> hygrometer that is. Every other hygrometer has a named
                    !> block like any gas, which is what lets its fit survive a
                    !> round trip at all.
                    if (gas == water_slot) cycle
                    if (len_trim(sa_tags(gas)) == 0) cycle
                    if (trim(adjustl(sa_tags(gas))) == trim(blockname)) then
                        slot = gas
                        exit
                    end if
                end do
            end if

            !> A hygrometer's block carries nine RH classes, a gas's
            !> twelve months. The header says which: `numerosity` is
            !> the count column only the RH tables have, and it has
            !> been in this format since before the records.
            if (index(dataline, 'numerosity') /= 0) then
                !> A hygrometer's own RH/cut-off coefficients, where the block
                !> states them. Without this the fit is read back only as nine
                !> class cut-offs, and the exponential the correction actually
                !> evaluates stays the primary's.
                if (slot > 0) then
                    call SpectralStampToken(dataline, 'exp=', exp_text)
                    if (len_trim(exp_text) > 0) then
                        do i = 1, len_trim(exp_text)
                            if (exp_text(i:i) == ',') exp_text(i:i) = ' '
                        end do
                        read(exp_text, *, iostat = exp_status) &
                            RegPar(slot, dum)%e1, RegPar(slot, dum)%e2, &
                            RegPar(slot, dum)%e3
                        if (exp_status /= 0) then
                            RegPar(slot, dum)%e1 = error
                            RegPar(slot, dum)%e2 = error
                            RegPar(slot, dum)%e3 = error
                        end if
                    end if
                end if
                do cls = RH10, RH90
                    read(udf, '(a)', iostat = read_status) dataline
                    if (read_status /= 0) exit
                    dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
                    if (slot > 0) then
                        read(dataline, *, iostat = read_status) &
                            RegPar(slot, cls)%Fn, RegPar(slot, cls)%fc
                    else
                        read(dataline, *, iostat = read_status) skipFn, skipfc
                    end if
                    if (read_status /= 0) exit
                end do
                if (read_status /= 0) exit
                cycle
            end if

            !> Read the twelve months, then map them onto this project's
            !> classes. The file is keyed by month and RegPar by class; storing
            !> a row straight into RegPar(slot, cls) conflated the two, which
            !> is right only for a single all-months group.
            monthFn = error
            monthfc = error
            do cls = JAN, DEC
                read(udf, '(a)', iostat = read_status) dataline
                if (read_status /= 0) exit
                dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
                if (slot > 0) then
                    read(dataline, *, iostat = read_status) &
                        monthFn(cls), monthfc(cls)
                else
                    !> A gas this project does not carry. Consume the block so
                    !> the file stays aligned rather than skipping it.
                    read(dataline, *, iostat = read_status) skipFn, skipfc
                end if
                if (read_status /= 0) exit
            end do
            if (read_status /= 0) exit
            if (slot > 0) call MonthlyRegParToClasses(slot, monthFn, monthfc)
        end do

        !> Short means a gas this project wants got no block, which is now
        !> answered by what was found rather than by how far the loop got.
        do gas = firstGas, lastGas
            if (GasSlotIsWater(gas)) cycle
            if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
            !> Classes again, not months - "no class of this gas got a fit".
            if (all(RegPar(gas, 1:MaxGasClasses)%fc == error)) short_file = .true.
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

        !> A file read from disk can be partial for the same reason a file
        !> written on the fly can: some gases were fitted and some were not.
        !> The writer says so as it writes; this is the other way in, and
        !> without it a stored partial assessment would mix in-situ and
        !> analytic corrections without a word.
        !> Counted over every configured gas, water included, rather than off
        !> short_file: that flag deliberately skips hygrometers, and water
        !> unfitted while the gases are fitted is exactly the common case.
        n_fitted = 0
        n_configured = 0
        do gas = firstGas, lastGas
            if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
            n_configured = n_configured + 1
            if (GasHasSpectralFit(gas)) n_fitted = n_fitted + 1
        end do
        if (n_fitted > 0 .and. n_fitted < n_configured) call ExceptionHandler(110)

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
        call LogSay(' Done.')
    else
        !> If the specified file was not found or is empty,
        !> switches to an analytic method
        EddyFlowProj%hf_meth = 'moncrieff_97'
        call ExceptionHandler(65)
    end if
end subroutine ReadSpectralAssessmentFile
