!***************************************************************************
! gas_slot_resolution.f90
! ---------------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
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
! \brief       Resolving a gas slot: which species holds it, what it is called
!              and how its columns are scaled.
! \author      Jonathan Muller
! \note        Every routine here answers for whichever gas slot it is asked
!              about, from that slot's record. This is the layer the rest of
!              the engine goes through instead of assuming a species sits at a
!              fixed index.
!
!              Was gas4_output_units.f90, from when it held the fourth slot's
!              unit handling alone. It has not been that for a long time -
!              PrimaryWaterSlot, StatsLayoutSlots, SpectralVarTags and
!              GasSlotFromDynMDTag are none of them about a fourth gas or
!              about units - and a file whose name describes a fixed slot is
!              an odd place to keep the code that removes fixed slots.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!***************************************************************************
!
! \brief       Mole-basis scale of a gas's FLUXNET columns, by species.
! \author      Jonathan Muller
! \note        The FLUXNET row carries no units line, and the FP-In format
!              fixes the basis per species: CO2 in umol mol-1, H2O in
!              mmol mol-1, every other species in nmol mol-1. Internally the
!              engine holds every trace gas in umol mol-1 and water in
!              mmol mol-1 (ConvertTraceGasUnits), so these factors take the
!              internal value to the column's basis and nothing else.
!
!              This replaces a gain hard-coded for whatever gases occupied the
!              ch4 and gas4 slots. That named a position rather than a
!              species, so the same gas was reported in nmol mol-1 from slot 7
!              and in umol mol-1 from slot 9. The species comes from the gas
!              record, which is the authority everywhere else in this effort.
!
!              Deliberately NOT derived from the input unit: the raw signal was
!              already normalised on input, so scaling by the input factor
!              again would count it twice. The full output is the place where
!              the project's own unit is honoured, because it carries a units
!              row to declare it.
!***************************************************************************
real(kind = dbl) function FluxnetGasScale(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    integer :: rec
    character(32) :: species

    !> Anything not named by a record is a trace gas by default.
    FluxnetGasScale = 1d3
    rec = gas_slot - firstGas + 1
    if (rec < 1 .or. rec > min(EddyFlowProj%gas_num, MaxNumGases)) return

    species = EddyFlowProj%gas(rec)%var
    call uppercase(species)
    select case (trim(adjustl(species)))
        case ('CO2', 'H2O'); FluxnetGasScale = 1d0
        case default;        FluxnetGasScale = 1d3
    end select
end function FluxnetGasScale

!***************************************************************************
!
! \brief       Scale of a gas's FLUXNET vertical-advection column, by species.
! \author      Jonathan Muller
! \note        Advection is w times the molar density, so it starts from the
!              mmol m-3 basis rather than the mole basis: one further factor of
!              1d3 for every species except water, whose target already is the
!              mmol basis. Reproduces the 1d3 / none / 1d6 three-way switch
!              this replaces, without naming a slot.
!***************************************************************************
real(kind = dbl) function FluxnetGasAdvScale(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: species
    integer :: rec
    real(kind = dbl), external :: FluxnetGasScale

    rec = gas_slot - firstGas + 1
    if (rec >= 1 .and. rec <= min(EddyFlowProj%gas_num, MaxNumGases)) then
        species = EddyFlowProj%gas(rec)%var
        call uppercase(species)
        if (trim(adjustl(species)) == 'H2O') then
            FluxnetGasAdvScale = 1d0
            return
        end if
    end if
    FluxnetGasAdvScale = FluxnetGasScale(gas_slot) * 1d3
end function FluxnetGasAdvScale

!***************************************************************************
!
! \brief       Whether a gas slot holds water.
! \author      Jonathan Muller
! \note        Water is carved out of the per-gas unit handling: it is the one
!              species held internally on the mmol basis, so the umol default
!              of GasFullOutputUnits would mislabel it. Answered from the gas
!              record and not from the slot number, because a second water
!              record sits well past the historical h2o slot.
!***************************************************************************
logical function GasSlotIsWater(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: species
    integer :: rec

    GasSlotIsWater = .false.
    rec = gas_slot - firstGas + 1
    if (rec < 1 .or. rec > min(EddyFlowProj%gas_num, MaxNumGases)) return

    species = EddyFlowProj%gas(rec)%var
    call uppercase(species)
    GasSlotIsWater = trim(adjustl(species)) == 'H2O'
end function GasSlotIsWater

!***************************************************************************
!
! \brief       Whether a gas has a fitted transfer function to correct it with.
! \author      Jonathan Muller
! \note        A gas is corrected in situ only if some class of its own got a
!              cut-off frequency. Water is binned by relative humidity and
!              every other gas by the month group, so the range to search
!              differs - asking for water's fit in the month range demands
!              something no hygrometer ever has.
!
!              Asked by both the writer, which decides whether the assessment
!              file has anything to say, and the diagnostics, which reports it.
!              The two disagreed before: the writer tested whether water's
!              ensemble spectra existed and the report tested whether each gas
!              had been fitted, so the file could be refused while the report
!              called four gases PASS.
!***************************************************************************
logical function GasHasSpectralFit(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    integer :: cls
    logical, external :: GasSlotIsWater

    GasHasSpectralFit = .false.
    if (gas_slot < firstGas .or. gas_slot > lastGas) return

    if (GasSlotIsWater(gas_slot)) then
        do cls = RH10, RH90
            if (RegPar(gas_slot, cls)%fc /= error) GasHasSpectralFit = .true.
        end do
    else
        do cls = 1, MaxGasClasses
            if (RegPar(gas_slot, cls)%fc /= error) GasHasSpectralFit = .true.
        end do
    end if
end function GasHasSpectralFit

!***************************************************************************
!
! \brief       The gas slot holding the site's primary water measurement.
! \author      Jonathan Muller
! \note        The hygrometer whose H, LE and ET carry the bare column names,
!              which is the first water record unless one is flagged. Every
!              hygrometer now gets its own family - see WaterOutSlots - so this
!              no longer decides which humidity the site *has*, only which one
!              a reader finds under the unsuffixed FLUXNET spelling. Gases that
!              need *their own* water for a WPL correction use
!              E2Col(gas)%moist_ref, as they already did.
!
!              Returns 0 when the project describes no water at all. Callers
!              must treat that as "not performed" - the quantity is `error`
!              and its column is not written - rather than computing from a
!              slot that holds something else.
!
!              Record-derived and not presence-derived, for the same reason
!              GasOutputLabel is: output headers are written before the first
!              data file is read, and FCC has no E2Col. Callers apply their
!              own availability test on top.
!
!              This replaces reading the literal `h2o` slot. That slot is
!              record two, which is water only by convention - and with
!              per-instrument water there is more than one water record, so
!              even the convention no longer identifies a unique slot.
!***************************************************************************
integer function PrimaryWaterSlot()
    use m_common_global_var
    implicit none
    integer, external :: DesignatedGasSlot

    !> Same designation the FLUXNET naming uses, so the bare H2O column and
    !> the site's LE and ET come from one hygrometer rather than two.
    PrimaryWaterSlot = DesignatedGasSlot('H2O')
end function PrimaryWaterSlot

!***************************************************************************
!
! \brief       Slot to gate the one-per-site water columns on.
! \author      Jonathan Muller
! \note        The full output and the QC details carry LE_strg, un_LE, LE_scf
!              and the water QC cells only when the site measures humidity.
!              That question used to be asked as OutVarPresent(h2o) - the
!              historical slot, which is water only by convention - so a
!              project declaring its water elsewhere could have those columns
!              decided by whichever species record two happened to hold.
!
!              Header and rows are written by different routines in different
!              executables, and they must agree exactly or the file shifts.
!              So the choice lives here, in one function both call, rather
!              than being spelled out at each of the seventeen gates.
!
!              Falls back to the historical slot when the project describes no
!              water at all, which is what those gates have always evaluated
!              in that case.
!***************************************************************************
!***************************************************************************
!
! \brief       Ordered slot list the statistics files are laid out on.
! \author      Jonathan Muller
! \note        The seven st1..st7 files write five per-slot families - mean,
!              var, st_dev, skw, kur - and their header names the variables
!              once. The two used to agree by coincidence: the writer looped
!              `u, pe`, and while E2NumVar was 14 that produced exactly the
!              twelve names the header lists. E2NumVar is now 102 - 64 gas
!              slots and 32 per-instrument cell slots - so the rows became
!              seven times wider than their own header and the files stopped
!              being readable. No fixture enables them, which is why it went
!              unseen.
!
!              Both sides now walk this list. It is the historical set
!              generalised: the anemometer, then one entry per *configured*
!              gas - by record count rather than by presence, so a gas named
!              without a column keeps its column of error codes exactly as
!              the fourth slot always did - then instrument 1's cell
!              temperature and pressure, then ambient. At four gases it
!              reproduces the historical twelve exactly.
!***************************************************************************
subroutine StatsLayoutSlots(slots, nslots)
    use m_common_global_var
    implicit none
    integer, intent(out) :: slots(E2NumVar)
    integer, intent(out) :: nslots
    integer :: gas

    slots = 0
    slots(1) = u
    slots(2) = v
    slots(3) = w
    slots(4) = ts
    nslots = 4

    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        nslots = nslots + 1
        slots(nslots) = gas
    end do

    !> tc and pi are instrument 1's cell block, which is where these two have
    !> always pointed; te and pe are ambient.
    nslots = nslots + 1; slots(nslots) = tc
    nslots = nslots + 1; slots(nslots) = pi
    nslots = nslots + 1; slots(nslots) = te
    nslots = nslots + 1; slots(nslots) = pe
end subroutine StatsLayoutSlots

!***************************************************************************
!
! \brief       Ordered gas slots the full output is laid out on.
! \author      Jonathan Muller
! \note        Same failure as StatsLayoutSlots, in the other output file.
!              Both row writers looped firstGas..lastGas while the header
!              named a smaller set, so at lastGas = 68 the row carried sixty
!              phantom gas blocks: fifteen fields each, nine hundred fields
!              past the end of a header describing a hundred and ninety four.
!              Everything after the gas block was shifted, and FCC parses the
!              ex record by comma count, so it consumed the row without
!              complaint.
!
!              The whole gas block, which is what both sides walk. The
!              OutVarPresent guard at each of the eight use sites does the
!              narrowing, and it is the *same* guard in the header and in the
!              row, so they cannot disagree.
!***************************************************************************
subroutine FullOutputGasSlots(slots, nslots)
    use m_common_global_var
    implicit none
    integer, intent(out) :: slots(GHGNumVar)
    integer, intent(out) :: nslots
    integer :: gas

    slots = 0
    nslots = 0

    do gas = firstGas, lastGas
        nslots = nslots + 1
        slots(nslots) = gas
    end do
end subroutine FullOutputGasSlots

!***************************************************************************
!
! \brief       The FLUXNET/essentials gas layout: one slot per column family.
! \author      Jonathan Muller
! \note        RP writes this file, FCC rewrites it, and ReadExRecord parses it
!              by field position. Three programs, one order - and until this
!              existed each derived it separately, RP from nFluxnetLayoutSlots
!              and the other two from min(gas_num, MaxNumGases). They agreed
!              only because the layout happened to be the record list; nothing
!              made them agree.
!
!              Position, not arithmetic. Callers must index this list rather
!              than assume position n is slot firstGas + n - 1: the order is
!              free to differ from record order, and a caller computing the
!              slot itself would read one gas's field into another's slot with
!              nothing to flag it.
!
!              CO2, H2O and CH4 lead, in that order, then every other record
!              in the order the project declares them. FLUXNET requires those
!              three variables, so their columns exist whatever the site
!              measures - a species no record names is carried by a slot past
!              the configured records, which holds no data and so reports the
!              error label throughout, exactly as a record with no column does.
!
!              The designation is by NAME and ignores %col: a project that
!              declares CH4 without a column still holds position three with
!              it. That is a different question from DesignatedGasSlot, which
!              requires a column because it answers which instrument supplies
!              FC and LE - a quantity, not a column heading.
!***************************************************************************
subroutine FluxnetLayoutGasSlots(slots, nslots)
    use m_common_global_var
    implicit none
    integer, intent(out) :: slots(GHGNumVar)
    integer, intent(out) :: nslots
    integer :: rec
    integer :: nrec
    integer :: nsynth
    integer :: placed(MaxNumGases)
    integer :: req
    integer :: designated
    character(32) :: species
    character(3) :: required(3)

    required = (/ 'CO2', 'H2O', 'CH4' /)

    slots = 0
    nslots = 0
    placed = 0
    nsynth = 0
    nrec = min(EddyFlowProj%gas_num, MaxNumGases)

    !> The three required variables first, one entry each.
    do req = 1, 3
        designated = 0
        do rec = 1, nrec
            species = EddyFlowProj%gas(rec)%var
            call uppercase(species)
            if (trim(adjustl(species)) /= trim(required(req))) cycle
            !> The first record of the species holds it until one the project
            !> flags takes it, so the choice matches the naming and FC/LE.
            if (designated == 0) designated = rec
            if (EddyFlowProj%gas(rec)%fluxnet_default == 1) then
                designated = rec
                exit
            end if
        end do

        if (designated > 0) then
            if (firstGas + designated - 1 > lastGas) cycle
            nslots = nslots + 1
            slots(nslots) = firstGas + designated - 1
            placed(designated) = 1
        else
            !> Nothing names it. Carry it on a slot past the configured
            !> records, which no record owns and nothing fills, so every field
            !> of it is written as the error label. Dropped rather than
            !> overrun when the project already fills the slot space: the
            !> arrays this indexes are exactly MaxNumGases wide.
            if (firstGas + nrec + nsynth > lastGas) cycle
            nslots = nslots + 1
            slots(nslots) = firstGas + nrec + nsynth
            nsynth = nsynth + 1
        end if
    end do

    !> Then every record the required block did not already take.
    do rec = 1, nrec
        if (placed(rec) == 1) cycle
        if (firstGas + rec - 1 > lastGas) exit
        nslots = nslots + 1
        slots(nslots) = firstGas + rec - 1
    end do
end subroutine FluxnetLayoutGasSlots

!***************************************************************************
!
! \brief       Whether a gas slot has an in-situ spectral correction.
! \author      Jonathan Muller
! \note        The in-situ (Fratini/Ibrom) corrections are derived from the
!              measured cospectra of CO2, water and CH4; there is no such
!              derivation for an arbitrary trace gas, which is why those go
!              through FCC with BPCF = 1 instead.
!
!              This was asked as OutVarPresent(gas4) - "is the fourth slot
!              occupied" - which both missed a COS on record five when record
!              four was empty, and mistook a second CO2 sitting in slot eight
!              for a gas that needs the FCC-only path.
!***************************************************************************
!***************************************************************************
!
! \brief       First gas slot the project configures, or the historical fifth.
! \author      Jonathan Muller
! \note        For the handful of places that need "an" analyser rather than a
!              particular one - the logger software version, say. They read
!              E2Col(co2), which is the first record only when CO2 happens to
!              be it.
!***************************************************************************
integer function FirstConfiguredGasSlot()
    use m_common_global_var
    implicit none
    integer :: gas

    !> Unreachable as a result: the loop below assigns before returning for
    !> any project with a configured gas, and a project with none is refused
    !> by ApplyGasRecords. It is the first slot, said as the first slot.
    FirstConfiguredGasSlot = firstGas
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (EddyFlowProj%gas(gas - firstGas + 1)%col <= 0) cycle
        FirstConfiguredGasSlot = gas
        return
    end do
end function FirstConfiguredGasSlot

!***************************************************************************
!
! \brief       First gas slot whose analyser model contains `fragment`, or 0.
! \author      Jonathan Muller
! \note        For the gates that name a *particular instrument* rather than a
!              particular gas: the LI-7700 multiplier block, and the AGC/RSSI
!              columns whose own heading says LI-7200 or LI-7500. Those asked
!              E2Col(ch4) and E2Col(co2) - slots seven and five - which is that
!              analyser only when the project happens to order its records so.
!              A site carrying both a 7200 and a 7500 had each one's column
!              decided by whichever firmware version slot five reported, and a
!              site whose first record sits on a third instrument had both
!              decided by something unrelated to either.
!
!              Distinct from FirstConfiguredGasSlot, which answers "any
!              analyser" for things like the logger version. This answers "that
!              analyser", and returns 0 when the site carries none - the arm the
!              slot lookup took when its instrument was unset.
!
!              Matching is by substring, as InterpretLicorDiagnostics has always
!              done it, because the model string carries a trailing revision.
!***************************************************************************
integer function GasSlotByInstrModel(fragment)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: fragment
    integer :: gas

    GasSlotByInstrModel = 0
    do gas = firstGas, lastGas
        if (index(E2Col(gas)%instr%model, fragment) == 0) cycle
        GasSlotByInstrModel = gas
        return
    end do
end function GasSlotByInstrModel

!***************************************************************************
!
! \brief       Firmware reported by the first analyser whose model matches.
! \author      Jonathan Muller
! \note        The companion to GasSlotByInstrModel, holding the one rule its
!              three callers would otherwise each have to restate: a site with
!              no such analyser reports version zero, which compares older than
!              any threshold. That is the arm an unpopulated
!              E2Col(co2)%instr%sw_ver took, so the absent-instrument case
!              keeps behaving as it did.
!
!              Two thresholds are asked of this, and they are different
!              questions: 5.3.0 changed how the diagnostic word *encodes* the
!              signal strength, and 6.0.0 changed what the column is *called*.
!***************************************************************************
function InstrSwVerFor(fragment) result(ver)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: fragment
    type(SwVerType) :: ver
    integer :: slot
    integer, external :: GasSlotByInstrModel

    ver = SwVerType(0, 0, 0)
    slot = GasSlotByInstrModel(fragment)
    if (slot > 0) ver = E2Col(slot)%instr%sw_ver
end function InstrSwVerFor

logical function HasInSituSpectralCorrection(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: species
    character(32), external :: GasOutputLabel

    HasInSituSpectralCorrection = .false.
    if (gas_slot < firstGas .or. gas_slot > lastGas) return

    species = GasOutputLabel(gas_slot)
    call lowercase(species)
    select case (trim(adjustl(species)))
        case ('co2', 'h2o', 'ch4'); HasInSituSpectralCorrection = .true.
    end select
end function HasInSituSpectralCorrection

!***************************************************************************
!
! \brief       Width and legend of the packed statistical-flag strings.
! \author      Jonathan Muller
! \note        Eight of the Vickers and Mahrt flags carry one digit per
!              variable, behind a leading filler digit. The full output's
!              units row spells out which variable each digit stands for:
!              `8u/v/w/ts/co2/h2o/ch4/<gas4>`, eight names.
!
!              The row emitted `'8' // CharHF%sr(2:FlagStrLen)`. FlagStrLen
!              is 1 + GHGNumVar = 69, so the cell carried sixty-nine
!              characters against a legend naming eight - the real flags
!              followed by sixty '9's for gas slots the project does not
!              configure. Readable by eye, unreadable by column.
!
!              Both sides now ask this. nvars is the anemometric block plus
!              one per *configured* gas, so at four gases the cell is nine
!              characters again, exactly as it was before the slot capacity
!              grew, and the legend names exactly what the digits are.
!
!              The legend is built from the same tags the column names use,
!              so a project measuring CO2 twice reads `co2` and `co2_2`
!              rather than two identical entries.
!***************************************************************************
subroutine StatisticalFlagVars(nvars, legend)
    use m_common_global_var
    implicit none
    integer, intent(out) :: nvars
    character(*), intent(out) :: legend
    character(64) :: tags(GHGNumVar)
    integer :: gas

    call SpectralVarTags(tags)

    legend = 'u/v/w/ts'
    nvars = ts - u + 1
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        nvars = nvars + 1
        legend = trim(legend) // '/' // trim(tags(gas))
    end do
end subroutine StatisticalFlagVars

!***************************************************************************
!
! \brief       Legend of the gas-only flag strings: time-lag hard and soft.
! \author      Jonathan Muller
! \note        These are still four digits wide. TestTimeLag packs its flags
!              into a base-10 integer of exactly four gas digits and says so
!              in a comment - generalising it needs that packing replaced by
!              PackFlagString, as the other six tests already have been.
!              Until then the legend names the first four gas records, from
!              their own tags rather than from the literals co2/h2o/ch4, so
!              at least it does not misname what is there.
!***************************************************************************
subroutine TimelagFlagLegend(nvars, legend)
    use m_common_global_var
    implicit none
    integer, intent(out) :: nvars
    character(*), intent(out) :: legend
    character(64) :: tags(GHGNumVar)
    integer :: gas

    call SpectralVarTags(tags)

    !> Every configured gas. This ran co2..gas4 to match a four-digit packing
    !> in TestTimeLag; that packing is a per-variable string now, so the legend
    !> follows the gases the test actually ran over.
    legend = ''
    nvars = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        nvars = nvars + 1
        if (len_trim(tags(gas)) > 0) then
            if (nvars == 1) then
                legend = trim(tags(gas))
            else
                legend = trim(legend) // '/' // trim(tags(gas))
            end if
        end if
    end do
end subroutine TimelagFlagLegend

!***************************************************************************
!
! \brief       The per-gas field suffixes of the dynamic metadata file.
! \author      Jonathan Muller
! \note        Fourteen fields per analyser, in the order the historical index
!              table lays them out (co2_irga_manufacturer = 20 through
!              co2_irga_tau = 33, then the same for h2o, ch4 and gas4).
!
!              measure_type is the odd one: it has no `irga_` in its name,
!              which is why the suffixes are listed rather than composed.
!
!              Named here so the header matcher and the reader cannot
!              disagree about how many fields there are or what they are
!              called - four unrolled blocks of fifteen statements each had
!              already drifted once, with the fourth gas's measure_type
!              landing in a separation field.
!***************************************************************************
subroutine DynMDGasFieldNames(names, nfields)
    use m_common_global_var
    implicit none
    character(*), intent(out) :: names(nDynMDGasFields)
    integer, intent(out) :: nfields

    names(1)  = '_irga_manufacturer'
    names(2)  = '_irga_model'
    names(3)  = '_measure_type'
    names(4)  = '_irga_northward_separation'
    names(5)  = '_irga_eastward_separation'
    names(6)  = '_irga_vertical_separation'
    names(7)  = '_irga_tube_length'
    names(8)  = '_irga_tube_diameter'
    names(9)  = '_irga_tube_flowrate'
    names(10) = '_irga_kw'
    names(11) = '_irga_ko'
    names(12) = '_irga_hpath_length'
    names(13) = '_irga_vpath_length'
    names(14) = '_irga_tau'
    nfields = nDynMDGasFields
end subroutine DynMDGasFieldNames

!***************************************************************************
!
! \brief       Gas slot a `<stem><suffix>` column name refers to, or 0.
! \author      Jonathan Muller
! \note        The drift subsystem reads two kinds of per-gas column whose
!              names are authored by the user: `<gas>_offset` and `<gas>_ref`
!              in the dynamic metadata file, and `<gas>_ref` again as a raw
!              data column holding reference counts.
!
!              Both sides used to spell the four names out - co2_ref = 9,
!              h2o_ref = 10, ch4_ref = 11, gas4_ref = 12, packed into the same
!              array as the offsets at 5..8 and read back through a `- 4`.
!              A fifth gas had no name it could be given.
!
!              Resolution order, and both halves matter:
!
!              1. The record's own label, as GasOutputLabel spells it,
!                 including the `_2` a second record of the same species
!                 gets. That is what lets a project name `n2o_ref` or
!                 `co2_2_ref` and be understood.
!              2. The legacy aliases co2/h2o/ch4/gas4, through
!                 HistoricGasSlot. Dynamic metadata files in the wild spell
!                 those, and renaming them would break every one - so they
!                 are accepted as *aliases for the first four slots*, exactly
!                 as LegacySpectralVarTag widens the (co)spectra reader.
!
!              One helper, because the dynamic metadata file and the raw
!              column scan must agree about what `n2o_ref` means. Spelled out
!              twice, they would not have to.
!***************************************************************************
integer function GasSlotFromDynMDTag(field, suffix)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: field
    character(*), intent(in) :: suffix
    character(64) :: stem
    character(64) :: tags(GHGNumVar)
    integer :: gas
    integer :: cut
    integer, external :: HistoricGasSlot
    logical, external :: IsHistoricGasVar

    GasSlotFromDynMDTag = 0

    cut = len_trim(field) - len_trim(suffix)
    if (cut <= 0) return
    if (field(cut + 1 : len_trim(field)) /= trim(suffix)) return
    stem = field(1:cut)
    call lowercase(stem)

    !> The record's own name first, so a project that measures two of a
    !> species can address the second one.
    call SpectralVarTags(tags)
    do gas = firstGas, lastGas
        if (len_trim(tags(gas)) == 0) cycle
        if (trim(tags(gas)) /= trim(stem)) cycle
        GasSlotFromDynMDTag = gas
        return
    end do

    !> Then the four legacy names, for files written before records existed.
    !> 'gas4' is not a species, so IsHistoricGasVar does not know it.
    if (trim(stem) == 'gas4') then
        GasSlotFromDynMDTag = histGas4
        return
    end if
    if (IsHistoricGasVar(stem)) GasSlotFromDynMDTag = HistoricGasSlot(stem)
end function GasSlotFromDynMDTag

!***************************************************************************
!
! \brief       The gas slot holding the site's primary carbon dioxide.
! \author      Jonathan Muller
! \note        The companion to PrimaryWaterSlot, and it exists for one
!              caller: conditional eddy covariance is defined on a CO2/water
!              *pair* (Zahn et al. 2022), so unlike everything else in this
!              file it does not generalise to N gases - it needs to know
!              which two. What it must not do is take that pair from slots
!              five and six, which are CO2 and water by convention only.
!
!              Returns 0 when the project describes no CO2. Callers treat
!              that as "not performed", the same way they treat no water.
!***************************************************************************
integer function PrimaryCarbonSlot()
    use m_common_global_var
    implicit none
    integer, external :: DesignatedGasSlot

    PrimaryCarbonSlot = DesignatedGasSlot('CO2')
end function PrimaryCarbonSlot

!***************************************************************************
!
! \brief       The record a project designates for a species, as a gas slot.
! \author      Jonathan Muller
! \note        A species can be measured more than once. The record whose
!              gas_<i>_fluxnet_default is set is the site's default for it:
!              the one whose FLUXNET columns carry the bare species name, and
!              the one the one-per-site quantities are computed from. Those
!              two must agree - were the naming to follow the flag and the
!              flux not, FC and the bare CO2 column would describe different
!              instruments, which is worse than either choice alone.
!
!              Unflagged, the first record of the species is designated. That
!              is what every project written before this flag existed already
!              gets, so adding the key changes no existing result.
!
!              Only records naming a column are eligible: a species the site
!              declares but does not measure has nothing to designate. Returns
!              0 when the project describes no such species, which callers
!              treat as "not performed".
!***************************************************************************
integer function DesignatedGasSlot(wantedVar)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: wantedVar
    integer :: gas
    integer :: rec
    character(32) :: species
    character(32) :: wanted

    wanted = wantedVar
    call uppercase(wanted)

    DesignatedGasSlot = 0
    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (EddyFlowProj%gas(rec)%col <= 0) cycle
        species = EddyFlowProj%gas(rec)%var
        call uppercase(species)
        if (trim(adjustl(species)) /= trim(adjustl(wanted))) cycle
        !> The first match holds the slot until a flagged one claims it.
        if (DesignatedGasSlot == 0) DesignatedGasSlot = gas
        if (EddyFlowProj%gas(rec)%fluxnet_default == 1) then
            DesignatedGasSlot = gas
            return
        end if
    end do
end function DesignatedGasSlot

!***************************************************************************
!
! \brief       Slot to gate the one-per-site carbon columns on.
! \author      Jonathan Muller
! \note        As PrimaryWaterOutSlot is to PrimaryWaterSlot: falls back to
!              the historical slot when the project describes no CO2, which
!              is what those gates evaluated in that case.
!***************************************************************************
integer function PrimaryCarbonOutSlot()
    use m_common_global_var
    implicit none
    integer, external :: PrimaryCarbonSlot

    PrimaryCarbonOutSlot = PrimaryCarbonSlot()
    if (PrimaryCarbonOutSlot < firstGas) PrimaryCarbonOutSlot = histGas1
end function PrimaryCarbonOutSlot

integer function PrimaryWaterOutSlot()
    use m_common_global_var
    implicit none
    integer, external :: PrimaryWaterSlot

    PrimaryWaterOutSlot = PrimaryWaterSlot()
    if (PrimaryWaterOutSlot < firstGas) PrimaryWaterOutSlot = histGas2
end function PrimaryWaterOutSlot

!***************************************************************************
!
! \brief       Every hygrometer the project describes, with its column suffix.
! \author      Jonathan Muller
! \note        H, LE and ET used to be one column each, computed from the
!              designated hygrometer, because a site was taken to have one
!              humidity however many hygrometers it carried. A site that
!              deliberately fields two got one answer and no way to see the
!              other: the second hygrometer produced its own water flux and
!              nothing else - no latent heat, no evapotranspiration, no
!              sensible heat, no stability of its own.
!
!              So every water record gets a family. The designated one keeps
!              the bare names - H, LE, ET - because those are the FLUXNET
!              spellings and a reader that knows nothing of this change still
!              finds what it expects, reading the same hygrometer it read
!              before. The rest are numbered.
!
!              The suffix is the water's *position in this list*, not a running
!              count of the numbered ones. With the designated hygrometer
!              second of three the families are H_1, H, H_3 - a gap where the
!              bare name sits, rather than H_2 meaning the first record and
!              H_3 the third. The number always names which hygrometer, which
!              is the only property worth having when a project is edited.
!
!              Slots AND suffixes together, from one call, because header and
!              row are written by different routines in different executables:
!              getting the list right and the naming wrong shifts the file
!              exactly as getting the list wrong does. There is nothing here to
!              recompute, so no caller has the chance to disagree.
!
!              Record-derived like PrimaryWaterSlot, for the same reason -
!              headers are written before the first data file is read, and FCC
!              has no E2Col. Only records naming a column are eligible.
!***************************************************************************
subroutine WaterOutSlots(slots, tags, nslots)
    use m_common_global_var
    implicit none
    integer, intent(out) :: slots(GHGNumVar)
    character(*), intent(out) :: tags(GHGNumVar)
    integer, intent(out) :: nslots
    integer :: gas
    integer :: rec
    integer :: i
    integer :: designated
    integer, external :: PrimaryWaterSlot
    logical, external :: GasSlotIsWater

    slots = 0
    tags = ''
    nslots = 0
    designated = PrimaryWaterSlot()

    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (EddyFlowProj%gas(rec)%col <= 0) cycle
        if (.not. GasSlotIsWater(gas)) cycle
        nslots = nslots + 1
        slots(nslots) = gas
    end do

    !> Suffixes only once the list is complete: whether the designated
    !> hygrometer is in it at all decides whether any name is left bare, and
    !> the first record cannot know that.
    do i = 1, nslots
        if (slots(i) == designated) then
            tags(i) = ''
        else
            write(tags(i), '(a,i0)') '_', i
        end if
    end do
end subroutine WaterOutSlots

!***************************************************************************
!
! \brief       Full-output column-name stems, one per configured gas slot.
! \author      Jonathan Muller
! \note        Returns 'co2_', 'h2o_', … indexed by gas slot, lowercased and
!              suffixed for repeats, so a project measuring CO2 on two
!              analysers gets 'co2_' and 'co2_2_' rather than two identical
!              column families. The FLUXNET row solves the same problem the
!              same way in SelectFluxnetGasSlots; without it here, widening
!              the full output would trade missing columns for duplicate ones.
!
!              Empty for a slot with no configured gas, which is what the
!              OutVarPresent guards at every use site already test for.
!***************************************************************************
subroutine FullOutputGasTags(tags)
    use m_common_global_var
    implicit none
    character(*), intent(out) :: tags(GHGNumVar)
    character(32) :: label
    character(32) :: seen(MaxNumGases)
    integer :: nseen
    integer :: gas, j, k, repeat, occurrence
    character(32), external :: GasOutputLabel

    !> Two passes. A species measured once keeps its bare name - cos_, n2o_ -
    !> and one measured more than once has *every* occurrence numbered:
    !> h2o_1_ and h2o_2_, not h2o_ and h2o_2_. Whether a name needs a number
    !> depends on the total count, which the first occurrence cannot know, so
    !> the labels are collected before any of them is written.
    !>
    !> Leaving the first occurrence bare made the pair asymmetric: h2o_flux
    !> read as the site's water flux when it was one of two, and a reader
    !> keying on it silently got whichever happened to be recorded first.
    tags = ''
    nseen = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        label = GasOutputLabel(gas)
        call lowercase(label)
        nseen = nseen + 1
        seen(nseen) = label
    end do

    do k = 1, nseen
        gas = firstGas + k - 1
        if (len_trim(seen(k)) == 0) cycle

        !> How many records name this species, and which of them this is.
        repeat = 0
        occurrence = 0
        do j = 1, nseen
            if (trim(seen(j)) /= trim(seen(k))) cycle
            repeat = repeat + 1
            if (j <= k) occurrence = occurrence + 1
        end do

        if (repeat == 1) then
            tags(gas) = trim(seen(k)) // '_'
        else
            write(tags(gas), '(a,i0,a)') trim(seen(k)) // '_', occurrence, '_'
        end if
    end do
end subroutine FullOutputGasTags

!***************************************************************************
!
! \brief       Variable names used in the (co)spectra files, per slot.
! \author      Jonathan Muller
! \note        The spectra and full-cospectra files name their columns
!              `var(<tag>)` and `cov(w_<tag>)`, and the writer and the reader
!              that imports them back for the in-situ corrections must agree
!              exactly - a name that does not match is simply not imported, and
!              the gas silently gets no correction factor. One helper, so they
!              cannot drift.
!
!              The four anemometer channels keep their literal names: those
!              are fixed measurements, not gases. **Every gas slot is named
!              from its record**, including the historical four. Pinning them
!              to co2/h2o/ch4/gas4 named a position rather than a species, so
!              on a project that orders its records differently the label lied
!              about which gas a column held - and two slots could end up
!              carrying the same name, which breaks the round trip through
!              GasIndexFromLabel.
!
!              Readers of files written before this keep working through an
!              alias set, since the old names are still what those files say.
!***************************************************************************
subroutine SpectralVarTags(tags)
    use m_common_global_var
    implicit none
    character(*), intent(out) :: tags(GHGNumVar)
    character(64) :: gas_tags(GHGNumVar)
    integer :: gas

    tags = ''
    tags(u)    = 'u'
    tags(v)    = 'v'
    tags(w)    = 'w'
    tags(ts)   = 'ts'

    call FullOutputGasTags(gas_tags)
    do gas = firstGas, lastGas
        !> FullOutputGasTags returns a stem with a trailing underscore, since
        !> the full output concatenates it directly onto a quantity name.
        if (len_trim(gas_tags(gas)) > 1) &
            tags(gas) = gas_tags(gas)(1:len_trim(gas_tags(gas)) - 1)
    end do
end subroutine SpectralVarTags

!***************************************************************************
!
! \brief       The legacy on-disk name of a gas slot, or blank if it has none.
! \author      Jonathan Muller
! \note        Files written before SpectralVarTags became record-derived name
!              the first four gas slots co2/h2o/ch4/gas4 regardless of what
!              those slots held. A reader matching columns by name has to
!              accept both spellings or it silently imports nothing for those
!              gases - and "no cospectrum" means "no correction applied",
!              which is a quiet loss rather than a loud one.
!
!              Only ever used to *widen* what a reader accepts. Nothing
!              writes these.
!***************************************************************************
character(32) function LegacySpectralVarTag(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot

    select case (gas_slot)
        case (histGas1);  LegacySpectralVarTag = 'co2'
        case (histGas2);  LegacySpectralVarTag = 'h2o'
        case (histGas3);  LegacySpectralVarTag = 'ch4'
        case (histGas4); LegacySpectralVarTag = 'gas4'
        case default; LegacySpectralVarTag = ''
    end select
end function LegacySpectralVarTag

!***************************************************************************
!
!> \brief       Human-readable species names for the spectral reports
!> \author      EddyFlow
!> \note        Every gas slot is named by its own record. No slot is assumed
!>              to hold a particular species: slots are assigned by record
!>              order (slot = firstGas + i - 1), so the constants co2, h2o
!>              and ch4 are aliases for records one to three and say nothing
!>              about what those records declare. This used to pin the first
!>              three to 'co2'/'h2o'/'ch4', which named the wrong gas on any
!>              project that ordered its records differently.
!>
!>              Free text throughout - the assessment reader skips these
!>              header lines rather than matching them.
!
!***************************************************************************
subroutine SpectralGasNames(names)
    use m_common_global_var
    implicit none
    character(*), intent(out) :: names(GHGNumVar)
    character(64) :: gas_tags(GHGNumVar)
    integer :: gas, rec

    names = ''
    call FullOutputGasTags(gas_tags)
    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        !> FullOutputGasTags returns a stem with a trailing underscore,
        !> since the full output concatenates it onto a quantity name. A gas
        !> configured without a column is named too - GasOutputLabel falls
        !> back to the species its record declares.
        if (len_trim(gas_tags(gas)) > 1) &
            names(gas) = gas_tags(gas)(1:len_trim(gas_tags(gas)) - 1)
    end do
end subroutine SpectralGasNames

!***************************************************************************
!
!> \brief       The species and analyser a spectral assessment block belongs to
!> \author      Jonathan Muller
!> \note        The block name alone cannot say which analyser a block is for.
!>              It is an ordinal over repeats of a label - CO2_1 is whichever
!>              CO2 record came first - so the moment a project lists its
!>              records in a different order, the same file hands CO2_1's
!>              transfer function to the other analyser. On CH-LAE that
!>              happened on an ordinary re-save: the first CO2 record went from
!>              the LI-7200 to the MIRO, and the two cells are 943 hPa and
!>              70 hPa apart.
!>
!>              So the block states what it is for. Written past `fc`, where
!>              every reader of this format stops - the block name is sliced at
!>              the word TFP - so a file carrying it still parses in a build
!>              that knows nothing about it, exactly as `groups=` does.
!
!***************************************************************************
subroutine SpectralBlockStamp(gas_slot, stamp)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(*), intent(out) :: stamp
    integer :: rec

    stamp = ''
    rec = gas_slot - firstGas + 1
    if (rec < 1 .or. rec > min(EddyFlowProj%gas_num, MaxNumGases)) return
    stamp = '   var=' // trim(adjustl(EddyFlowProj%gas(rec)%var)) // &
            ' instr=' // trim(adjustl(EddyFlowProj%gas(rec)%instr))
end subroutine SpectralBlockStamp

!***************************************************************************
!
!> \brief       The gas slot a stamped assessment block names, or 0
!> \author      Jonathan Muller
!> \note        Answers only when the line carries both tokens and exactly one
!>              record matches them. Anything else - an unstamped file, a
!>              species this project does not measure, or two records that
!>              genuinely declare the same species on the same analyser - is
!>              left to the caller's name match, which is what every file
!>              written before the stamp relies on.
!
!***************************************************************************
integer function SlotFromSpectralStamp(dataline)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: dataline
    character(64) :: want_var, want_instr
    integer :: gas, rec, hits

    SlotFromSpectralStamp = 0
    call SpectralStampToken(dataline, 'var=', want_var)
    call SpectralStampToken(dataline, 'instr=', want_instr)
    if (len_trim(want_var) == 0 .or. len_trim(want_instr) == 0) return

    hits = 0
    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (trim(adjustl(EddyFlowProj%gas(rec)%var)) /= trim(want_var)) cycle
        if (trim(adjustl(EddyFlowProj%gas(rec)%instr)) /= trim(want_instr)) cycle
        hits = hits + 1
        SlotFromSpectralStamp = gas
    end do
    if (hits /= 1) SlotFromSpectralStamp = 0
end function SlotFromSpectralStamp

!***************************************************************************
!
!> \brief       The value of one `key=value` token on a block header line
!> \author      Jonathan Muller
!> \note        Case-insensitive on the value, since the writer uppercases the
!>              block name and a hand-edited file may uppercase anything else.
!
!***************************************************************************
subroutine SpectralStampToken(dataline, key, value)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: dataline
    character(*), intent(in) :: key
    character(*), intent(out) :: value
    integer :: at, stop_at

    value = ''
    at = index(dataline, key)
    if (at == 0) return
    at = at + len_trim(key)
    stop_at = at
    do while (stop_at <= len_trim(dataline))
        if (dataline(stop_at:stop_at) == ' ') exit
        stop_at = stop_at + 1
    end do
    if (stop_at <= at) return
    value = dataline(at:stop_at - 1)
    call lowercase(value)
end subroutine SpectralStampToken

!***************************************************************************
!
! \brief       Full-output scales and labels for every configured gas slot.
! \author      Jonathan Muller
! \note        One call fills the header and the row writer alike, so the
!              units a column is labelled with and the factor its value is
!              scaled by cannot drift apart - the failure this replaces was
!              exactly that, with a per-slot label and no per-slot scale.
!
!              Water takes fixed mmol labels and unit scales rather than the
!              unit lookup, matching the hard-coded water arms this replaces.
!***************************************************************************
subroutine GasFullOutputUnitsAll(flux_sc, dens_sc, &
    flux_label, conc_label, mixr_label, dens_label)
    use m_common_global_var
    implicit none
    real(kind = dbl), intent(out) :: flux_sc(GHGNumVar), dens_sc(GHGNumVar)
    character(*), intent(out) :: flux_label(GHGNumVar), conc_label(GHGNumVar)
    character(*), intent(out) :: mixr_label(GHGNumVar), dens_label(GHGNumVar)
    integer :: gas
    logical, external :: GasSlotIsWater
    character(32), external :: GasUnitIn

    flux_sc = 1d0
    dens_sc = 1d0
    flux_label = ''
    conc_label = ''
    mixr_label = ''
    dens_label = ''

    do gas = firstGas, lastGas
        if (GasSlotIsWater(gas)) then
            flux_label(gas) = '[mmol+1s-1m-2]'
            conc_label(gas) = '[mmol+1mol_a-1]'
            mixr_label(gas) = '[mmol+1mol_d-1]'
            dens_label(gas) = '[mmol+1m-3]'
        else
            call GasFullOutputUnits(GasUnitIn(gas), flux_sc(gas), dens_sc(gas), &
                flux_label(gas), conc_label(gas), mixr_label(gas), dens_label(gas))
        end if
    end do
end subroutine GasFullOutputUnitsAll

!***************************************************************************
!
! \brief       Full-output scales and labels for one gas, from its input unit.
! \author      Jonathan Muller
! \note        The full output carries a units row, so unlike the FLUXNET row
!              it reports each gas in the basis the project configured. The
!              body was always generic; only its name and its single caller
!              tied it to the fourth slot, which is why gases 5+ were absent
!              from that file rather than merely mislabelled in it.
!
!              Water is not passed through here. Its internal basis is mmol
!              rather than umol, so the default arm's umol labels would be
!              wrong for it; the callers keep water's own fixed labels, the
!              same carve-out water has everywhere else in this work.
!***************************************************************************
subroutine GasFullOutputUnits(unit_in, flux_scale, dens_scale, &
    flux_label, conc_label, mixr_label, dens_label)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: unit_in
    real(kind = dbl), intent(out) :: flux_scale
    real(kind = dbl), intent(out) :: dens_scale
    character(*), intent(out) :: flux_label
    character(*), intent(out) :: conc_label
    character(*), intent(out) :: mixr_label
    character(*), intent(out) :: dens_label

    select case (trim(adjustl(unit_in)))
        case ('ppb', 'nmol_mol', 'nmol/mol')
            flux_scale = 1d3
            dens_scale = 1d6
            flux_label = '[nmol+1s-1m-2]'
            conc_label = '[nmol+1mol_a-1]'
            mixr_label = '[nmol+1mol_d-1]'
            dens_label = '[nmol+1m-3]'
        case ('pmol_mol', 'pmol/mol')
            flux_scale = 1d6
            dens_scale = 1d9
            flux_label = '[pmol+1s-1m-2]'
            conc_label = '[pmol+1mol_a-1]'
            mixr_label = '[pmol+1mol_d-1]'
            dens_label = '[pmol+1m-3]'
        case default
            flux_scale = 1d0
            dens_scale = 1d0
            flux_label = '[' // char(194) // char(181) // 'mol+1s-1m-2]'
            conc_label = '[' // char(194) // char(181) // 'mol+1mol_a-1]'
            mixr_label = '[' // char(194) // char(181) // 'mol+1mol_d-1]'
            dens_label = '[mmol+1m-3]'
    end select
end subroutine GasFullOutputUnits

!***************************************************************************
!
! \brief       Every Conditional Eddy Covariance pairing, resolved, with the
!              suffix its columns carry.
! \author      Jonathan Muller
! \note        A pairing is one CO2 channel, one water channel, and any
!              further species partitioned in the octants those two define.
!              Both channels are needed whichever fluxes the pairing is asked
!              to partition, because the octants are built from the signs of
!              w', q' and c' together - "carbon only" still needs the water.
!
!              Stated pairings win. An absent cec_num means the project has
!              never said, and the layout answers instead: one pairing per CO2
!              channel, each with the water on its own analyser. That is the
!              pairing a site means when it says nothing, and it is the one
!              choice that does not silently mix two analysers' time lags and
!              spectral responses into the series that define the octants.
!
!              Columns are suffixed with the analyser the CO2 sits on rather
!              than a bare number, so a reader can tell which instrument a
!              partition describes without consulting the project file. Two
!              pairings on one analyser - a site with two water channels on it
!              - are told apart by an occurrence number after that.
!***************************************************************************
subroutine CecPairs(pairs, npairs)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(out) :: pairs(MaxNumCecPairs)
    integer, intent(out) :: npairs

    integer :: i
    integer :: j
    integer :: k
    integer :: rec
    integer :: slot
    integer :: gas
    integer :: occurrence
    integer :: repeat
    character(48) :: base
    integer, external :: PrimaryCarbonSlot
    integer, external :: CecWaterOnAnalyserOf
    character(32), external :: GasOutputLabel
    logical, save :: crossPairWarned(MaxNumCecPairs) = .false.
    logical, save :: cecWplWarned = .false.

    do i = 1, MaxNumCecPairs
        pairs(i)%meth = 0
        pairs(i)%carbon_slot = 0
        pairs(i)%water_slot = 0
        pairs(i)%n_extra = 0
        pairs(i)%extra_slot = 0
        pairs(i)%tag = ''
    end do
    npairs = 0
    if (EddyFlowProj%do_cec == 0) return

    if (EddyFlowProj%cec_num > 0) then
        do i = 1, min(EddyFlowProj%cec_num, MaxNumCecPairs)
            if (EddyFlowProj%cec_pair(i)%meth == 0) cycle

            slot = CecSlotOfRecord(EddyFlowProj%cec_pair(i)%carbon, 'CO2')
            if (slot == 0) slot = PrimaryCarbonSlot()
            if (slot < firstGas) cycle

            npairs = npairs + 1
            pairs(npairs)%meth = EddyFlowProj%cec_pair(i)%meth
            pairs(npairs)%carbon_slot = slot
            pairs(npairs)%water_slot = &
                CecSlotOfRecord(EddyFlowProj%cec_pair(i)%water, 'H2O')
            if (pairs(npairs)%water_slot == 0) &
                pairs(npairs)%water_slot = CecWaterOnAnalyserOf(slot)

            do k = 1, MaxNumCecExtra
                rec = EddyFlowProj%cec_pair(i)%extra(k)
                if (rec <= 0) cycle
                if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) cycle
                if (EddyFlowProj%gas(rec)%col <= 0) cycle
                gas = firstGas + rec - 1
                !> A pairing's own two channels are already targets one and
                !> two. Naming one again would partition it twice under two
                !> names, so it is dropped rather than duplicated.
                if (gas == pairs(npairs)%carbon_slot) cycle
                if (gas == pairs(npairs)%water_slot) cycle
                pairs(npairs)%n_extra = pairs(npairs)%n_extra + 1
                pairs(npairs)%extra_slot(pairs(npairs)%n_extra) = gas
            end do
        end do
    else
        call AutoCecPairs(pairs, npairs)
    end if

    !> Drop a pairing that cannot be run, and one that repeats an earlier one.
    j = 0
    do i = 1, npairs
        if (pairs(i)%carbon_slot < firstGas .or. pairs(i)%water_slot < firstGas) cycle
        if (CecPairIsRepeated(pairs, j, pairs(i)%carbon_slot, pairs(i)%water_slot)) cycle
        j = j + 1
        if (j /= i) pairs(j) = pairs(i)
    end do
    do i = j + 1, npairs
        pairs(i)%meth = 0
        pairs(i)%carbon_slot = 0
        pairs(i)%water_slot = 0
        pairs(i)%n_extra = 0
        pairs(i)%extra_slot = 0
        pairs(i)%tag = ''
    end do
    npairs = j

    !> Once per run, not once per period: this is called by the header
    !> builders and by both row writers, so an unguarded message is one per
    !> half-hour per pairing and the log the warning exists to be read in is
    !> buried under it. Same guard, and the same reason, as Warning(106).
    do i = 1, npairs
        if (pairs(i)%carbon_slot < firstGas .or. pairs(i)%water_slot < firstGas) cycle
        if (CecSameAnalyser(pairs(i)%carbon_slot, pairs(i)%water_slot)) cycle
        if (crossPairWarned(i)) cycle
        crossPairWarned(i) = .true.
        write(*, '(a,i0,a)') '  Warning(113)> CEC pairing ', i, ': ' &
            // trim(GasOutputLabel(pairs(i)%carbon_slot)) // ' with the water on ' &
            // trim(GasOutputLabel(pairs(i)%water_slot)) // '.'
        write(ulog, '(a,i0,a)') '  Warning(113)> CEC pairing ', i, ': ' &
            // trim(GasOutputLabel(pairs(i)%carbon_slot)) // ' with the water on ' &
            // trim(GasOutputLabel(pairs(i)%water_slot)) // '.'
        call ExceptionHandler(113)
    end do

    !> Same guard and the same reason as the pairing warning above. Asked once
    !> the list is known to be non-empty, so a project with CEC on but nothing
    !> to partition does not also complain about the density correction.
    if (npairs > 0 .and. .not. EddyFlowProj%wpl .and. .not. cecWplWarned) then
        cecWplWarned = .true.
        call ExceptionHandler(114)
    end if

    !> Suffixes last, in two passes: whether a name needs an occurrence number
    !> depends on how many pairings share the analyser, which the first of them
    !> cannot know.
    do i = 1, npairs
        base = CecPairBaseTag(pairs(i)%carbon_slot, i)
        repeat = 0
        occurrence = 0
        do k = 1, npairs
            if (trim(CecPairBaseTag(pairs(k)%carbon_slot, k)) /= trim(base)) cycle
            repeat = repeat + 1
            if (k <= i) occurrence = occurrence + 1
        end do
        if (repeat == 1) then
            pairs(i)%tag = '_' // trim(base)
        else
            write(pairs(i)%tag, '(a,i0)') '_' // trim(base) // '_', occurrence
        end if
    end do

contains

    !> Do the two slots sit on the same analyser? Asked of the project rather
    !> than of E2Col, like everything else here, because the header builders
    !> want an answer before the first period is imported.
    logical function CecSameAnalyser(a, b)
        integer, intent(in) :: a
        integer, intent(in) :: b
        character(32) :: ia, ib
        integer :: ra, rb

        CecSameAnalyser = .true.
        ra = a - firstGas + 1
        rb = b - firstGas + 1
        if (ra < 1 .or. ra > min(EddyFlowProj%gas_num, MaxNumGases)) return
        if (rb < 1 .or. rb > min(EddyFlowProj%gas_num, MaxNumGases)) return
        ia = EddyFlowProj%gas(ra)%instr
        ib = EddyFlowProj%gas(rb)%instr
        call lowercase(ia)
        call lowercase(ib)
        !> An unnamed analyser is not evidence of two, so it is not warned
        !> about: 'none' and 'other' are what a project says when it has not
        !> told us which instrument a channel is on.
        if (len_trim(ia) == 0 .or. trim(ia) == 'none' .or. trim(ia) == 'other') return
        if (len_trim(ib) == 0 .or. trim(ib) == 'none' .or. trim(ib) == 'other') return
        CecSameAnalyser = trim(ia) == trim(ib)
    end function CecSameAnalyser

    !> The gas slot a 1-based record index names, if it names the wanted
    !> species and a column. Zero for "not stated", which the caller resolves.
    integer function CecSlotOfRecord(rec, wanted)
        integer, intent(in) :: rec
        character(*), intent(in) :: wanted
        character(32) :: species

        CecSlotOfRecord = 0
        if (rec <= 0) return
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) return
        if (EddyFlowProj%gas(rec)%col <= 0) return
        species = EddyFlowProj%gas(rec)%var
        call uppercase(species)
        if (trim(adjustl(species)) /= trim(adjustl(wanted))) return
        CecSlotOfRecord = firstGas + rec - 1
    end function CecSlotOfRecord

    logical function CecPairIsRepeated(seen, nseen, csl, wsl)
        type(CECResolvedPairType), intent(in) :: seen(MaxNumCecPairs)
        integer, intent(in) :: nseen
        integer, intent(in) :: csl
        integer, intent(in) :: wsl
        integer :: n

        CecPairIsRepeated = .false.
        do n = 1, nseen
            if (seen(n)%carbon_slot == csl .and. seen(n)%water_slot == wsl) then
                CecPairIsRepeated = .true.
                return
            end if
        end do
    end function CecPairIsRepeated

    !> The analyser the CO2 sits on, from the project rather than from E2Col:
    !> the header writers ask for these before the first period is imported,
    !> and E2Col is not populated until one is. A pairing whose analyser has no
    !> name falls back to its ordinal, which is unlovely but unambiguous.
    character(48) function CecPairBaseTag(csl, ordinal)
        integer, intent(in) :: csl
        integer, intent(in) :: ordinal
        integer :: r

        CecPairBaseTag = ''
        r = csl - firstGas + 1
        if (r >= 1 .and. r <= min(EddyFlowProj%gas_num, MaxNumGases)) &
            CecPairBaseTag = EddyFlowProj%gas(r)%instr
        call lowercase(CecPairBaseTag)
        if (len_trim(CecPairBaseTag) == 0 &
            .or. trim(CecPairBaseTag) == 'none' &
            .or. trim(CecPairBaseTag) == 'other') &
            write(CecPairBaseTag, '(a,i0)') 'pair', ordinal
    end function CecPairBaseTag

end subroutine CecPairs

!***************************************************************************
!
! \brief       The pairings a project means when it states none.
! \author      Jonathan Muller
! \note        One per CO2 channel, each with the water on the same analyser.
!              Falling back to the designated hygrometer when an analyser
!              carries none is a stand-in, not a measurement of that cell, so
!              it is warned about - the octants are built from q' and c'
!              together, and across two analysers those two arrive at
!              different time lags and through different spectral responses.
!***************************************************************************
subroutine AutoCecPairs(pairs, npairs)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(out) :: pairs(MaxNumCecPairs)
    integer, intent(out) :: npairs

    integer :: gas
    integer :: rec
    character(32) :: species
    integer, external :: CecWaterOnAnalyserOf

    npairs = 0
    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (EddyFlowProj%gas(rec)%col <= 0) cycle
        species = EddyFlowProj%gas(rec)%var
        call uppercase(species)
        if (trim(adjustl(species)) /= 'CO2') cycle
        if (npairs >= MaxNumCecPairs) exit

        npairs = npairs + 1
        pairs(npairs)%meth = 1
        pairs(npairs)%carbon_slot = gas
        pairs(npairs)%water_slot = CecWaterOnAnalyserOf(gas)
        pairs(npairs)%n_extra = 0
        pairs(npairs)%extra_slot = 0
        pairs(npairs)%tag = ''
    end do
end subroutine AutoCecPairs

!***************************************************************************
!
! \brief       The water channel on the same analyser as a given gas slot.
! \author      Jonathan Muller
! \note        Falls back to the designated hygrometer, and says so, because a
!              pairing across two analysers is a real answer with a real caveat
!              rather than an error.
!***************************************************************************
integer function CecWaterOnAnalyserOf(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot

    integer :: gas
    integer :: rec
    integer :: own
    character(32) :: analyser
    character(32) :: other
    character(32) :: species
    integer, external :: PrimaryWaterSlot

    CecWaterOnAnalyserOf = 0
    if (gas_slot < firstGas .or. gas_slot > lastGas) return

    own = gas_slot - firstGas + 1
    if (own < 1 .or. own > min(EddyFlowProj%gas_num, MaxNumGases)) return
    analyser = EddyFlowProj%gas(own)%instr
    call lowercase(analyser)

    if (len_trim(analyser) > 0 .and. trim(analyser) /= 'none' &
        .and. trim(analyser) /= 'other') then
        do gas = firstGas, lastGas
            rec = gas - firstGas + 1
            if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
            if (EddyFlowProj%gas(rec)%col <= 0) cycle
            species = EddyFlowProj%gas(rec)%var
            call uppercase(species)
            if (trim(adjustl(species)) /= 'H2O') cycle
            other = EddyFlowProj%gas(rec)%instr
            call lowercase(other)
            if (trim(other) /= trim(analyser)) cycle
            CecWaterOnAnalyserOf = gas
            return
        end do
    end if

    !> No water on that analyser, so the designated hygrometer stands in. The
    !> octants are built from q' and c' together, and across two analysers
    !> those two arrive at different time lags and through different spectral
    !> responses - a real answer with a real caveat, so say so.
    CecWaterOnAnalyserOf = PrimaryWaterSlot()
end function CecWaterOnAnalyserOf

!***************************************************************************
!
! \brief       The custom column carrying a gas's own signal-strength
!              diagnostic, or 0 if its analyser reports none.
! \author      Jonathan Muller
! \note        Matched on the instrument's metadata block number, not on a
!              model substring. The conditional eddy covariance screen used to
!              take the first column named AGC or RSSI on anything matching
!              li7200 or li7500 and apply it to both gases of the pair - so a
!              site running two analysers screened one analyser's gas with the
!              other's diagnostic, and a site running neither, a quantum
!              cascade laser measuring carbonyl sulfide for instance, was not
!              screened at all.
!***************************************************************************
integer function CecSignalColumnFor(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot

    integer :: j

    CecSignalColumnFor = 0
    if (gas_slot < firstGas .or. gas_slot > lastGas) return
    if (.not. allocated(UserCol)) return
    if (E2Col(gas_slot)%instr%slot <= 0) return

    do j = 1, NumUserVar
        if (UserCol(j)%var /= 'AGC' .and. UserCol(j)%var /= 'RSSI') cycle
        if (UserCol(j)%instr%slot /= E2Col(gas_slot)%instr%slot) cycle
        CecSignalColumnFor = j
        return
    end do
end function CecSignalColumnFor

!***************************************************************************
!
! \brief       Is that diagnostic a received-signal strength, where high is
!              clean, or an automatic gain control, where high is dirty?
! \author      Jonathan Muller
! \note        The two are stored under the same name and compared against the
!              same number, and they run in opposite directions. A column the
!              user labelled RSSI is one. A column labelled AGC is one too from
!              analyser firmware 5.3.0, which is where LI-COR changed what the
!              diagnostic word encodes - the same boundary interpret_diagnostics
!              uses to decode it.
!
!              An analyser that states no firmware compares older than any
!              threshold, so it is read as a true AGC. That is the arm an
!              unpopulated sw_ver has always taken elsewhere in the engine, and
!              it is the safer reading: treating a real AGC as an RSSI keeps
!              exactly the samples that should go.
!***************************************************************************
logical function CecSignalIsRssi(gas_slot, user_col)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    integer, intent(in) :: user_col

    logical, external :: CompareSwVer
    type(SwVerType), external :: SwVerFromString

    CecSignalIsRssi = .false.
    if (.not. allocated(UserCol)) return
    if (user_col < 1 .or. user_col > NumUserVar) return
    if (UserCol(user_col)%var == 'RSSI') then
        CecSignalIsRssi = .true.
        return
    end if
    if (UserCol(user_col)%var /= 'AGC') return
    if (gas_slot < firstGas .or. gas_slot > lastGas) return

    CecSignalIsRssi = CompareSwVer(E2Col(gas_slot)%instr%sw_ver, &
        SwVerFromString('5.3.0'))
end function CecSignalIsRssi

!***************************************************************************
!
! \brief       The scalars one pairing partitions, in output order.
! \author      Jonathan Muller
! \note        Water first, carbon second, extras after, always. Header
!              builders, row writers, the essentials record and this module all
!              walk the same list from the same helper, because a target list
!              that is right in one place and wrong in another shifts a file
!              rather than failing to build.
!***************************************************************************
subroutine CecTargetSlots(pair, slots, ntarget)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(in) :: pair
    integer, intent(out) :: slots(MaxNumCecTargets)
    integer, intent(out) :: ntarget

    integer :: k

    slots = 0
    slots(cecTargetWater) = pair%water_slot
    slots(cecTargetCarbon) = pair%carbon_slot
    ntarget = 2
    do k = 1, pair%n_extra
        if (pair%extra_slot(k) <= 0) cycle
        if (ntarget >= MaxNumCecTargets) exit
        ntarget = ntarget + 1
        slots(ntarget) = pair%extra_slot(k)
    end do
end subroutine CecTargetSlots

!***************************************************************************
!
! \brief       The full-output columns one pairing contributes, named and
!              with their units, in the order the row writers emit them.
! \author      Jonathan Muller
! \note        One helper for both, because the header and the row are built
!              in different files and, for the FCC executable, a different
!              program. A name list that is right in one and wrong in the other
!              does not fail to build; it shifts every column after it.
!
!              The order is: the pairing's targets, water first and carbon
!              second and extras after, then the octant statistics the ratios
!              were computed from. Water and carbon are named for what they are
!              - evaporation, transpiration, respiration, photosynthesis - and
!              anything else generically, because the octants sort by pathway
!              rather than by medium and "soil" would be a guess about a
!              species the method knows nothing about.
!***************************************************************************
subroutine CecOutputColumns(pair, flux_label, names, units, nnames)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(in) :: pair
    character(*), intent(in) :: flux_label(GHGNumVar)
    character(*), intent(out) :: names(:)
    character(*), intent(out) :: units(:)
    integer, intent(out) :: nnames

    integer :: k
    integer :: ntarget
    integer :: slots(MaxNumCecTargets)
    character(64) :: stem
    character(64) :: gas_tags(GHGNumVar)
    character(48) :: tag

    names = ''
    units = ''
    nnames = 0
    if (pair%meth == 0) return

    call CecTargetSlots(pair, slots, ntarget)
    call FullOutputGasTags(gas_tags)
    tag = pair%tag

    do k = 1, ntarget
        if (slots(k) < firstGas .or. slots(k) > lastGas) cycle
        if (k == cecTargetWater) then
            if (pair%meth /= 1 .and. pair%meth /= 2) cycle
            call Emit('E_cec' // trim(tag), flux_label(slots(k)))
            call Emit('Tr_cec' // trim(tag), flux_label(slots(k)))
            call Emit('E_cec_ET' // trim(tag), '[mm+1hour-1]')
            call Emit('Tr_cec_ET' // trim(tag), '[mm+1hour-1]')
            call Emit('ET_cec' // trim(tag), flux_label(slots(k)))
            call Emit('r_ET_cec' // trim(tag), '[#]')
            call Emit('qc_cec_h2o' // trim(tag), '[#]')
        else if (k == cecTargetCarbon) then
            if (pair%meth /= 1 .and. pair%meth /= 3) cycle
            call Emit('Reco_cec' // trim(tag), flux_label(slots(k)))
            call Emit('P_cec' // trim(tag), flux_label(slots(k)))
            call Emit('NEE_cec' // trim(tag), flux_label(slots(k)))
            call Emit('r_Fc_cec' // trim(tag), '[#]')
            call Emit('qc_cec_co2' // trim(tag), '[#]')
        else
            stem = gas_tags(slots(k))
            if (len_trim(stem) > 1) stem = stem(1:len_trim(stem) - 1)
            call Emit(trim(stem) // '_nonstomatal_cec' // trim(tag), &
                flux_label(slots(k)))
            call Emit(trim(stem) // '_stomatal_cec' // trim(tag), &
                flux_label(slots(k)))
            call Emit(trim(stem) // '_total_cec' // trim(tag), &
                flux_label(slots(k)))
            call Emit('r_' // trim(stem) // '_cec' // trim(tag), '[#]')
            call Emit('qc_cec_' // trim(stem) // trim(tag), '[#]')
        end if
    end do

    !> How many points each octant held, and what fraction of the record that
    !> was. These decide whether a ratio is worth anything, and they used to be
    !> written only to the essentials file - so the one output a user reads
    !> could not tell a confident partition from a badly sampled one.
    call Emit('cec_n_o1' // trim(tag), '[#]')
    call Emit('cec_n_o2' // trim(tag), '[#]')
    call Emit('cec_frac_o1' // trim(tag), '[#]')
    call Emit('cec_frac_o2' // trim(tag), '[#]')

contains

    subroutine Emit(name, unit)
        character(*), intent(in) :: name
        character(*), intent(in) :: unit

        if (nnames >= size(names)) return
        nnames = nnames + 1
        names(nnames) = name
        units(nnames) = unit
    end subroutine Emit

end subroutine CecOutputColumns

!***************************************************************************
!
! \brief       The values behind CecOutputColumns, in the same order.
! \author      Jonathan Muller
! \note        Deliberately next to CecOutputColumns and deliberately walking
!              the same branches in the same order: these two are one contract
!              expressed twice, and a row that emits its values in a different
!              order than its header names them is a file whose every later
!              column is misread rather than a build that fails.
!
!              `is_int` says which are counts and status codes, so the writers
!              format them as integers without having to know which is which.
!
!              The scale is the caller's, for the same reason the labels are
!              the caller's in CecOutputColumns: RP resolves a gas's unit from
!              the metadata columns it read, FCC from its own copy taken at
!              startup, and a shared helper asking GasUnitIn is right only in
!              RP. That is what had FCC labelling carbonyl sulfide in umol and
!              printing it unscaled - a thousandfold error in the one output a
!              reader actually reads.
!***************************************************************************
subroutine CecRowValues(pair, descriptor, flux, flux_sc, values, is_int, nvalues)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(in) :: pair
    type(CECDescriptorType), intent(in) :: descriptor
    type(CECFluxType), intent(in) :: flux
    real(kind = dbl), intent(in) :: flux_sc(GHGNumVar)
    real(kind = dbl), intent(out) :: values(:)
    logical, intent(out) :: is_int(:)
    integer, intent(out) :: nvalues

    integer :: k
    integer :: ntarget
    integer :: slots(MaxNumCecTargets)
    real(kind = dbl) :: sc

    values = error
    is_int = .false.
    nvalues = 0
    if (pair%meth == 0) return

    call CecTargetSlots(pair, slots, ntarget)

    do k = 1, ntarget
        if (slots(k) < firstGas .or. slots(k) > lastGas) cycle
        sc = flux_sc(slots(k))
        if (k == cecTargetWater) then
            if (pair%meth /= 1 .and. pair%meth /= 2) cycle
            call EmitScaled(flux%comp(k)%nonstomatal, sc)
            call EmitScaled(flux%comp(k)%stomatal, sc)
            call EmitScaled(flux%E_cec_ET, 1d0)
            call EmitScaled(flux%Tr_cec_ET, 1d0)
            call EmitScaled(flux%comp(k)%total, sc)
            call EmitScaled(descriptor%target(k)%r, 1d0)
            call EmitCount(flux%comp(k)%status)
        else if (k == cecTargetCarbon) then
            if (pair%meth /= 1 .and. pair%meth /= 3) cycle
            call EmitScaled(flux%comp(k)%nonstomatal, sc)
            call EmitScaled(flux%comp(k)%stomatal, sc)
            call EmitScaled(flux%comp(k)%total, sc)
            call EmitScaled(descriptor%target(k)%r, 1d0)
            call EmitCount(flux%comp(k)%status)
        else
            call EmitScaled(flux%comp(k)%nonstomatal, sc)
            call EmitScaled(flux%comp(k)%stomatal, sc)
            call EmitScaled(flux%comp(k)%total, sc)
            call EmitScaled(descriptor%target(k)%r, 1d0)
            call EmitCount(flux%comp(k)%status)
        end if
    end do

    call EmitCount(descriptor%n_O1)
    call EmitCount(descriptor%n_O2)
    call EmitScaled(descriptor%frac_O1, 1d0)
    call EmitScaled(descriptor%frac_O2, 1d0)

contains

    subroutine EmitScaled(value, scale)
        real(kind = dbl), intent(in) :: value
        real(kind = dbl), intent(in) :: scale

        if (nvalues >= size(values)) return
        nvalues = nvalues + 1
        if (value == error) then
            values(nvalues) = error
        else
            values(nvalues) = value * scale
        end if
        is_int(nvalues) = .false.
    end subroutine EmitScaled

    subroutine EmitCount(value)
        integer, intent(in) :: value

        if (nvalues >= size(values)) return
        nvalues = nvalues + 1
        values(nvalues) = dble(value)
        is_int(nvalues) = .true.
    end subroutine EmitCount

end subroutine CecRowValues

!***************************************************************************
!
! \brief       One pairing's descriptor as the essentials row carries it.
! \author      Jonathan Muller
! \note        The transport from RP to FCC. Three writers emit this - the
!              normal row, the skipped-period row, and FCC's echo of both - and
!              ReadExRecord consumes it, so it is built once here.
!
!              The widths are fixed and named: nine fields for the pairing,
!              then six for each of the targets the NINTH of those announces.
!              The count comes from the pairing rather than from the
!              descriptor, because a period whose extraction bailed still has
!              to occupy its own width in the row.
!***************************************************************************
subroutine CecExRowValues(pair, descriptor, values, is_int, nvalues)
    use m_common_global_var
    implicit none
    type(CECResolvedPairType), intent(in) :: pair
    type(CECDescriptorType), intent(in) :: descriptor
    real(kind = dbl), intent(out) :: values(:)
    logical, intent(out) :: is_int(:)
    integer, intent(out) :: nvalues

    integer :: k
    integer :: ntarget
    integer :: slots(MaxNumCecTargets)

    values = error
    is_int = .false.
    nvalues = 0

    call CecTargetSlots(pair, slots, ntarget)

    call EmitInt(pair%meth)
    call EmitInt(pair%carbon_slot)
    call EmitInt(pair%water_slot)
    call EmitInt(descriptor%n_valid)
    call EmitInt(descriptor%n_O1)
    call EmitInt(descriptor%n_O2)
    call EmitReal(descriptor%frac_O1)
    call EmitReal(descriptor%frac_O2)
    call EmitInt(ntarget)

    do k = 1, ntarget
        call EmitInt(slots(k))
        call EmitReal(descriptor%target(k)%f_O1)
        call EmitReal(descriptor%target(k)%f_O2)
        call EmitReal(descriptor%target(k)%r)
        call EmitInt(descriptor%target(k)%status)
        call EmitInt(merge(1, 0, descriptor%target(k)%valid))
    end do

contains

    subroutine EmitInt(value)
        integer, intent(in) :: value

        if (nvalues >= size(values)) return
        nvalues = nvalues + 1
        values(nvalues) = dble(value)
        is_int(nvalues) = .true.
    end subroutine EmitInt

    subroutine EmitReal(value)
        real(kind = dbl), intent(in) :: value

        if (nvalues >= size(values)) return
        nvalues = nvalues + 1
        values(nvalues) = value
        is_int(nvalues) = .false.
    end subroutine EmitReal

end subroutine CecExRowValues
