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
! \brief       The gas slot holding the site's primary water measurement.
! \author      Jonathan Muller
! \note        A site has one humidity, one latent heat flux and one
!              evapotranspiration however many hygrometers it carries, and
!              those come from the first water record. Gases that need *their
!              own* water for a WPL correction use E2Col(gas)%moist_ref
!              instead; this is only for the one-per-site quantities.
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
    integer :: gas
    logical, external :: GasSlotIsWater

    PrimaryWaterSlot = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (.not. GasSlotIsWater(gas)) cycle
        if (EddyFlowProj%gas(gas - firstGas + 1)%col <= 0) cycle
        PrimaryWaterSlot = gas
        return
    end do
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
!              The full output has two header branches. The dynamic one names
!              one block per present gas; the fix_out_format one is a literal
!              naming exactly co2, h2o, ch4 and the fourth slot. Both row
!              writers, however, loop firstGas..lastGas and - this is the part
!              that bites - their `elseif (fix_out_format)` arms emit
!              placeholder fields for slots with *no* gas at all. While
!              lastGas was 8 the two agreed. At 68 the row carries sixty
!              phantom gas blocks: fifteen fields each, nine hundred fields
!              past the end of a header that describes a hundred and ninety
!              four. Everything after the gas block is shifted, and FCC parses
!              the ex record by comma count, so it consumes the row without
!              complaint. No fixture set fix_out_format, which is why it went
!              unseen.
!
!              fix_out_format is a compatibility mode: it promises the fixed
!              EddyPro 7.x column set, so it returns exactly the four
!              historical slots whatever the project configures - present or
!              not, since the row fills absent ones with the error label and
!              the header names them unconditionally. A fifth gas cannot be
!              represented in that format and is dropped from this file; the
!              caller warns once. Widening it instead would break the
!              compatibility the flag exists to provide.
!
!              Otherwise the whole gas block, which is what both sides walk
!              today. The OutVarPresent guard at each of the eight use sites
!              does the narrowing, and it is the *same* guard in the header
!              and in the row, so they cannot disagree.
!***************************************************************************
subroutine FullOutputGasSlots(slots, nslots)
    use m_common_global_var
    implicit none
    integer, intent(out) :: slots(GHGNumVar)
    integer, intent(out) :: nslots
    integer :: gas

    slots = 0
    nslots = 0

    if (EddyFlowProj%fix_out_format) then
        do gas = histGas1, histGas4
            nslots = nslots + 1
            slots(nslots) = gas
        end do
        return
    end if

    do gas = firstGas, lastGas
        nslots = nslots + 1
        slots(nslots) = gas
    end do
end subroutine FullOutputGasSlots

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
    integer :: gas
    character(32) :: species

    PrimaryCarbonSlot = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (EddyFlowProj%gas(gas - firstGas + 1)%col <= 0) cycle
        species = EddyFlowProj%gas(gas - firstGas + 1)%var
        call uppercase(species)
        if (trim(adjustl(species)) /= 'CO2') cycle
        PrimaryCarbonSlot = gas
        return
    end do
end function PrimaryCarbonSlot

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
