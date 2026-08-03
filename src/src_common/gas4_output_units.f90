!***************************************************************************
! gas4_output_units.f90
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
! \brief       Per-gas output scales, labels and column-name stems.
! \author      Jonathan Muller
! \note        The file name is historical: none of this is specific to the
!              fourth slot any more. Every routine here answers for whichever
!              gas slot it is asked about, from the gas record.
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
        do gas = co2, gas4
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
        GasSlotFromDynMDTag = gas4
        return
    end if
    if (IsHistoricGasVar(stem)) GasSlotFromDynMDTag = HistoricGasSlot(stem)
end function GasSlotFromDynMDTag

integer function PrimaryWaterOutSlot()
    use m_common_global_var
    implicit none
    integer, external :: PrimaryWaterSlot

    PrimaryWaterOutSlot = PrimaryWaterSlot()
    if (PrimaryWaterOutSlot < firstGas) PrimaryWaterOutSlot = h2o
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
    integer :: gas, k, repeat
    character(32), external :: GasOutputLabel

    tags = ''
    nseen = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        label = GasOutputLabel(gas)
        call lowercase(label)
        if (len_trim(label) == 0) cycle

        repeat = 1
        do k = 1, nseen
            if (trim(seen(k)) == trim(label)) repeat = repeat + 1
        end do
        nseen = nseen + 1
        seen(nseen) = label

        if (repeat == 1) then
            tags(gas) = trim(label) // '_'
        else
            write(tags(gas), '(a,i0,a)') trim(label) // '_', repeat, '_'
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
        case (co2);  LegacySpectralVarTag = 'co2'
        case (h2o);  LegacySpectralVarTag = 'h2o'
        case (ch4);  LegacySpectralVarTag = 'ch4'
        case (gas4); LegacySpectralVarTag = 'gas4'
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
