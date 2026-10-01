!***************************************************************************
! sa_rates.f90
! ------------
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
! \brief       Spectral assessment per acquisition frequency.
! \author      Jonathan Muller
! \note        A project's raw files need not all be at one rate: a GHG project
!              can go from 10 to 20 Hz, and an analyser can change its own rate
!              while the file's stays - a QCLS at 1 Hz and later at 0.33 Hz in
!              a 20 Hz file. Each gas's spectra then mean different things in
!              different periods: the noise floor, the Nyquist frequency and the
!              frequencies a transfer function can be fitted over all depend on
!              the rate. So each gas is assessed separately at each of its
!              rates, and each period is corrected with the result for its own.
!
!              A period's RATE CONFIGURATION is its file rate plus every gas's
!              effective rate - the analyser's own where it states one, capped
!              at the file's, as RP's ColumnAcFreq has it. The assessment runs
!              one pass per configuration. A gas's ensemble in a pass pools
!              every configuration in which that gas has the same rate, so an
!              analyser that never changed is fitted once from all its periods,
!              not split between passes. Results are kept per gas and per rate
!              of that gas (RATE SLOTS, fastest first); the Ibrom et al. (2007)
!              correction-factor model, which comes from w'T', per file rate.
!
!              With one configuration nothing here is used beyond the slot
!              lookup: the assessment runs exactly as it always has.
!***************************************************************************
module m_sa_rates
    use m_fx_global_var
    implicit none
    save

    integer, parameter :: MaxRateConfigs = 8
    integer, parameter :: MaxRateSlots = 8

    !> Distinct rate configurations over the whole ex file
    integer :: nRateConfigs = 0
    real(kind = dbl) :: ConfigFileRate(MaxRateConfigs) = -1d0
    real(kind = dbl) :: ConfigGasRate(GHGNumVar, MaxRateConfigs) = -1d0
    integer :: ConfigPeriods(MaxRateConfigs) = 0
    logical :: ConfigOverflowWarned = .false.

    !> Each gas's distinct rates, fastest first, and the periods at each
    integer :: nGasRates(GHGNumVar) = 0
    real(kind = dbl) :: GasRate(GHGNumVar, MaxRateSlots) = -1d0
    integer :: GasRatePeriods(GHGNumVar, MaxRateSlots) = 0

    !> The file's distinct rates, fastest first, for the Ibrom model
    integer :: nFileRates = 0
    real(kind = dbl) :: FileRate(MaxRateSlots) = -1d0

    !> More than one configuration: the assessment runs per configuration
    logical :: MultiRateSA = .false.
    !> The configuration being assessed. 1 outside the per-rate passes, which
    !> for a single-configuration project is the only one there is.
    integer :: SAPassConfig = 1

    !> Per-configuration spectral sums, filled beside MeanBinSpec
    type(MeanSpectraType), allocatable :: ConfigBinSpec(:, :, :)

    !> Results per gas and rate slot: the whole RegPar row, which carries a
    !> water gas's RH exponential in its dum element. SlotCnt is the ensemble
    !> size each class was fitted from, for the file's numerosity column.
    type(RegParType) :: RegParR(GHGNumVar, MaxGasClasses, MaxRateSlots)
    integer :: SlotCnt(GHGNumVar, MaxGasClasses, MaxRateSlots) = 0
    logical :: SlotFilled(GHGNumVar, MaxRateSlots) = .false.
    !> The slot a rate slot is actually corrected with, after fallbacks.
    !> 0 when no slot of that gas is usable.
    integer :: GasRateUse(GHGNumVar, MaxRateSlots) = 0
    real(kind = dbl) :: UnParR(2, MaxRateSlots) = -9999d0
    real(kind = dbl) :: StParR(2, MaxRateSlots) = -9999d0
    logical :: FileSlotFilled(MaxRateSlots) = .false.
    integer :: FileRateUse(MaxRateSlots) = 0

    !> Output control for OutputSpectralAssessmentResults: during the passes
    !> only a configuration's own averaged-spectra files are written, under a
    !> suffix naming its rates; the assessment file is written once, after.
    logical :: SAOutTxt = .true.
    character(256) :: SAOutSuffix = ''

contains

    !***************************************************************************
    !> Two rates are the same rate
    !***************************************************************************
    logical function SameRate(r1, r2)
        real(kind = dbl), intent(in) :: r1, r2
        SameRate = abs(r1 - r2) <= 1d-6 * max(abs(r1), abs(r2))
    end function SameRate

    !***************************************************************************
    !> The rate a gas's column was sampled at in this record: its analyser's,
    !> capped at the file's - RP's ColumnAcFreq - or the file's when the
    !> analyser states none.
    !***************************************************************************
    real(kind = dbl) function EffectiveGasRate(lEx, gas)
        type(ExType), intent(in) :: lEx
        integer, intent(in) :: gas

        EffectiveGasRate = lEx%ac_freq
        if (lEx%gas_instr(gas)%ac_freq > 0d0) then
            if (lEx%ac_freq > 0d0) then
                EffectiveGasRate = min(lEx%gas_instr(gas)%ac_freq, lEx%ac_freq)
            else
                EffectiveGasRate = lEx%gas_instr(gas)%ac_freq
            end if
        end if
    end function EffectiveGasRate

    !***************************************************************************
    !> The Nyquist frequency of `gas` in the configuration being assessed, or
    !> a negative value when its rate is not known.
    !***************************************************************************
    real(kind = dbl) function GasNyquist(gas)
        integer, intent(in) :: gas
        GasNyquist = -1d0
        if (SAPassConfig < 1 .or. SAPassConfig > max(nRateConfigs, 1)) return
        if (ConfigGasRate(gas, SAPassConfig) > 0d0) &
            GasNyquist = ConfigGasRate(gas, SAPassConfig) / 2d0
    end function GasNyquist

    !***************************************************************************
    !> A gas's rate in this record for telling configurations apart: negative
    !> - unknown - when the gas has no data in the period. RP resolves an
    !> absent column's rate to the file's, and a QCLS channel that was not
    !> logged for a half-hour would otherwise open a configuration of its own
    !> at the sonic's rate.
    !***************************************************************************
    real(kind = dbl) function RecordGasRate(lEx, gas)
        type(ExType), intent(in) :: lEx
        integer, intent(in) :: gas

        RecordGasRate = -1d0
        if (.not. lEx%var_present(gas)) return
        RecordGasRate = EffectiveGasRate(lEx, gas)
    end function RecordGasRate

    !***************************************************************************
    !> Whether gas slot `gas` is one the project configures
    !***************************************************************************
    logical function ConfiguredGas(gas)
        integer, intent(in) :: gas
        ConfiguredGas = gas >= firstGas .and. gas <= lastGas .and. &
            gas - firstGas + 1 <= min(EddyFlowProj%gas_num, MaxNumGases)
    end function ConfiguredGas

    !***************************************************************************
    !> Insert rate r into a fastest-first list, returning its slot
    !***************************************************************************
    subroutine InsertRate(list, n, r, slot)
        real(kind = dbl), intent(inout) :: list(MaxRateSlots)
        integer, intent(inout) :: n
        real(kind = dbl), intent(in) :: r
        integer, intent(out) :: slot
        integer :: i, j

        do i = 1, n
            if (SameRate(list(i), r)) then
                slot = i
                return
            end if
        end do
        if (n == MaxRateSlots) then
            !> Full: the nearest rate stands in
            slot = NearestSlot(list, n, r)
            return
        end if
        !> Keep fastest first
        i = n + 1
        do j = 1, n
            if (r > list(j)) then
                i = j
                exit
            end if
        end do
        list(i + 1:n + 1) = list(i:n)
        list(i) = r
        n = n + 1
        slot = i
    end subroutine InsertRate

    !***************************************************************************
    !> The slot for rate r: the same rate, else the nearest faster one, else
    !> the nearest slower one. 0 when the list is empty.
    !***************************************************************************
    integer function NearestSlot(list, n, r)
        real(kind = dbl), intent(in) :: list(MaxRateSlots)
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: r
        integer :: i

        NearestSlot = 0
        if (n <= 0) return
        do i = 1, n
            if (SameRate(list(i), r)) then
                NearestSlot = i
                return
            end if
        end do
        !> Fastest first, so the last faster one is the nearest faster
        do i = n, 1, -1
            if (list(i) > r) then
                NearestSlot = i
                return
            end if
        end do
        !> None faster: the fastest of the slower ones is the first
        NearestSlot = 1
    end function NearestSlot

    integer function GasRateSlot(gas, r)
        integer, intent(in) :: gas
        real(kind = dbl), intent(in) :: r
        GasRateSlot = NearestSlot(GasRate(gas, :), nGasRates(gas), r)
    end function GasRateSlot

    integer function FileRateSlot(r)
        real(kind = dbl), intent(in) :: r
        FileRateSlot = NearestSlot(FileRate, nFileRates, r)
    end function FileRateSlot

    !***************************************************************************
    !> Record the rates of one valid ex record. Called for every one, by
    !> InitExVars.
    !***************************************************************************
    subroutine RegisterRecordRates(lEx)
        type(ExType), intent(in) :: lEx
        integer :: c, gas, slot
        real(kind = dbl) :: r(GHGNumVar)

        if (lEx%ac_freq <= 0d0) return

        r = -1d0
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            r(gas) = RecordGasRate(lEx, gas)
        end do

        c = MatchConfig(lEx%ac_freq, r)
        if (c > 0) then
            !> Matched on a gas this configuration had not seen yet: it has one
            !> now. Two configurations this makes identical are merged at the end.
            do gas = firstGas, lastGas
                if (r(gas) > 0d0 .and. ConfigGasRate(gas, c) <= 0d0) &
                    ConfigGasRate(gas, c) = r(gas)
            end do
        else
            if (nRateConfigs < MaxRateConfigs) then
                nRateConfigs = nRateConfigs + 1
                c = nRateConfigs
                ConfigFileRate(c) = lEx%ac_freq
                ConfigGasRate(:, c) = r
            else
                !> A ninth arrangement of rates is merged into the nearest one
                !> by file rate. Said once; the periods are still processed.
                if (.not. ConfigOverflowWarned) then
                    call LogSay(' Warning: more than 8 combinations of acquisition &
                        &frequencies; further ones are assessed with the nearest.')
                    ConfigOverflowWarned = .true.
                end if
                c = NearestConfig(lEx%ac_freq)
            end if
        end if
        ConfigPeriods(c) = ConfigPeriods(c) + 1

        !> Periods per gas rate are counted once the lists are final - an
        !> insertion here shifts the slots after it
        call InsertRate(FileRate, nFileRates, lEx%ac_freq, slot)
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            if (r(gas) <= 0d0) cycle
            call InsertRate(GasRate(gas, :), nGasRates(gas), r(gas), slot)
        end do
    end subroutine RegisterRecordRates

    !***************************************************************************
    !> The configuration a set of rates belongs to. A rate that is unknown on
    !> either side - a gas with no data in the period, or one the configuration
    !> has not seen - does not tell configurations apart. An exact match is
    !> preferred over one that relies on that.
    !***************************************************************************
    integer function MatchConfig(fr, r)
        real(kind = dbl), intent(in) :: fr
        real(kind = dbl), intent(in) :: r(GHGNumVar)
        integer :: c, pass

        MatchConfig = 0
        do pass = 1, 2
            do c = 1, nRateConfigs
                if (.not. SameRate(ConfigFileRate(c), fr)) cycle
                if (RatesAgree(ConfigGasRate(:, c), r, pass == 1)) then
                    MatchConfig = c
                    return
                end if
            end do
        end do
    end function MatchConfig

    !> Two sets of gas rates agree: every gas at the same rate, or - unless
    !> `exact` - unknown on one side
    logical function RatesAgree(r1, r2, exact)
        real(kind = dbl), intent(in) :: r1(GHGNumVar), r2(GHGNumVar)
        logical, intent(in) :: exact
        integer :: gas

        RatesAgree = .false.
        do gas = firstGas, lastGas
            if (r1(gas) <= 0d0 .and. r2(gas) <= 0d0) cycle
            if (r1(gas) <= 0d0 .or. r2(gas) <= 0d0) then
                if (exact) return
                cycle
            end if
            if (.not. SameRate(r1(gas), r2(gas))) return
        end do
        RatesAgree = .true.
    end function RatesAgree

    integer function NearestConfig(fr)
        real(kind = dbl), intent(in) :: fr
        integer :: c
        NearestConfig = 1
        do c = 1, nRateConfigs
            if (abs(ConfigFileRate(c) - fr) < abs(ConfigFileRate(NearestConfig) - fr)) &
                NearestConfig = c
        end do
    end function NearestConfig

    !***************************************************************************
    !> The configuration of one record
    !***************************************************************************
    integer function RateConfigOf(lEx)
        type(ExType), intent(in) :: lEx
        integer :: gas
        real(kind = dbl) :: r(GHGNumVar)

        RateConfigOf = 1
        if (nRateConfigs <= 1) return
        r = -1d0
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            r(gas) = RecordGasRate(lEx, gas)
        end do
        RateConfigOf = MatchConfig(lEx%ac_freq, r)
        if (RateConfigOf == 0) RateConfigOf = NearestConfig(lEx%ac_freq)
    end function RateConfigOf

    !***************************************************************************
    !> After InitExVars: order the configurations fastest first - so that
    !> configuration 1, whose averaged spectra keep today's file names, is the
    !> fastest - and decide whether the assessment runs per configuration.
    !***************************************************************************
    subroutine FinaliseRateConfigs()
        integer :: c, c2, gas, p
        real(kind = dbl) :: fr, gr(GHGNumVar)

        if (nRateConfigs == 0) then
            nRateConfigs = 1
            ConfigFileRate(1) = FCCMetadata%ac_freq
        end if
        !> Configurations that learnt a gas's rate after they were opened can
        !> have become the same one
        c = 1
        do while (c < nRateConfigs)
            c2 = c + 1
            do while (c2 <= nRateConfigs)
                if (SameRate(ConfigFileRate(c), ConfigFileRate(c2)) .and. &
                    RatesAgree(ConfigGasRate(:, c), ConfigGasRate(:, c2), .true.)) then
                    ConfigPeriods(c) = ConfigPeriods(c) + ConfigPeriods(c2)
                    ConfigFileRate(c2:nRateConfigs - 1) = ConfigFileRate(c2 + 1:nRateConfigs)
                    ConfigGasRate(:, c2:nRateConfigs - 1) = ConfigGasRate(:, c2 + 1:nRateConfigs)
                    ConfigPeriods(c2:nRateConfigs - 1) = ConfigPeriods(c2 + 1:nRateConfigs)
                    nRateConfigs = nRateConfigs - 1
                else
                    c2 = c2 + 1
                end if
            end do
            c = c + 1
        end do
        !> Insertion sort, descending by file rate then by summed gas rates
        do c = 2, nRateConfigs
            fr = ConfigFileRate(c)
            gr = ConfigGasRate(:, c)
            p = ConfigPeriods(c)
            c2 = c - 1
            do while (c2 >= 1)
                if (.not. Faster(fr, gr, ConfigFileRate(c2), ConfigGasRate(:, c2))) exit
                ConfigFileRate(c2 + 1) = ConfigFileRate(c2)
                ConfigGasRate(:, c2 + 1) = ConfigGasRate(:, c2)
                ConfigPeriods(c2 + 1) = ConfigPeriods(c2)
                c2 = c2 - 1
            end do
            ConfigFileRate(c2 + 1) = fr
            ConfigGasRate(:, c2 + 1) = gr
            ConfigPeriods(c2 + 1) = p
        end do
        MultiRateSA = nRateConfigs > 1
        do gas = firstGas, lastGas
            if (nGasRates(gas) == 0 .and. FCCMetadata%ac_freq > 0d0) then
                nGasRates(gas) = 1
                GasRate(gas, 1) = FCCMetadata%ac_freq
            end if
        end do
        if (nFileRates == 0) then
            nFileRates = 1
            FileRate(1) = FCCMetadata%ac_freq
        end if

        GasRatePeriods = 0
        do c = 1, nRateConfigs
            do gas = firstGas, lastGas
                if (.not. ConfiguredGas(gas)) cycle
                if (ConfigGasRate(gas, c) <= 0d0) cycle
                p = GasRateSlot(gas, ConfigGasRate(gas, c))
                if (p > 0) GasRatePeriods(gas, p) = GasRatePeriods(gas, p) + ConfigPeriods(c)
            end do
        end do

        RegParR%Fn = error
        RegParR%fc = error
        RegParR%f2 = error
        RegParR%e1 = error
        RegParR%e2 = error
        RegParR%e3 = error
        UnParR = error
        StParR = error
    end subroutine FinaliseRateConfigs

    logical function Faster(fr1, gr1, fr2, gr2)
        real(kind = dbl), intent(in) :: fr1, gr1(GHGNumVar), fr2, gr2(GHGNumVar)
        if (.not. SameRate(fr1, fr2)) then
            Faster = fr1 > fr2
        else
            Faster = sum(max(gr1, 0d0)) > sum(max(gr2, 0d0)) * (1d0 + 1d-9)
        end if
    end function Faster

    !***************************************************************************
    !> Suffix naming a configuration's rates, for its averaged-spectra files
    !> and its line in the log: the file rate, then each distinct slower rate
    !> a gas is at - _20Hz_1Hz for a 20 Hz file with a 1 Hz analyser in it.
    !***************************************************************************
    character(256) function ConfigSuffix(c)
        integer, intent(in) :: c
        integer :: c2

        ConfigSuffix = BareSuffix(c)
        !> Two configurations with the same set of rates on different gases
        !> would share a name
        do c2 = 1, nRateConfigs
            if (c2 == c) cycle
            if (BareSuffix(c2) == ConfigSuffix) then
                write(ConfigSuffix, '(a, a, i0)') trim(BareSuffix(c)), '_', c
                exit
            end if
        end do
    end function ConfigSuffix

    !> The file rate, then each distinct slower rate of a gas, fastest first
    character(256) function BareSuffix(c)
        integer, intent(in) :: c
        integer :: gas, k, n
        real(kind = dbl) :: seen(GHGNumVar)
        real(kind = dbl) :: r
        logical :: new

        BareSuffix = '_' // trim(RateText(ConfigFileRate(c))) // 'Hz'
        n = 0
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            r = ConfigGasRate(gas, c)
            if (r <= 0d0 .or. SameRate(r, ConfigFileRate(c))) cycle
            new = .true.
            do k = 1, n
                if (SameRate(seen(k), r)) new = .false.
            end do
            if (new) then
                n = n + 1
                seen(n) = r
            end if
        end do
        !> Fastest first
        do k = 1, n
            r = maxval(seen(k:n))
            seen(maxloc(seen(k:n), 1) + k - 1) = seen(k)
            seen(k) = r
            BareSuffix = trim(BareSuffix) // '_' // trim(RateText(r)) // 'Hz'
        end do
    end function BareSuffix

    !***************************************************************************
    !> A rate as text: 20, 0.333 - three decimals, trailing zeros dropped
    !***************************************************************************
    character(16) function RateText(r)
        real(kind = dbl), intent(in) :: r
        integer :: i

        write(RateText, '(f0.3)') r
        if (RateText(1:1) == '.') RateText = '0' // RateText(1:15)
        i = len_trim(RateText)
        do while (i > 1 .and. RateText(i:i) == '0')
            RateText(i:i) = ' '
            i = i - 1
        end do
        if (RateText(i:i) == '.') RateText(i:i) = ' '
    end function RateText

    !***************************************************************************
    !> A gas's rates as a token value: 20.000,10.000
    !***************************************************************************
    character(160) function RateList(gas)
        integer, intent(in) :: gas
        integer :: k
        character(16) :: s

        RateList = ''
        do k = 1, nGasRates(gas)
            write(s, '(f0.3)') GasRate(gas, k)
            if (s(1:1) == '.') s = '0' // s(1:15)
            if (k > 1) RateList = trim(RateList) // ','
            RateList = trim(RateList) // trim(s)
        end do
    end function RateList

    character(160) function FileRateList()
        integer :: k
        character(16) :: s

        FileRateList = ''
        do k = 1, nFileRates
            write(s, '(f0.3)') FileRate(k)
            if (s(1:1) == '.') s = '0' // s(1:15)
            if (k > 1) FileRateList = trim(FileRateList) // ','
            FileRateList = trim(FileRateList) // trim(s)
        end do
    end function FileRateList

    !***************************************************************************
    !> Fill MeanBinSpec with configuration c's ensemble. For each gas, every
    !> configuration in which that gas has the same rate as in c is pooled, so
    !> a gas that never changed rate gets the same ensemble in every pass.
    !> Un-normalised sums, like the accumulation it replaces.
    !***************************************************************************
    subroutine LoadConfigEnsemble(c, nbins)
        integer, intent(in) :: c
        integer, intent(in) :: nbins
        integer :: c2, gas, cls, bin

        MeanBinSpec = NullMeanSpec
        dMeanBinSpec = NullMeanSpec
        do gas = firstGas, lastGas
            do c2 = 1, nRateConfigs
                if (ConfigGasRate(gas, c) > 0d0 .or. ConfigGasRate(gas, c2) > 0d0) then
                    if (.not. SameRate(ConfigGasRate(gas, c2), ConfigGasRate(gas, c))) cycle
                end if
                do cls = 1, MaxGasClasses
                    do bin = 1, nbins
                        MeanBinSpec(bin, cls)%cnt(gas) = MeanBinSpec(bin, cls)%cnt(gas) &
                            + ConfigBinSpec(bin, cls, c2)%cnt(gas)
                        MeanBinSpec(bin, cls)%fnum(gas) = MeanBinSpec(bin, cls)%fnum(gas) &
                            + ConfigBinSpec(bin, cls, c2)%fnum(gas)
                        MeanBinSpec(bin, cls)%fn(gas) = MeanBinSpec(bin, cls)%fn(gas) &
                            + ConfigBinSpec(bin, cls, c2)%fn(gas)
                        MeanBinSpec(bin, cls)%of(gas) = MeanBinSpec(bin, cls)%of(gas) &
                            + ConfigBinSpec(bin, cls, c2)%of(gas)
                        MeanBinSpec(bin, cls)%ts(gas) = MeanBinSpec(bin, cls)%ts(gas) &
                            + ConfigBinSpec(bin, cls, c2)%ts(gas)
                    end do
                end do
            end do
        end do
    end subroutine LoadConfigEnsemble

    !***************************************************************************
    !> Keep the pass's results: each gas into the slot of its rate in this
    !> configuration, the Ibrom model into the slot of the file rate when this
    !> pass fitted it.
    !***************************************************************************
    subroutine SaveAssessment(c, nbins, ibrom_fitted)
        integer, intent(in) :: c
        integer, intent(in) :: nbins
        logical, intent(in) :: ibrom_fitted
        integer :: gas, k, cls

        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            !> A gas with no data in any of this configuration's periods has
            !> no rate here, and nothing was fitted for it
            if (ConfigGasRate(gas, c) <= 0d0) cycle
            k = GasRateSlot(gas, ConfigGasRate(gas, c))
            if (k == 0) cycle
            RegParR(gas, :, k) = RegPar(gas, :)
            do cls = 1, MaxGasClasses
                SlotCnt(gas, cls, k) = MeanBinSpec(max(nbins / 2, 1), cls)%cnt(gas)
            end do
            SlotFilled(gas, k) = .true.
        end do
        if (ibrom_fitted) then
            k = FileRateSlot(ConfigFileRate(c))
            if (k > 0) then
                UnParR(:, k) = UnPar
                StParR(:, k) = StPar
                FileSlotFilled(k) = .true.
            end if
        end if
    end subroutine SaveAssessment

    logical function SlotUsable(gas, k)
        integer, intent(in) :: gas, k
        logical, external :: GasSlotIsWater
        integer :: cls

        SlotUsable = .false.
        if (k < 1 .or. k > nGasRates(gas)) return
        if (.not. SlotFilled(gas, k)) return
        if (GasSlotIsWater(gas)) then
            do cls = RH10, RH90
                if (RegParR(gas, cls, k)%fc /= error) SlotUsable = .true.
            end do
        else
            do cls = 1, MaxGasClasses
                if (RegParR(gas, cls, k)%fc /= error) SlotUsable = .true.
            end do
        end if
    end function SlotUsable

    !***************************************************************************
    !> Decide which slot each rate slot is corrected with. An unusable one
    !> takes the gas's next faster usable rate - a faster rate resolves the
    !> cut-off better - else the next slower, else none (the analytic method).
    !> A gas never borrows another gas's result. Every substitution is listed.
    !***************************************************************************
    subroutine ResolveAssessmentFallbacks()
        integer :: gas, k, j
        character(64) :: tags(GHGNumVar)

        call SpectralGasNames(tags)
        GasRateUse = 0
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            do k = 1, nGasRates(gas)
                if (SlotUsable(gas, k)) then
                    GasRateUse(gas, k) = k
                    cycle
                end if
                do j = k - 1, 1, -1
                    if (SlotUsable(gas, j)) then
                        GasRateUse(gas, k) = j
                        exit
                    end if
                end do
                if (GasRateUse(gas, k) == 0) then
                    do j = k + 1, nGasRates(gas)
                        if (SlotUsable(gas, j)) then
                            GasRateUse(gas, k) = j
                            exit
                        end if
                    end do
                end if
                if (nGasRates(gas) > 1) then
                    if (GasRateUse(gas, k) > 0) then
                        call LogSay('  ' // trim(tags(gas)) // ' at ' &
                            // trim(RateText(GasRate(gas, k))) &
                            // ' Hz: no usable spectral assessment of its own,' &
                            // ' using the ' // trim(RateText(GasRate(gas, GasRateUse(gas, k)))) &
                            // ' Hz result.')
                    else
                        call LogSay('  ' // trim(tags(gas)) // ' at ' &
                            // trim(RateText(GasRate(gas, k))) &
                            // ' Hz: no usable spectral assessment at any rate.')
                    end if
                end if
            end do
        end do

        FileRateUse = 0
        do k = 1, nFileRates
            if (FileSlotUsable(k)) then
                FileRateUse(k) = k
                cycle
            end if
            do j = k - 1, 1, -1
                if (FileSlotUsable(j)) then
                    FileRateUse(k) = j
                    exit
                end if
            end do
            if (FileRateUse(k) == 0) then
                do j = k + 1, nFileRates
                    if (FileSlotUsable(j)) then
                        FileRateUse(k) = j
                        exit
                    end if
                end do
            end if
        end do
    end subroutine ResolveAssessmentFallbacks

    logical function FileSlotUsable(k)
        integer, intent(in) :: k
        FileSlotUsable = FileSlotFilled(k) .and. UnParR(1, k) /= error &
            .and. StParR(1, k) /= error
    end function FileSlotUsable

    !***************************************************************************
    !> Before correcting one period: its own rate's results into RegPar, UnPar
    !> and StPar, which every correction routine reads.
    !***************************************************************************
    subroutine LoadAssessment(lEx)
        type(ExType), intent(in) :: lEx
        integer :: gas, k, primary
        integer, external :: PrimaryWaterSlot

        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            k = GasRateSlot(gas, EffectiveGasRate(lEx, gas))
            if (k > 0) k = GasRateUse(gas, k)
            if (k > 0) then
                RegPar(gas, :) = RegParR(gas, :, k)
            else
                RegPar(gas, :)%Fn = error
                RegPar(gas, :)%fc = error
                RegPar(gas, :)%f2 = error
            end if
        end do
        !> The project-wide RH exponential is the primary hygrometer's
        primary = PrimaryWaterSlot()
        if (primary >= firstGas) then
            RegPar(dum, dum)%e1 = RegPar(primary, dum)%e1
            RegPar(dum, dum)%e2 = RegPar(primary, dum)%e2
            RegPar(dum, dum)%e3 = RegPar(primary, dum)%e3
        end if
        k = FileRateSlot(lEx%ac_freq)
        if (k > 0) k = FileRateUse(k)
        if (k > 0) then
            UnPar = UnParR(:, k)
            StPar = StParR(:, k)
        else
            UnPar = error
            StPar = error
        end if
    end subroutine LoadAssessment

    !***************************************************************************
    !> Reading an assessment file back: the numbers on one row, after its `=`
    !***************************************************************************
    subroutine ParseNumbers(text, vals, nmax, nv)
        character(*), intent(in) :: text
        integer, intent(in) :: nmax
        real(kind = dbl), intent(out) :: vals(nmax)
        integer, intent(out) :: nv
        integer :: i, j, ios
        character(len(text)) :: t

        t = text
        do i = 1, len(t)
            if (t(i:i) == ',') t(i:i) = ' '
        end do
        nv = 0
        vals = error
        i = 1
        do while (i <= len_trim(t) .and. nv < nmax)
            if (t(i:i) == ' ') then
                i = i + 1
                cycle
            end if
            j = index(t(i:), ' ')
            if (j == 0) then
                j = len_trim(t) + 1
            else
                j = i + j - 1
            end if
            nv = nv + 1
            read(t(i:j - 1), *, iostat = ios) vals(nv)
            if (ios /= 0) then
                nv = nv - 1
                exit
            end if
            i = j
        end do
    end subroutine ParseNumbers

    !***************************************************************************
    !> The `rates=` list on a header or label line; nf = 0 when there is none
    !***************************************************************************
    subroutine ParseRatesToken(line, fr, nf)
        character(*), intent(in) :: line
        real(kind = dbl), intent(out) :: fr(MaxRateSlots)
        integer, intent(out) :: nf
        character(1024) :: txt

        fr = -1d0
        nf = 0
        call SpectralStampToken(line, 'rates=', txt)
        if (len_trim(txt) == 0) return
        call ParseNumbers(txt, fr, MaxRateSlots, nf)
    end subroutine ParseRatesToken

    !> The file column that serves a gas's rate slot k
    integer function FileColumnFor(rate, fr, nf)
        real(kind = dbl), intent(in) :: rate
        real(kind = dbl), intent(in) :: fr(MaxRateSlots)
        integer, intent(in) :: nf
        FileColumnFor = 1
        if (nf <= 0) return
        FileColumnFor = NearestSlot(fr, nf, rate)
        if (FileColumnFor < 1) FileColumnFor = 1
    end function FileColumnFor

    !***************************************************************************
    !> One RH row of a hygrometer's table - Fn fc numerosity per file column -
    !> into every rate slot of that hygrometer
    !***************************************************************************
    subroutine StoreRHRow(g, cls, text, fr, nf)
        integer, intent(in) :: g, cls
        character(*), intent(in) :: text
        real(kind = dbl), intent(in) :: fr(MaxRateSlots)
        integer, intent(in) :: nf
        real(kind = dbl) :: vals(3 * MaxRateSlots)
        integer :: nv, k, j

        if (g < firstGas .or. g > lastGas) return
        call ParseNumbers(text, vals, size(vals), nv)
        do k = 1, max(nGasRates(g), 1)
            j = FileColumnFor(GasRate(g, k), fr, nf)
            if (3 * j - 1 > nv) j = 1
            RegParR(g, cls, k)%Fn = vals(3 * j - 2)
            RegParR(g, cls, k)%fc = vals(3 * j - 1)
            SlotFilled(g, k) = .true.
        end do
    end subroutine StoreRHRow

    !***************************************************************************
    !> A hygrometer's RH exponential - `exp_by_rate=a,b,c/a,b,c` if the file
    !> has one per rate, else the single `e` given - into its rate slots
    !***************************************************************************
    subroutine StoreExp(g, line, e1, e2, e3, fr, nf)
        integer, intent(in) :: g
        character(*), intent(in) :: line
        real(kind = dbl), intent(in) :: e1, e2, e3
        real(kind = dbl), intent(in) :: fr(MaxRateSlots)
        integer, intent(in) :: nf
        character(1024) :: txt
        real(kind = dbl) :: vals(3 * MaxRateSlots)
        integer :: nv, k, j, i

        if (g < firstGas .or. g > lastGas) return
        nv = 0
        call SpectralStampToken(line, 'exp_by_rate=', txt)
        if (len_trim(txt) > 0) then
            do i = 1, len_trim(txt)
                if (txt(i:i) == '/') txt(i:i) = ','
            end do
            call ParseNumbers(txt, vals, size(vals), nv)
        end if
        do k = 1, max(nGasRates(g), 1)
            j = FileColumnFor(GasRate(g, k), fr, nf)
            if (nv >= 3 * j) then
                RegParR(g, dum, k)%e1 = vals(3 * j - 2)
                RegParR(g, dum, k)%e2 = vals(3 * j - 1)
                RegParR(g, dum, k)%e3 = vals(3 * j)
            else
                RegParR(g, dum, k)%e1 = e1
                RegParR(g, dum, k)%e2 = e2
                RegParR(g, dum, k)%e3 = e3
            end if
        end do
    end subroutine StoreExp

    !***************************************************************************
    !> A gas's twelve month rows, as Fn/fc per file column, into its rate
    !> slots. The file is keyed by month and the slots by class, so each
    !> column goes through MonthlyRegParToClasses, which writes RegPar; the
    !> caller's own column-1 mapping is put back afterwards.
    !***************************************************************************
    subroutine StoreMonthBlock(g, monthFn, monthfc, fr, nf)
        integer, intent(in) :: g
        real(kind = dbl), intent(in) :: monthFn(12, MaxRateSlots)
        real(kind = dbl), intent(in) :: monthfc(12, MaxRateSlots)
        real(kind = dbl), intent(in) :: fr(MaxRateSlots)
        integer, intent(in) :: nf
        type(RegParType) :: keep(MaxGasClasses)
        integer :: k, j

        if (g < firstGas .or. g > lastGas) return
        keep = RegPar(g, :)
        do k = 1, max(nGasRates(g), 1)
            j = FileColumnFor(GasRate(g, k), fr, nf)
            call MonthlyRegParToClasses(g, monthFn(:, j), monthfc(:, j))
            RegParR(g, :, k)%Fn = RegPar(g, :)%Fn
            RegParR(g, :, k)%fc = RegPar(g, :)%fc
            SlotFilled(g, k) = .true.
        end do
        RegPar(g, :) = keep
    end subroutine StoreMonthBlock

    !***************************************************************************
    !> Ibrom's c1 c2 rows: one pair per file-rate column, into the file-rate
    !> slots
    !***************************************************************************
    subroutine StoreIbromRows(untext, sttext, fr, nf)
        character(*), intent(in) :: untext, sttext
        real(kind = dbl), intent(in) :: fr(MaxRateSlots)
        integer, intent(in) :: nf
        real(kind = dbl) :: uv(2 * MaxRateSlots), sv(2 * MaxRateSlots)
        integer :: nu, ns, k, j

        call ParseNumbers(untext, uv, size(uv), nu)
        call ParseNumbers(sttext, sv, size(sv), ns)
        do k = 1, max(nFileRates, 1)
            j = FileColumnFor(FileRate(k), fr, nf)
            if (2 * j > min(nu, ns)) j = 1
            UnParR(1, k) = uv(2 * j - 1)
            UnParR(2, k) = uv(2 * j)
            StParR(1, k) = sv(2 * j - 1)
            StParR(2, k) = sv(2 * j)
            FileSlotFilled(k) = .true.
        end do
    end subroutine StoreIbromRows

    !***************************************************************************
    !> Warning(119) and the list of rates it refers to: once, and only when
    !> some gas has more than one rate.
    !***************************************************************************
    subroutine WarnMultiRateAssessment()
        integer :: gas, k
        character(64) :: tags(GHGNumVar)
        character(1024) :: line
        character(16) :: n
        logical :: any_multi

        any_multi = .false.
        do gas = firstGas, lastGas
            if (ConfiguredGas(gas) .and. nGasRates(gas) > 1) any_multi = .true.
        end do
        if (.not. any_multi) return

        call SpectralGasNames(tags)
        call LogSayList('')
        do gas = firstGas, lastGas
            if (.not. ConfiguredGas(gas)) cycle
            if (nGasRates(gas) < 2) cycle
            line = ' Warning(119)> ' // trim(tags(gas)) // ':'
            do k = 1, nGasRates(gas)
                write(n, '(i0)') GasRatePeriods(gas, k)
                if (k > 1) line = trim(line) // ','
                line = trim(line) // ' ' // trim(RateText(GasRate(gas, k))) &
                    // ' Hz (' // trim(n) // ' periods)'
            end do
            call LogSayList(trim(line))
        end do
        call ExceptionHandler(119)
    end subroutine WarnMultiRateAssessment
end module m_sa_rates
