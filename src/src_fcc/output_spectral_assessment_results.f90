!***************************************************************************
! output_spectral_assessment_results.f90
! --------------------------------------
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
! \brief       Write results of spectral assessment on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OutputSpectralAssessmentResults(nbins)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nbins
    !> local variables
    integer, external :: CreateDir
    integer :: i
    integer :: pick
    integer :: goodj
    integer :: cls
    integer :: gas
    integer :: month
    integer :: open_status
    integer :: mkdir_status
    real(kind = dbl), external :: func
    real(kind = dbl), external :: kaimal
    character(128) :: Filename
    character(64) :: sa_tags(GHGNumVar)
    character(64) :: sa_name
    integer :: cosp_slots(1 + MaxNumGases)
    integer :: n_cosp
    integer :: k
    character(PathLen) :: FilePath
    character(PathLen) :: SpecDir
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    logical :: proceed
    include '../src_common/interfaces_1.inc'

    !> Create output directory
    mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')
    SpecDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirSpecAn // slash
    mkdir_status = CreateDir('"' // SpecDir(1:len_trim(SpecDir)) // '"')

    !> Select one frequency vector which does not contain only -9999.
    !> If none, it means all spectra are -9999. Basically, spectral
    !> assessment failed.
    goodj = ierror
    ol: do cls = RH10, RH90
        il: do i = 1, nbins - 1
                if (MeanBinSpec(i, cls)%fn(h2o) /= error) then
                    goodj = cls
                    exit ol
                end if
        end do il
    end do ol

    !> Species names for every block header written below. Hoisted: the
    !> passive-gas spectra file further down is written under its own
    !> condition and is reachable without the assessment block.
    call SpectralGasNames(sa_tags)

    !> The cospectrum column set, built once and walked by both the headers
    !> and the row writers below. They used to be a literal naming w/T and
    !> four gases against a loop over w_ts..w_gas4; on a project with more
    !> gases the rows grew and the header did not, so the extra columns were
    !> unlabelled. A gas slot is its own w_ index, so the slot list serves
    !> both. Water is included: unlike the transfer-function blocks, a
    !> cospectrum with w is meaningful for every gas.
    n_cosp = 1
    cosp_slots(1) = w_ts
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        n_cosp = n_cosp + 1
        cosp_slots(n_cosp) = gas
    end do

    !> SPECTRAL ASSESSMENT
    if (FCCsetup%do_spectral_assessment) then

        if (goodj == ierror) then
            call ExceptionHandler(76)
        else
            write(*,'(a)') ' Writing spectral assessment results on file.. '

            !> Transfer function parameters
            Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) // SA_FilePadding  &
                // Timestamp_FilePadding // TxtExt
            FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Transfer_function_parameters_(TFP)_for_&
                &IIR-shaped_filter_(see_Ibrom_et_al._2007_AFM).'
            write(udf,'(a)') 'fc:_IIR_cut-off_frequency'
            write(udf,'(a)') 'Fn:_normalization_parameter'
            write(udf,'(a)') 'Water_vapour_TFP_are_calculated_for_9_RH_classes.'
            write(udf,'(a)') 'Other_gases_TFP_are_calculated_on_a_monthly_base_&
                &(currently_all_months_together_).'
            write(udf,'(a)') '-----------------------------------------------------&
                &-----------------------------'
            write(udf,'(a)') 'Water vapour TFP              Fn          fc    numerosity'
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class   5 - 15% = ', &
                RegPar(h2o, RH10)%Fn, RegPar(h2o, RH10)%fc, &
                MeanBinSpec(nbins/2, RH10)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  15 - 25% = ', &
                RegPar(h2o, RH20)%Fn, RegPar(h2o, RH20)%fc, &
                MeanBinSpec(nbins/2, RH20)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  25 - 35% = ', &
                RegPar(h2o, RH30)%Fn, RegPar(h2o, RH30)%fc, &
                MeanBinSpec(nbins/2, RH30)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  35 - 45% = ', &
                RegPar(h2o, RH40)%Fn, RegPar(h2o, RH40)%fc, &
                MeanBinSpec(nbins/2, RH40)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  45 - 55% = ', &
                RegPar(h2o, RH50)%Fn, RegPar(h2o, RH50)%fc, &
                MeanBinSpec(nbins/2, RH50)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  55 - 65% = ', &
                RegPar(h2o, RH60)%Fn, RegPar(h2o, RH60)%fc, &
                MeanBinSpec(nbins/2, RH60)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  65 - 75% = ', &
                RegPar(h2o, RH70)%Fn, RegPar(h2o, RH70)%fc, &
                MeanBinSpec(nbins/2, RH70)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  75 - 85% = ', &
                RegPar(h2o, RH80)%Fn, RegPar(h2o, RH80)%fc, &
                MeanBinSpec(nbins/2, RH80)%cnt(h2o)
            write(udf,'(a, 2(f11.5,1x), i13)') 'RH class  85 - 95% = ', &
                RegPar(h2o, RH90)%Fn, RegPar(h2o, RH90)%fc, &
                MeanBinSpec(nbins/2, RH90)%cnt(h2o)
            write(udf,'(a)') ''

            !> One block per configured gas but water, whose cutoffs are the
            !> RH-class table above. The reader consumes these blocks
            !> positionally over the same range, so the two must agree on the
            !> count; the header line names the gas so the file stays readable.
            do gas = firstGas, lastGas
                if (GasSlotIsWater(gas)) cycle
                if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
                sa_name = sa_tags(gas)
                call uppercase(sa_name)
                write(udf,'(a)') trim(sa_name) // '            TFP            &
                    &Fn          fc'

                if (FCCsetup%SA%class(gas, JAN) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'January            = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JAN))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JAN))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'January            = ', error, error
                end if
                if (FCCsetup%SA%class(gas, FEB) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'February           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, FEB))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, FEB))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'February           = ', error, error
                end if
                if (FCCsetup%SA%class(gas, MAR) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'March              = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, MAR))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, MAR))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'March              = ', error, error
                end if
                if (FCCsetup%SA%class(gas, APR) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'April              = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, APR))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, APR))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'April              = ', error, error
                end if
                if (FCCsetup%SA%class(gas, MAY) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'May                = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, MAY))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, MAY))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'May                = ', error, error
                end if
                if (FCCsetup%SA%class(gas, JUN) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'June               = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JUN))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JUN))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'June               = ', error, error
                end if
                if (FCCsetup%SA%class(gas, JUL) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'July               = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, JUL))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, JUL))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'July               = ', error, error
                end if
                if (FCCsetup%SA%class(gas, AUG) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'August             = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, AUG))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, AUG))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'August             = ', error, error
                end if
                if (FCCsetup%SA%class(gas, SEP) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'September          = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, SEP))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, SEP))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'September          = ', error, error
                end if
                if (FCCsetup%SA%class(gas, OCT) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'October            = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, OCT))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, OCT))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'October            = ', error, error
                end if
                if (FCCsetup%SA%class(gas, NOV) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'November           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, NOV))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, NOV))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'November           = ', error, error
                end if
                if (FCCsetup%SA%class(gas, DEC) /= 0) then
                    write(udf,'(a, 2(f11.5,1x))') 'December           = ', &
                        RegPar(gas, FCCsetup%SA%class(gas, DEC))%Fn, &
                        RegPar(gas, FCCsetup%SA%class(gas, DEC))%fc
                else
                    write(udf,'(a, 2(f11.5,1x))') 'December           = ', error, error
                end if
                write(udf, *)
            end do

            !> Exponential fit f_co vs. RH
            write(udf,'(a)') 'RH/fc_exponential_fit_parameters_for_water_vapour&
                &_spectral_corrections'
            write(udf,'(a)') '-----------------------------------'
            write(udf,'(a)') '         exp1         exp2         exp3'
            write(udf,'(3(f13.6))') RegPar(dum, dum)%e1, &
                RegPar(dum, dum)%e2, RegPar(dum, dum)%e3
            write(udf,'(a)') ''
            write(udf,'(a)') ''

            write(udf, '(a)') 'High-pass_correction_factor_model_parameters'
            write(udf, '(a)') 'Model: CF = [c1 * u / (c2 + f_co) + 1] after_&
                &Ibrom_et_al_(2007_AFM)'
            write(udf, '(a)') '---------------------------------------------&
                &----------------------'
            write(udf, '(a)') '                   c1          c2'
            write(udf, '(a, 2(f11.7,1x))') 'unstable = ',UnPar(1), UnPar(2)
            write(udf, '(a, 2(f11.7,1x))') 'stable   = ',StPar(1), StPar(2)
            close(udf)
        end if
    end if

    !> ENSEMBLE AVERAGED SPECTRA
    if (FCCsetup%do_spectral_assessment .or. EddyFlowProj%out_avrg_spec) then
        !> Average H2O spectra, sorted in RH classes, and predicted spectra
        !> (RHS of eq. 6 in Ibrom et al. 2007, AFM)
        if (goodj == ierror) then
            call ExceptionHandler(77)
        else
            !> Initialize file
            Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                // H2OAvrg_FilePadding // Timestamp_FilePadding // CsvExt
            FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Binned_average_and_predicted_H2O_spectra_sorted_by_RH-class.'
            write(udf,'(a)') ',RH=0.1,,,,RH=0.2,,,,RH=0.3,,,,RH=0.4,,,,RH=0.5,,,,RH=0.6&
                &,,,,RH=0.7,,,,RH=0.8,,,,RH=0.9'

            !> Add number of spectra per class
            dataline = ''
            call AddDatum(dataline, '', separator)
            do cls = RH10, RH90
                call WriteDatumInt(MeanBinSpec(1, cls)%cnt(h2o), datum, &
                    EddyFlowProj%err_label)
                call AddDatum(dataline, 'n_=_' // datum(1:len_trim(datum)), &
                    separator)
                call AddDatum(dataline, '', separator)
                call AddDatum(dataline, '', separator)
                call AddDatum(dataline, '', separator)
            end do
            write(udf,'(a)') dataline(1:len_trim(dataline) - 1)
            write(udf,'(a)') 'nat_freq,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)&
                            &,avrg_sp(T),avrg_sp(h2o),denoised_avrg_sp(h2o),pred_sp(h2o)'

            do i = 1, nbins - 1
                call clearstr(dataline)
                if (MeanBinSpec(i, goodj)%fn(h2o) /= error) then
                    call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o), &
                        datum, EddyFlowProj%err_label)
                    call AddDatum(dataline, datum, separator)
                    do cls = RH10, RH90
                        if (MeanBinSpecAvailable(cls, h2o)) then
                            !> Natural frequency
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * MeanBinSpec(i, cls)%ts(h2o), datum, &
                                EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Ensemble averaged spectrum
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * MeanBinSpec(i, cls)%of(h2o), datum, &
                                EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Denoised ensemble averaged spectrum
                            call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(h2o) &
                                * dMeanBinSpec(i, cls)%of(h2o), datum, &
                                EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                            !> Modelled spectrum
                            call WriteDatumFloat(RegPar(h2o, cls)%Fn &
                                * (1d0 / (1d0 + (MeanBinSpec(i, goodj)%fn(h2o) &
                                / RegPar(h2o, cls)%fc)**2 )) &
                                * MeanBinSpec(i, cls)%ts(h2o) &
                                * MeanBinSpec(i, goodj)%fn(h2o), datum, &
                                EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)

                        else
                            call AddDatum(dataline, &
                                trim(adjustl(EddyFlowProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyFlowProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyFlowProj%err_label)), separator)
                            call AddDatum(dataline, &
                                trim(adjustl(EddyFlowProj%err_label)), separator)
                        end if
                    end do
                    write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                end if
            end do
            close(udf)

            !> =====================================================================
            !> Average CO2/CH4/GAS4 spectra, sorted by month, and predicted spectra (RHS of
            !> eq. 6 in Ibrom et al. (2007, AFM)

            !> Select one frequency vector which does not contain only -9999.
            !> If none, it means all spectra are -9999, so -9999 is written.
            !> Basically, run failed.
            call GetFnIndex(MeanBinSpec, &
                size(MeanBinSpec, 1), size(MeanBinSpec, 2), goodj, pick)

            !> Write output file if valid goodj and pick were found
            if (goodj > 0 .and. goodj < MaxGasClasses &
                .and. pick > 0 .and. pick < gas4) then
                Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) // PASGAS_Avrg_FilePadding  &
                    // Timestamp_FilePadding // CsvExt
                FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
                open(udf, file = FilePath, iostat = open_status)
                if (open_status /= 0) call ExceptionHandler(64)
                write(udf,'(a)') 'Binned_average_and_predicted_spectra_for_passive_gases'
                !> Add number of spectra per class
                dataline = ''
                call AddDatum(dataline, '', separator)
                do month = JAN, JAN
                    do gas = firstGas, lastGas
                        if (gas - firstGas + 1 > &
                            min(EddyFlowProj%gas_num, MaxNumGases)) exit
                        if (.not. GasSlotIsWater(gas)) then
                            if (FCCsetup%SA%class(gas, month) /= 0) then
                                call WriteDatumInt(MeanBinSpec(1, FCCsetup%SA%class(gas, month))%cnt(gas) &
                                    , datum, EddyFlowProj%err_label)
                            else
                                call WriteDatumInt(0, datum, EddyFlowProj%err_label)
                            end if
                            call AddDatum(dataline, 'n_=_' // datum(1:len_trim(datum)), separator)
                            call AddDatum(dataline, '', separator)
                            call AddDatum(dataline, '', separator)
                            call AddDatum(dataline, '', separator)
                        end if
                    end do
                end do
                write(udf,'(a)') dataline(1:len_trim(dataline) - 1)

                !> Header generated over the same range as the row writer
                !> below, four columns per gas. It used to be a literal
                !> naming co2, ch4 and the fourth gas, so on a project with
                !> more gases the rows grew and the header did not - the
                !> columns past the third were unlabelled and misread.
                dataline = 'nat_freq'
                do gas = firstGas, lastGas
                    if (gas - firstGas + 1 > &
                        min(EddyFlowProj%gas_num, MaxNumGases)) exit
                    if (GasSlotIsWater(gas)) cycle
                    sa_name = sa_tags(gas)
                    dataline = trim(dataline) // ',avrg_sp(T),avrg_sp(' &
                        // trim(sa_name) // '),denoised_avrg_sp(' &
                        // trim(sa_name) // '),pred_sp(' &
                        // trim(sa_name) // ')'
                end do
                write(udf,'(a)') trim(dataline)

                do i = 1, nbins - 1
                    call clearstr(dataline)
                    if (MeanBinSpec(i, goodj)%fn(pick) /= error) then
                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick), datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        do month = JAN, JAN
                            do gas = firstGas, lastGas
                                if (gas - firstGas + 1 > &
                                    min(EddyFlowProj%gas_num, MaxNumGases)) exit
                                if (GasSlotIsWater(gas)) cycle
                                if (FCCsetup%SA%class(gas, month) /= 0) then
                                    if (MeanBinSpecAvailable(FCCsetup%SA%class(gas, month), gas))then
                                        !> Natural frequency
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * MeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%ts(gas), datum, EddyFlowProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Ensemble averaged spectrum
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * MeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%of(gas), datum, EddyFlowProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Denoised ensemble averaged spectrum
                                        call WriteDatumFloat(MeanBinSpec(i, goodj)%fn(pick) * dMeanBinSpec(i, &
                                            FCCsetup%SA%class(gas, month))%of(gas), datum, EddyFlowProj%err_label)
                                        call AddDatum(dataline, datum, separator)
                                        !> Modelled spectrum
                                        call WriteDatumFloat(RegPar(gas, FCCsetup%SA%class(gas, month))%fn &
                                            * (1d0 / (1d0 + (MeanBinSpec(i, goodj)%fn(pick) &
                                            / RegPar(gas, FCCsetup%SA%class(gas, month))%fc)**2 )) &
                                            * MeanBinSpec(i, FCCsetup%SA%class(gas, month))%ts(gas) &
                                            * MeanBinSpec(i, goodj)%fn(pick), datum, EddyFlowProj%err_label)
                                        call AddDatum(dataline, datum, separator)

                                    else
                                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                    end if
                                else
                                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                                end if
                            end do
                        end do
                        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                    end if
                end do
                close(udf)
            else
                call ExceptionHandler(92)
            end if
        end if
    end if

    !> ENSEMBLE AVERAGED COSPECTRA
    if (EddyFlowProj%out_avrg_cosp) then

        !> =====================================================================
        !> Ensemble cospectra by time of day
        !> Select one frequency vector which does not contain only -9999.
        !> If none, it means all spectra are -9999, so -9999 is written.
        !> Basically, run failed.
        call GetFnIndex(MeanBinCosp, &
            size(MeanBinCosp, 1), size(MeanBinCosp, 2), goodj, pick)


        if (goodj == ierror .or. pick == ierror) then
            call ExceptionHandler(75)
        else
            Filename = trim(adjustl(EddyFlowProj%id)) &
                // Cosp_FilePadding // Timestamp_FilePadding // CsvExt
            FilePath = trim(adjustl(SpecDir)) // trim(adjustl(Filename))
            open(udf, file = FilePath, iostat = open_status)
            if (open_status /= 0) call ExceptionHandler(64)

            write(udf,'(a)') 'Binned_average_cospectra_sorted_by_time_of_day.'
            write(udf,'(a)') ',00:00-02:59,,,,,03:00-5:59,,,,,06:00-08:59&
                &,,,,,09:00-11:59,,,,,12:00-14:59,,,,,15:00-17:59,,,,,&
                &18:00-20:59,,,,,21:00-23:59'

            !> Add number of cospectra per class
            dataline = ''
            call AddDatum(dataline, '', separator)
            do cls = 1, 8
                do k = 1, n_cosp
                    gas = cosp_slots(k)
                    call WriteDatumInt(MeanBinCosp(1, cls)%cnt(gas), &
                        datum, EddyFlowProj%err_label)
                    call AddDatum(dataline, 'n_=_' // trim(adjustl(datum)), &
                        separator)
                end do
            end do
            write(udf,'(a)') dataline(1:len_trim(dataline) - 1)
            !> Add header piece, one group per time-of-day class and one
            !> column per cospectrum slot, walked in the same order as the
            !> rows below.
            dataline = 'nat_freq'
            do cls = 1, 8
                do k = 1, n_cosp
                    if (cosp_slots(k) == w_ts) then
                        sa_name = 'T'
                    else
                        sa_name = sa_tags(cosp_slots(k))
                    end if
                    dataline = trim(dataline) // ',avrg_cosp(w/' &
                        // trim(sa_name) // ')'
                end do
            end do
            write(udf,'(a)') trim(dataline)

            do i = 1, nbins - 1
                call clearstr(dataline)
                if (MeanBinCosp(i, goodj)%fn(pick) /= error) then
                    call WriteDatumFloat(MeanBinCosp(i, goodj)%fn(pick), &
                        datum, EddyFlowProj%err_label)
                    call AddDatum(dataline, datum, separator)
                    do cls = 1, 8
                        do k = 1, n_cosp
                            gas = cosp_slots(k)
                            if (MeanBinCospAvailable(cls, gas))then
                                call WriteDatumFloat(MeanBinCosp(i, goodj)%fn(pick) &
                                    * MeanBinCosp(i, cls)%of(gas), datum, &
                                        EddyFlowProj%err_label)
                                call AddDatum(dataline, datum, separator)
                            else
                                call AddDatum(dataline, &
                                    trim(adjustl(EddyFlowProj%err_label)), &
                                    separator)
                            end if
                        end do
                    end do
                    write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
                end if
            end do
            close(udf)
        end if

        !> =====================================================================
        !> Stability-sorted cospectra and models,
        !> only if at least one fit succeed
        proceed = .false.
        do gas = w_ts, w_co2
            if ((MassPar(gas, unstable)%a0   /= error .or. &
                MassPar(gas, unstable)%fpeak /= error .or. &
                MassPar(gas, unstable)%mu    /= error) .and. &
                (MassPar(gas, stable)%a0     /= error .or. &
                MassPar(gas, stable)%fpeak   /= error .or. &
                MassPar(gas, stable)%mu      /= error)) then
                proceed = .true.
                exit
            end if
        end do

        !> If not one fit went good, alert on output and exit routine
        if (.not. proceed) then
            call ExceptionHandler(45)
            return
        end if

        Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) // Stability_FilePadding  &
            // Timestamp_FilePadding // CsvExt
        FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))

        open(udf, file = FilePath, iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(64)

        write(udf,'(a)') 'Ensemble_cospectra,fitted_Massman_cospectra_and_Kaimal_cospectra.'
        write(udf,'(a)') 'Massman_model:'
        write(udf,'(a)') 'n*Co(ws)/cov(ws)=a0*(fn/fpeak)/(1+(fn/fpeak)^(2mu))^(1.167/mu)'
        write(udf,'(a)') '(note_that_slope_parameter_"m"_is_fixed_to_m=0.75)'
        write(udf,'(a)') 'a0=normalization_factor'
        write(udf,'(a)') 'fpeak=frequency_at_which_cospectrum_attains_highest_value'
        write(udf,'(a)') 'mu=broadness_factor'
        write(udf,'(a)') ''
        write(udf,'(a)') 'For_Kaimal_ideal_cospectra:see_e.g._Moncrieff_et_al._(1997_JoH)_Eqs.12-16.'
        write(udf,'(a)') '(for_stable_stratifications_two_extreme_ideal_cospectra_are_reported_&
            &corresponding_to_z/L=0.01_(slightly_stable)_and_to_z/L=10.0_(very_stable).'
        write(udf,'(a)') '-----------------------------------------------------------------------'
        write(udf,'(a)') 'Massman_model_fit_parameters_for_this_run:'
        !> One pair of blocks per cospectrum slot. These were five
        !> hand-written pairs naming CO2, H2O, CH4 and whatever g4lab had
        !> parsed out of the FLUXNET header, so a project with more gases
        !> reported fit parameters for four of them and no way to tell which.
        do k = 1, n_cosp
            gas = cosp_slots(k)
            if (gas == w_ts) then
                sa_name = 'T'
            else
                sa_name = sa_tags(gas)
                call uppercase(sa_name)
            end if
            write(udf,'(a)') 'w/' // trim(sa_name) // ' (unstable):'
            write(udf,'(a, f10.4)') 'a0,', MassPar(gas, unstable)%a0
            write(udf,'(a, f10.4)') 'fpeak,', MassPar(gas, unstable)%fpeak
            write(udf,'(a, f10.4)') 'mu,', MassPar(gas, unstable)%mu
            write(udf,'(a)') 'w/' // trim(sa_name) // ' (stable):'
            write(udf,'(a, f10.4)') 'a0,', MassPar(gas, stable)%a0
            write(udf,'(a, f10.4)') 'fpeak,', MassPar(gas, stable)%fpeak
            write(udf,'(a, f10.4)') 'mu,', MassPar(gas, stable)%mu
        end do
        write(udf,'(a, f10.4)') '-----------------------------------------------------------------------'


        write(udf,'(a)') 'unstable_(-650<L<0),,,,,,,,,,,,,,,,,,,,,,,,,stable(0<L<1000)'
        !> Add header piece: the unstable groups, then the stable ones,
        !> each walking the same slot list as the rows below.
        dataline = ''
        do k = 1, n_cosp
            gas = cosp_slots(k)
            if (gas == w_ts) then
                sa_name = 'T'
            else
                sa_name = sa_tags(gas)
            end if
            dataline = trim(dataline) // 'fn,avrg_cosp(w/' // trim(sa_name) &
                // '),fit_cosp(w/' // trim(sa_name) // '),kaimal_cosp,,'
        end do
        do k = 1, n_cosp
            gas = cosp_slots(k)
            if (gas == w_ts) then
                sa_name = 'T'
            else
                sa_name = sa_tags(gas)
            end if
            dataline = trim(dataline) // 'fn,avrg_cosp(w/' // trim(sa_name) &
                // '),fit_cosp(w/' // trim(sa_name) &
                // '),kaimal_cosp_zL_0.01,kaimal_cosp_zL_10.0,,'
        end do
        write(udf,'(a)') dataline(1:len_trim(dataline) - 2)


        do i = 1, ndkf
            call clearstr(dataline)
            !> Unstable
            do k = 1, n_cosp
                gas = cosp_slots(k)
                if (MeanStabCospAvailable(unstable, gas))then
                    if (MeanStabilityCosp(i, unstable)%fn(gas) /= error .and. MeanStabilityCosp(i, unstable)%fn(gas) /= 0d0) then
                        call WriteDatumFloat(MeanStabilityCosp(i, unstable)%fn(gas), datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ensemble cospectrum
                        if (MeanStabilityCosp(i, unstable)%cnt(gas) > 20) then
                            call WriteDatumFloat(MeanStabilityCosp(i, unstable)%of(gas), datum, EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        else
                            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        end if
                        !> Fitted model cospectrum
                        call WriteDatumFloat(dexp(func(MeanStabilityCosp(i, unstable)%fn(gas), &
                            MassPar(gas, unstable)%a0, MassPar(gas, unstable)%fpeak, &
                            MassPar(gas, unstable)%mu)), datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ideal cospectrum
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, unstable)%fn(gas), -9999d0, 'unstable'), &
                            datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call AddDatum(dataline, '', separator)
                    else
                        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, '', separator)
                    end if
                else
                    call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, '', separator)
                end if
            end do

            !> Stable
            do k = 1, n_cosp
                gas = cosp_slots(k)
                if (MeanStabCospAvailable(stable, gas))then
                    if (MeanStabilityCosp(i, stable)%fn(gas) /= error .and. MeanStabilityCosp(i, stable)%fn(gas) /= 0d0) then
                            call WriteDatumFloat(MeanStabilityCosp(i, stable)%fn(gas), datum, EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        !> Ensemble cospectrum
                        if (MeanStabilityCosp(i, stable)%cnt(gas) > 10) then
                            call WriteDatumFloat(MeanStabilityCosp(i, stable)%of(gas), datum, EddyFlowProj%err_label)
                            call AddDatum(dataline, datum, separator)
                        else
                            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
                        end if
                        !> Fitted model cospectrum
                        call WriteDatumFloat(dexp(func(MeanStabilityCosp(i, stable)%fn(gas), &
                            MassPar(gas, stable)%a0, MassPar(gas, stable)%fpeak, &
                            MassPar(gas, stable)%mu)), datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        !> Ideal cospectrum
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, stable)%fn(gas), 1d-2, 'stable'), &
                            datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call WriteDatumFloat(kaimal(MeanStabilityCosp(i, stable)%fn(gas), 1d1 , 'stable'), &
                            datum, EddyFlowProj%err_label)
                        call AddDatum(dataline, datum, separator)
                        call AddDatum(dataline, '', separator)
                    else
                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                        call AddDatum(dataline, '', separator)
                    end if
                else
                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, trim(adjustl(EddyFlowProj%err_label)), separator)
                    call AddDatum(dataline, '', separator)
                end if
            end do
            write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
        end do
        close(udf)
    end if
    write(*,'(a)') ' Done.'

end subroutine OutputSpectralAssessmentResults

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GetFnIndex(LocSpec, nrow, ncol, goodj, pick)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    type (MeanSpectraType), intent(in) :: LocSpec(nrow, ncol)
    integer, intent(out) :: goodj
    integer, intent(out) :: pick
    !> local variables
    integer :: i
    integer :: cls


    goodj = ierror
    pick = ierror

    ol: do cls = 1, ncol
        il: do i = 1, nrow - 1
                if (LocSpec(i, cls)%fn(ts) > 0d0) then
                    goodj = cls
                    pick = ts
                    exit ol
                end if
        end do il
    end do ol
    if (goodj == ierror) then
        ol1: do cls = 1, ncol
            il1: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(co2) > 0d0) then
                        goodj = cls
                        pick = co2
                        exit ol1
                    end if
            end do il1
        end do ol1
    end if
    if (goodj == ierror) then
        ol2: do cls = 1, ncol
            il2: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(ch4) > 0d0) then
                        goodj = cls
                        pick = ch4
                        exit ol2
                    end if
            end do il2
        end do ol2
    end if
    if (goodj == ierror) then
        ol3: do cls = 1, ncol
            il3: do i = 1, nrow - 1
                    if (LocSpec(i, cls)%fn(gas4) > 0d0) then
                        goodj = cls
                        pick = gas4
                        exit ol3
                    end if
            end do il3
        end do ol3
    end if
end subroutine GetFnIndex
