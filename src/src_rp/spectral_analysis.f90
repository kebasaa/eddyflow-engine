!***************************************************************************
! spectral_analysis.f90
! ---------------------
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
! \brief       Calculate spectra and co-spectra \n
!              and write them on output file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SpectralAnalysis(date, time, bf, Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    !> local variables
    integer :: i
    integer :: var
    integer :: bcnt(Meth%spec%nbins)
    real(kind = dbl) :: bnf(Meth%spec%nbins)
    real(kind = dbl) :: AuxSet(N, M)
    !> Tapering and FourierTransform both work in place, and Set is
    !> intent(in). The binned pass therefore transforms a copy - it used to
    !> transform Set itself, which the compiler only allowed because both have
    !> implicit interfaces, and which left the caller's array destroyed.
    real(kind = dbl) :: WorkSet(N, M)
    real(kind = dbl) :: nf(N/2)
    real(kind = dbl) :: sumw, auxsumw
    character(13) :: Datestring
    type (SpectralType) :: Spectrum(N/2 + 1)
    type (SpectralType) :: AuxSpectrum(N/2 + 1)
    type (SpectralType) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType) :: Cospectrum(N/2 + 1)
    type (SpectralType) :: AuxCospectrum(N/2 + 1)
    type (SpectralType) :: BinnedCospectrum(Meth%spec%nbins)
    type (SpectralType) :: Ogive(N/2 + 1)
    type (SpectralType) :: CoOgive(N/2 + 1)
    type (SpectralType) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType) :: BinnedCoOgive(Meth%spec%nbins)
    logical :: DoSpectrum(GHGNumVar)
    logical :: DoCospectrum(GHGNumVar)
    logical :: proceed

    !> Determine which spectra and cospectra can be calculated
    call DetectFeasibleSpectraAndCospectra(DoSpectrum, DoCospectrum)

    !> Determine if at least one (co)spectrum has to be written on output
    proceed = .false.
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .or. RPsetup%out_full_cosp(var)) then
            proceed = .true.
            exit
        end if
    end do

    !> If binned (co)spectra or at least one full (co)spectrum are requested,
    !> perform all related calculations
    call LogSay('  Calculating (co)spectra..')
    Datestring = date(1:4) // date(6:7) // date(9:10) &
               // '-' // time(1:2) // time(4:5)

    !> Define the "natural frequency" vector nf, extending
    !> from f_min = 1/N --> f_max = AcFreq/2 (= Nyquist frequency)
    !> it is period-dependent
    do i = 1, N/2
        nf(i) = dble(i) * Metadata%ac_freq / dble(N)
    end do

    !> Use "Aux" variables to calculate degraded covariances and to
    !> output full co-spectrum wT. "Aux" variables are not tapered
    AuxSet = Set
    !> (do not) taper dataset (squared window makes no tapering), but needed to
    !> calculate "window squared and summed"
    call Tapering('squared', AuxSet, N, M, auxsumw)

    !> Fft and calculate cospectra
    call FourierTransform(AuxSet, N, M)

    call AllCospectra(AuxSet, auxsumw, AuxSpectrum, AuxCospectrum, &
        DoSpectrum, DoCospectrum, N, M)

    !> Rebuild whatever runs slower than the file, from its own samples rather
    !> than from the interpolation. Unnormalised, like everything on this path.
    !>
    !> Before the ogives, not after: an ogive is the integral of the spectrum
    !> beside it, so rebuilding afterwards would leave the two files disagreeing
    !> about the same column.
    call SlowColumnSpectra(Set, N, M, 'squared', .false., &
        DoSpectrum, DoCospectrum, AuxSpectrum, AuxCospectrum)

    !> Calculate ogives if requested
    if (RPsetup%out_bin_og) &
        call AllOgives(AuxSpectrum, AuxCospectrum, DoSpectrum, DoCospectrum, Ogive, CoOgive, N)

    !> Nothing a column could not have measured reaches the file. Before the
    !> write and before the degraded covariances, which read w_ts - the sonic,
    !> which is never capped.
    call CapSpectraAtColumnNyquist(nf, N/2, AuxSpectrum, AuxCospectrum)

    !> Write out full co-spectra as requested
    if (proceed) &
        call WriteOutFullCoSpectra(Datestring, nf, &
            AuxSpectrum, AuxCospectrum, DoSpectrum, DoCospectrum, N)

    !> Output degraded covariances calculated by integrating filtered cospectra
    if (DoCospectrum(w_ts)) then
        call DegradedCovariances(nf, AuxCospectrum, N)
    else
        Essentials%degH(:) = error
    end if

    !> If binned (co-)spectra are requested, perform all related calculations
    if (RPsetup%out_bin_sp .or. RPsetup%out_bin_og) then
        !> From here, tapering is applied.
        !> Taper dataset
        WorkSet = Set
        call Tapering(RPsetup%tap_win, WorkSet, N, M, sumw)
        !> Fft and calculate cospectra
        call FourierTransform(WorkSet, N, M)
        call AllCospectra(WorkSet, sumw, Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N, M)
        !> Normalization by variances and covariances
        call NormalizeCoSpectra(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N)
        !> After the normalisation, not before: a rebuilt column is normalised
        !> by the variance of its OWN series inside the routine, and passing it
        !> through NormalizeCoSpectra as well would divide it twice by two
        !> different numbers.
        call SlowColumnSpectra(Set, N, M, RPsetup%tap_win, .true., &
            DoSpectrum, DoCospectrum, Spectrum, Cospectrum)
        !> Binned cospectra session
        if (RPsetup%out_bin_sp) then
            !> Exponential binning of frequencies, spectra and co-spectra
            call ExpAvrgCospectra(bf, nf, Spectrum, Cospectrum, N, bnf &
                , BinnedSpectrum, BinnedCospectrum, bcnt)
            call CapSpectraAtColumnNyquist(bnf, Meth%spec%nbins, &
                BinnedSpectrum, BinnedCospectrum)
            !> Write co-spectra on output file in csv format
            call WriteOutBinnedCoSpectra(Datestring, bnf, bcnt, BinnedSpectrum, BinnedCospectrum &
                , DoSpectrum, DoCospectrum)
        end if
        !> Ogive session
        if (RPsetup%out_bin_og) then
            !> Exponential binning of ogives
            call ExpAvrgOgives(bf, nf, Ogive, CoOgive, N, bnf &
                , BinnedOgive, BinnedCoOgive, bcnt)
            !> Write co-ogives on output file in csv format
            call WriteOutBinnedOgives(Datestring, bnf, bcnt, BinnedOgive, BinnedCoOgive &
                , DoSpectrum, DoCospectrum)
        end if
    end if
    call LogSay('  Done.')
end subroutine SpectralAnalysis

!***************************************************************************
!
! \brief       Blank whatever a column could not have measured.
! \author      Jonathan Muller, ETH Zurich
! \note        A column sampled below the file's rate resolves nothing above
!              its OWN Nyquist. FixDatasetForSpectra interpolates its missing
!              rows up onto the fast grid before the transform, so without this
!              every bin up to the station's Nyquist carries a number - and
!              every one above the column's own Nyquist is an artefact of that
!              interpolation rather than a measurement.
!
!              Blanked to the error code and not to zero: a spectral density of
!              zero is a claim about the data, and this is the absence of one.
!              A column at the file's own rate is never touched, so a
!              single-rate site is unaffected.
! \sa          column_sampling.f90
!***************************************************************************
subroutine CapSpectraAtColumnNyquist(freq, nfreq, Spectrum, Cospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nfreq
    real(kind = dbl), intent(in) :: freq(nfreq)
    type (SpectralType), intent(inout) :: Spectrum(nfreq)
    type (SpectralType), intent(inout) :: Cospectrum(nfreq)
    !> local variables
    integer :: i
    integer :: j
    real(kind = dbl) :: nyquist
    real(kind = dbl), external :: ColumnAcFreq

    do j = u, GHGNumVar
        if (ColumnAcFreq(j) >= Metadata%ac_freq) cycle
        nyquist = ColumnAcFreq(j) / 2d0
        do i = 1, nfreq
            !> A binned axis carries the error code for a bin no natural
            !> frequency fell into, which is already the answer this gives.
            if (freq(i) == error) cycle
            if (freq(i) > nyquist) then
                Spectrum(i)%of(j) = error
                !> The cospectral index space is the same space as u..lastGas.
                Cospectrum(i)%of(j) = error
            end if
        end do
    end do
end subroutine CapSpectraAtColumnNyquist

!***************************************************************************
!
! \brief       Rebuild a slower column's (co)spectrum on its own grid.
! \author      Jonathan Muller, ETH Zurich
! \note        The full-rate pass transforms every column on the file's row
!              grid, where a slower column has been interpolated up from one
!              real sample in `stride`. Below its own Nyquist that spectrum is
!              the column's own convolved with the interpolation kernel -
!              attenuated, and increasingly so towards the cut-off. This takes
!              the real samples instead.
!
!              Decimating by an integer factor over the same span leaves the
!              natural-frequency grid ALONE: nf(i) = i * F / N becomes
!              i * (F/k) / (N/k), the same number. So the result drops into the
!              same arrays at the same indices, and the shared axis, the
!              binning and the writers need no notion of a second rate. The
!              two grids coincide exactly when k divides N and to within
!              (k-1)/N - a twentieth of a per cent on a half-hour at 10 Hz -
!              when it does not.
!
!              The cospectrum is a property of the pair, not of the gas: `w` is
!              re-sampled onto the gas's grid, at the same instants, and how it
!              is taken there follows the instrument. Point-sampled by default,
!              because a point-sampled gas paired against an averaged `w`
!              biases the covariance; block-averaged when the instrument says
!              it integrates over its interval, which is the same thing the
!              instrument did to the gas.
! \sa          column_sampling.f90, fix_dataset_for_spectra.f90
!***************************************************************************
subroutine SlowColumnSpectra(Set, N, M, tap_win, normalise, &
    DoSpectrum, DoCospectrum, Spectrum, Cospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    character(*), intent(in) :: tap_win
    logical, intent(in) :: normalise
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    real(kind = dbl), intent(in) :: Set(N, M)
    type (SpectralType), intent(inout) :: Spectrum(N/2 + 1)
    type (SpectralType), intent(inout) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: i
    integer :: j
    integer :: row
    integer :: lo
    integer :: stride
    integer :: nd
    integer :: phase
    real(kind = dbl) :: freq
    real(kind = dbl) :: sumw
    real(kind = dbl) :: var_gas
    real(kind = dbl) :: cov_wgas
    real(kind = dbl) :: mean_gas
    real(kind = dbl) :: mean_w
    real(kind = dbl), allocatable :: DecSet(:, :)
    real(kind = dbl), allocatable :: raw_w(:)
    real(kind = dbl), allocatable :: raw_gas(:)
    real(kind = dbl), allocatable :: dspec(:)
    real(kind = dbl), allocatable :: dcosp(:)
    real(kind = dbl), external :: ColumnAcFreq

    do j = u, GHGNumVar
        if (j == w) cycle
        if (.not. DoSpectrum(j) .and. .not. DoCospectrum(j)) cycle

        freq = ColumnAcFreq(j)
        if (freq >= Metadata%ac_freq) cycle

        stride = nint(Metadata%ac_freq / freq)
        if (stride <= 1) cycle

        !> The column's samples sit at a fixed offset inside each interval;
        !> SpecPhase records it in E2Primes rows and SpecRowOffset is where
        !> this set starts in those rows. Sampling at the wrong offset would
        !> read interpolated blends of two real samples and quietly attenuate
        !> the very thing this exists to measure.
        phase = modulo(SpecPhase(j) - SpecRowOffset, stride)

        !> Counted from the phase, so the last sample lands on or before row N
        !> and there is nothing to clamp. Clamping would have repeated the last
        !> row - a fabricated sample, in the one routine whose whole purpose is
        !> to use only real ones.
        !>
        !> Even, because OneSidedPowerSpectrum walks the transform in pairs.
        nd = (N - phase) / stride
        nd = nd - mod(nd, 2)
        !> Too few samples to say anything about a spectrum. Left as the
        !> full-rate pass produced it, and capped at this column's Nyquist by
        !> CapSpectraAtColumnNyquist either way.
        if (nd < 16) cycle

        allocate(raw_w(nd), raw_gas(nd), DecSet(nd, 2))
        allocate(dspec(nd/2 + 1), dcosp(nd/2 + 1))

        do i = 1, nd
            row = phase + 1 + (i - 1) * stride
            raw_gas(i) = Set(row, j)
            if (E2Col(j)%instr%integrates) then
                !> Averaged over the interval the sample closes, which is what
                !> the instrument itself did to the gas.
                lo = max(1, row - stride + 1)
                raw_w(i) = sum(Set(lo:row, w)) / dble(row - lo + 1)
            else
                raw_w(i) = Set(row, w)
            end if
        end do

        !> Normalised by what THIS series varies by, not by the full-rate
        !> statistics: the two differ, and dividing by the wrong one would
        !> rescale the whole curve.
        mean_gas = sum(raw_gas) / dble(nd)
        mean_w = sum(raw_w) / dble(nd)
        var_gas = sum((raw_gas - mean_gas)**2) / dble(nd)
        cov_wgas = sum((raw_w - mean_w) * (raw_gas - mean_gas)) / dble(nd)

        DecSet(1:nd, 1) = raw_w(1:nd)
        DecSet(1:nd, 2) = raw_gas(1:nd)
        call Tapering(tap_win, DecSet, nd, 2, sumw)
        call FourierTransform(DecSet, nd, 2)

        if (DoSpectrum(j)) then
            call OneSidedPowerSpectrum(DecSet(:, 2), DecSet(:, 2), &
                freq, sumw, dspec, nd)
            if (normalise .and. var_gas /= 0d0) dspec = dspec / var_gas
            Spectrum(1:nd/2 + 1)%of(j) = dspec(1:nd/2 + 1)
            Spectrum(nd/2 + 2:N/2 + 1)%of(j) = error
        end if

        if (DoCospectrum(j)) then
            call OneSidedPowerSpectrum(DecSet(:, 1), DecSet(:, 2), &
                freq, sumw, dcosp, nd)
            if (normalise .and. cov_wgas /= 0d0) dcosp = dcosp / cov_wgas
            !> The cospectral index space is the same space as u..lastGas.
            Cospectrum(1:nd/2 + 1)%of(j) = dcosp(1:nd/2 + 1)
            Cospectrum(nd/2 + 2:N/2 + 1)%of(j) = error
        end if

        deallocate(raw_w, raw_gas, DecSet, dspec, dcosp)
    end do
end subroutine SlowColumnSpectra



!***************************************************************************
!
! \brief       Determine which spectra and cospectra can be evaluated, given \n
!              available variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DetectFeasibleSpectraAndCospectra(DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(inout) :: DoSpectrum(GHGNumVar)
    logical, intent(inout) :: DoCospectrum(GHGNumVar)


    !> Decision based on variables availability
    !> Spectra
    DoSpectrum = .false.
    where (SpecCol(u:GHGNumVar)%present)
        DoSpectrum(u:GHGNumVar) = .true.
    end where

    !> Cospectra
    DoCospectrum = .false.
    if (SpecCol(w)%present) then
        where (SpecCol(u:GHGNumVar)%present)
            DoCospectrum(u:GHGNumVar) = .true.
        end where
    end if
end subroutine DetectFeasibleSpectraAndCospectra

!***************************************************************************
!
! \brief       Calculate cospectral densities, starting from un-normalised Fourier \n
!              coefficients, as calculated by rfftf.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AllCospectra(Set, sumw, Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: sumw
    real(kind = dbl), intent(in) :: Set(N, M)
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(out) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(out) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: j
    real(kind = dbl) :: xx(N)
    real(kind = dbl) :: yy(N)

    call LogSayNoAdv('   Cospectral densities..')

    !> Both are intent(out), so a slot the loops below skip is left undefined.
    !> The binning downstream reads the whole range, so "skipped" has to be a
    !> value and not whatever the stack held - error, meaning not performed.
    !> One subscript per part reference: `Spectrum(:)%of(:)` would be two
    !> nonzero-rank part references, which Fortran does not allow.
    do j = u, GHGNumVar
        Spectrum(1:N/2 + 1)%of(j) = error
        Cospectrum(1:N/2 + 1)%of(j) = error
    end do

    !> spectra
    do j = u, GHGNumVar
        if (DoSpectrum(j)) then
            xx(1:N) = Set(1:N, j)
            call OneSidedPowerSpectrum(xx, xx, Metadata%ac_freq, sumw, Spectrum%of(j), N)
        end if
    end do

    !> Cospectra
    do j = w_u, w_lastGas
        if (j /= w_w) then
            if (DoCospectrum(j)) then
                xx(1:N) = Set(1:N, w)
                yy(1:N) = Set(1:N, j)
                call OneSidedPowerSpectrum(xx, yy, Metadata%ac_freq, sumw, Cospectrum%of(j), N)
            end if
        end if
    end do
    call LogSay(' Done.')
end subroutine AllCospectra

!***************************************************************************
!
! \brief       Calculate cospectral densities, starting from un-normalised Fourier \n
!              coefficients, as calculated by rfftf.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AllOgives(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, Ogive, CoOgive, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(in) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(in) :: Cospectrum(N/2 + 1)
    Type(SpectralType), intent(out) :: Ogive(N/2 + 1)
    Type(SpectralType), intent(out) :: CoOgive(N/2 + 1)
    !> local variables
    integer :: i
    integer :: j
    real(kind = dbl)  :: df

    call LogSayNoAdv('   Ogives..')

    !> Same reason as in AllCospectra: intent(out), and the binning reads the
    !> whole range, so a slot no loop fills must say "not performed".
    do j = u, GHGNumVar
        Ogive(1:N/2 + 1)%of(j) = error
        CoOgive(1:N/2 + 1)%of(j) = error
    end do

    df = Metadata%ac_freq / dble(N)

    !> Ogive
    do j = u, GHGNumVar
        if (DoSpectrum(j)) then
            Ogive(N/2+1)%Of(j) = Spectrum(N/2+1)%of(j)
            do i = N/2, 1, -1
                Ogive(i)%of(j) = Ogive(i+1)%of(j) + Spectrum(i)%of(j) * df
            end do
        end if
    end do

    !> CoOgive
    do j = w_u, w_lastGas
        if (j /= w_w) then
            if (DoCospectrum(j)) then
                CoOgive(N/2+1)%Of(j) = Cospectrum(N/2+1)%of(j)
                do i = N/2, 1, -1
                    CoOgive(i)%Of(j) = CoOgive(i+1)%of(j) + Cospectrum(i)%of(j)  * df
                end do
            end if
        end if
    end do

    !> Normalize
    do j = u, GHGNumVar
        if (DoSpectrum(j) .and. Stats%Cov(j, j) /= 0d0 .and. Stats%Cov(j, j) /= error) &
            Ogive(1:N/2 + 1)%of(j) = Ogive(1:N/2 + 1)%of(j) / Stats%Cov(j, j)
    end do

    do j = w_u, w_lastGas
    if (DoCospectrum(j) .and. Stats%Cov(w, j) /= 0d0 .and. Stats%Cov(w, j) /= error) &
        CoOgive(1:N/2 + 1)%of(j) = CoOgive(1:N/2 + 1)%of(j) / Stats%Cov(w, j)
    end do

    call LogSay(' Done.')
end subroutine AllOgives

!***************************************************************************
!
! \brief       Normalize calculated (co)spectra with the relevant (co)variance
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NormalizeCoSpectra(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(inout) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(inout) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: j

    do j = u, GHGNumVar
        if (DoSpectrum(j) .and. Stats%Cov(j, j) /= 0d0 .and. Stats%Cov(j, j) /= error) &
            Spectrum(1:N/2 + 1)%of(j) = Spectrum(1:N/2 + 1)%of(j) / Stats%Cov(j, j)
    end do

    do j = w_u, w_lastGas
    if (DoCospectrum(j) .and. Stats%Cov(w, j) /= 0d0 .and. Stats%Cov(w, j) /= error) &
        Cospectrum(1:N/2 + 1)%of(j) = Cospectrum(1:N/2 + 1)%of(j) / Stats%Cov(w, j)
    end do
end subroutine NormalizeCoSpectra

!***************************************************************************
!
! \brief       Average (co)spectra in the binned frequency range
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExpAvrgCospectra(bf, nf, Spectrum, Cospectrum, N, bin_nf, &
    BinnedSpectrum, BinnedCospectrum, bin_cnt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    real(kind = dbl), intent(in) :: nf(N/2)
    integer, intent(out) :: bin_cnt(Meth%spec%nbins)
    real(kind = dbl), intent(out) :: bin_nf(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedCospectrum(Meth%spec%nbins)
    type (SpectralType), intent(inout) :: Spectrum(N/2 + 1)
    type (SpectralType), intent(inout) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: i = 0
    integer :: j = 0


    !> Align co-spectra to frequencies before averaging
    do i = 1, N/2
        Spectrum(i) = Spectrum(i + 1)
        Cospectrum(i) = Cospectrum(i + 1)
    end do
    Spectrum(N/2 + 1)%of(:) = 0d0
    Cospectrum(N/2 + 1)%of(:) = 0d0

    call LogSayNoAdv('   Binning spectra and cospectra..')
    !> average variables in the exp-spaced ranges

    do i = 1, Meth%spec%nbins
        bin_cnt(i) = 0
        bin_nf(i) = 0.d0
        BinnedSpectrum(i)%of(u:GHGNumVar) = 0d0
        BinnedCospectrum(i)%of(w_u:w_lastGas) = 0d0
        do j = 1, N/2
            if(nf(j) > bf(i) .and. nf(j) <= bf(i + 1)) then
                bin_nf(i) = bin_nf(i) + nf(j)
                BinnedSpectrum(i)%of(u:GHGNumVar) = &
                    BinnedSpectrum(i)%of(u:GHGNumVar) + Spectrum(j)%of(u:GHGNumVar)
                BinnedCospectrum(i)%of(w_u:w_lastGas) = &
                    BinnedCospectrum(i)%of(w_u:w_lastGas) + Cospectrum(j)%of(w_u:w_lastGas)
                bin_cnt(i) = bin_cnt(i) + 1
            end if
        end do
        if(bin_cnt(i) /= 0) then
            bin_nf(i)    = bin_nf(i) / dble(bin_cnt(i))
            BinnedSpectrum(i)%of(u:GHGNumVar) = &
                BinnedSpectrum(i)%of(u:GHGNumVar) / dble(bin_cnt(i))
            BinnedCospectrum(i)%of(w_u:w_lastGas) = &
                BinnedCospectrum(i)%of(w_u:w_lastGas) / dble(bin_cnt(i))
        else
            bin_nf(i)    = error
            BinnedSpectrum(i)%of(u:GHGNumVar) = error
            BinnedCospectrum(i)%of(w_u:w_lastGas) = error
        end if
    end do
    call LogSay(' Done.')
end subroutine ExpAvrgCospectra

!***************************************************************************
!
! \brief       Average ogives in the binned frequency range
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExpAvrgOgives(bf, nf, Ogive, CoOgive, N, bin_nf, &
    BinnedOgive, BinnedCoOgive, bin_cnt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    real(kind = dbl), intent(in) :: nf(N/2)
    integer, intent(out) :: bin_cnt(Meth%spec%nbins)
    real(kind = dbl), intent(out) :: bin_nf(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedCoOgive(Meth%spec%nbins)
    type (SpectralType), intent(inout) :: Ogive(N/2 + 1)
    type (SpectralType), intent(inout) :: CoOgive(N/2 + 1)
    !> local variables
    integer :: i = 0
    integer :: j = 0

    !> Align co-ogives to frequencies before averaging
    do i = 1, N/2
        Ogive(i) = Ogive(i + 1)
        CoOgive(i) = CoOgive(i + 1)
    end do
    Ogive(N/2 + 1)%of(:) = 0d0
    CoOgive(N/2 + 1)%of(:) = 0d0

    call LogSayNoAdv('   Binning ogives..')
    !> average variables in the exp-spaced ranges
    do i = 1, Meth%spec%nbins
        bin_cnt(i) = 0
        bin_nf(i) = 0.d0
        BinnedOgive(i)%of(u:GHGNumVar) = 0d0
        BinnedCoOgive(i)%of(w_u:w_lastGas) = 0d0
        do j = 1, N/2
            if(nf(j) >= bf(i) .and. nf(j) < bf(i + 1)) then
                bin_nf(i) = bin_nf(i) + nf(j)
                BinnedOgive(i)%of(u:GHGNumVar) = &
                    BinnedOgive(i)%of(u:GHGNumVar) + Ogive(j)%of(u:GHGNumVar)
                BinnedCoOgive(i)%of(w_u:w_lastGas) = &
                    BinnedCoOgive(i)%of(w_u:w_lastGas) + CoOgive(j)%of(w_u:w_lastGas)
                bin_cnt(i) = bin_cnt(i) + 1
            end if
        end do
        if(bin_cnt(i) /= 0) then
            bin_nf(i)    = bin_nf(i)    / dble(bin_cnt(i))
            BinnedOgive(i)%of(u:GHGNumVar) = &
                BinnedOgive(i)%of(u:GHGNumVar) / dble(bin_cnt(i))
            BinnedCoOgive(i)%of(w_u:w_lastGas) = &
                BinnedCoOgive(i)%of(w_u:w_lastGas) / dble(bin_cnt(i))
        else
            bin_nf(i)    = error
            BinnedOgive(i)%of(u:GHGNumVar) = error
            BinnedCoOgive(i)%of(w_u:w_lastGas) = error
        end if
    end do
    call LogSay(' Done.')
end subroutine ExpAvrgOgives

!***************************************************************************
!
! \brief       Write binned (co)spectra on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBinnedCoSpectra(String, bnf, bcnt, BinnedSpectrum, BinnedCospectrum &
    , DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bcnt(Meth%spec%nbins)
    real(kind = dbl), intent(in) :: bnf(Meth%spec%nbins)
    character(*), intent(in) :: String
    type (SpectralType), intent(in) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType), intent(in) :: BinnedCospectrum(Meth%spec%nbins)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    !> local variables
    integer :: i
    integer :: j
    integer :: k
    integer :: n_spec, n_cosp
    integer :: spec_slots(E2NumVar), cosp_slots(E2NumVar)
    character(64) :: e2sg(E2NumVar)
    character(PathLen) :: BinCospectraPath
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''
    include '../src_common/interfaces.inc'

    !> Column names from the records, the same helper the reader uses so the
    !> two cannot drift. This used to be a literal naming u..ch4 with only
    !> slot 8 substituted - from SpecCol, a third spelling, uppercase where
    !> the other files are lower - and 18 fixed columns, so a project with
    !> more than four gases wrote none of them here at all. That is what kept
    !> the on-the-fly spectral assessment four-gas: this file is its input.
    call SpectralVarTags(e2sg)
    n_spec = 0
    do j = u, ts
        n_spec = n_spec + 1
        spec_slots(n_spec) = j
    end do
    n_cosp = 0
    do j = w_u, w_ts
        if (j == w_w) cycle
        n_cosp = n_cosp + 1
        cosp_slots(n_cosp) = j
    end do
    do j = firstGas, lastGas
        if (j - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        n_spec = n_spec + 1
        spec_slots(n_spec) = j
        n_cosp = n_cosp + 1
        cosp_slots(n_cosp) = j
    end do

    !> Open output file for binned co-spectra
    BinCospectraPath = BinCospectraDir(1:len_trim(BinCospectraDir)) // String &
                 // BinCospec_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = BinCospectraPath, encoding = 'utf-8')
    write(udf, *)           'normalised_and_exponentially_binned_(co)spectra'
    write(udf, *)           '-----------------------------------------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]_->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_natural_frequency_*_(co)spectrum_/_(co)variance'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', Ambient%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a, i4)')   'number_of_bins_=_', Meth%spec%nbins
    write(udf, '(a, a)')    'tapering_window_=_', RPsetup%tap_win(1:len_trim(RPsetup%tap_win))
    dataline = '#_freq,natural_frequency,normalized_frequency'
    do k = 1, n_spec
        dataline = trim(dataline) // ',f_nat*spec(' &
            // trim(e2sg(spec_slots(k))) // ')/var(' &
            // trim(e2sg(spec_slots(k))) // ')'
    end do
    do k = 1, n_cosp
        dataline = trim(dataline) // ',f_nat*cospec(w_' &
            // trim(e2sg(cosp_slots(k))) // ')/cov(w_' &
            // trim(e2sg(cosp_slots(k))) // ')'
    end do
    write(udf, '(a)') trim(dataline)

    !> Write to output file in csv style
    do i = 1, Meth%spec%nbins
        call clearstr(dataline)
        call WriteDatumInt(bcnt(i), datum, '-9999')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(bnf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (bnf(i) /= error .and. Ambient%WS /= 0d0) then
            call WriteDatumFloat(bnf(i) * (SpecCol(u)%Instr%height - Metadata%d) / Ambient%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if

        !> Spectra
        do k = 1, n_spec
            j = spec_slots(k)
            if (DoSpectrum(j) .and. bnf(i) /= error .and. BinnedSpectrum(i)%of(j) /= error) then
                call WriteDatumFloat(bnf(i) * BinnedSpectrum(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        !> Cospectra
        do k = 1, n_cosp
            j = cosp_slots(k)
            if (DoCospectrum(j) .and. bnf(i) /= error &
                .and. BinnedCospectrum(i)%of(j) /= error) then
                call WriteDatumFloat(bnf(i) * BinnedCospectrum(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
end subroutine WriteOutBinnedCoSpectra

!***************************************************************************
!
! \brief       Write binned ogives on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBinnedOgives(String, bnf, bcnt, BinnedOgive, BinnedCoOgive &
    , DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bcnt(Meth%spec%nbins)
    real(kind = dbl), intent(in) :: bnf(Meth%spec%nbins)
    character(*), intent(in) :: String
    type (SpectralType), intent(in) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType), intent(in) :: BinnedCoOgive(Meth%spec%nbins)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    !> local variables
    integer :: i
    integer :: j
    integer :: k
    integer :: n_spec, n_cosp
    integer :: spec_slots(E2NumVar), cosp_slots(E2NumVar)
    character(64) :: e2sg(E2NumVar)
    character(PathLen) :: BinOgivesPath
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''
    include '../src_common/interfaces.inc'

    !> Same slot list and the same helper as the binned (co)spectra above.
    call SpectralVarTags(e2sg)
    n_spec = 0
    do j = u, ts
        n_spec = n_spec + 1
        spec_slots(n_spec) = j
    end do
    n_cosp = 0
    do j = w_u, w_ts
        if (j == w_w) cycle
        n_cosp = n_cosp + 1
        cosp_slots(n_cosp) = j
    end do
    do j = firstGas, lastGas
        if (j - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        n_spec = n_spec + 1
        spec_slots(n_spec) = j
        n_cosp = n_cosp + 1
        cosp_slots(n_cosp) = j
    end do

    !> Open output file for binned co-spectra
    BinOgivesPath = BinOgivesDir(1:len_trim(BinOgivesDir)) // String &
                 // BinOgives_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = BinOgivesPath, encoding = 'utf-8')
    write(udf, '(a)')       'exponentially_binned_ogives'
    write(udf, '(a)')       '---------------------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]_->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_ogive'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', Ambient%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a, i4)')   'number_of_bins_=_', Meth%spec%nbins
    write(udf, '(a, a)')    'tapering_window_=_', RPsetup%tap_win(1:len_trim(RPsetup%tap_win))
    dataline = '#_freq,natural_frequency,normalized_frequency'
    do k = 1, n_spec
        dataline = trim(dataline) // ',og(' // trim(e2sg(spec_slots(k))) // ')'
    end do
    do k = 1, n_cosp
        dataline = trim(dataline) // ',og(w_' &
            // trim(e2sg(cosp_slots(k))) // ')'
    end do
    write(udf, '(a)') trim(dataline)

    !> Write to output file in csv style
    do i = 1, Meth%spec%nbins
        call clearstr(dataline)

        call WriteDatumInt(bcnt(i), datum, '-9999')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(bnf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (bnf(i) /= error .and. Ambient%WS /= 0d0) then
            call WriteDatumFloat(bnf(i) * (SpecCol(u)%Instr%height - Metadata%d) / Ambient%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if

        !> Ogives
        do k = 1, n_spec
            j = spec_slots(k)
            if (DoSpectrum(j)) then
                call WriteDatumFloat(BinnedOgive(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        !> Co-ogives
        do k = 1, n_cosp
            j = cosp_slots(k)
            if (DoCospectrum(j)) then
                call WriteDatumFloat(BinnedCoOgive(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
end subroutine WriteOutBinnedOgives

!***************************************************************************
!
! \brief       Write full (co)spectra on output file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFullCoSpectra(String, nf, Spectrum, Cospectrum, &
    DoSpectrum, DoCospectrum, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: nf(N/2)
    Type(SpectralType), intent(in) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(in) :: Cospectrum(N/2 + 1)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    character(*), intent(in) :: String
    !> local variables
    integer :: i
    integer :: var
    character(PathLen) :: CospectraPath
    character(LongOutstringLen) :: dataline
    character(LongOutstringLen) :: dataline1
    character(LongOutstringLen) :: dataline2
    character(LongOutstringLen) :: dataline3
    character(DatumLen) :: datum = ''
    !> Wide enough for a record-derived species tag; the historical eight are
    !> three or four characters, but a disambiguated one such as `cos_2` is not.
    character(64) :: e2sg(GHGNumVar)
    include '../src_common/interfaces.inc'

    call LogSayNoAdv('   Writing requested full (co)spectra on output file..')

    !> Column names per slot. Shared with the reader that imports this file
    !> back for the in-situ corrections: a name the two spell differently is
    !> simply not imported, and that gas silently gets no correction factor.
    !> Slots 1-8 keep their literal names, `gas4` included - those strings are
    !> the shipped file format.
    call SpectralVarTags(e2sg)

    !> Open output file for binned co-spectra
    CospectraPath = CospectraDir(1:len_trim(CospectraDir)) // String &
                 // Cospec_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = CospectraPath, encoding = 'utf-8')
    write(udf, '(a)')       'Full_(co)spectra'
    write(udf, '(a)')       '----------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_natural_frequency_*_(co)spectrum_/_(co)variance'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', Ambient%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a)')       'tapering_window_=_SQUARED_(no_tapering_forced_by_EddyFlow.&
                            &_Tapering_is_only_applied_for_binned_(co)spectra)'

    !> First writes (co)variances (header + numbers) and (co)spectra labels
    call clearstr(dataline1)
    call clearstr(dataline2)
    call clearstr(dataline3)
    call AddDatum(dataline1, ',', separator)
    call AddDatum(dataline2, ',', separator)
    call AddDatum(dataline3, 'natural_frequency,normalized_frequency', separator)
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .and. DoSpectrum(var)) then
            call AddDatum(dataline1, 'var(' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
            write(datum, *) Stats%cov(var, var)
            call AddDatum(dataline2, datum, separator)
            call AddDatum(dataline3, 'f_nat*spec(' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
        end if
    end do
    do var = 1, GHGNumVar
        if (RPsetup%out_full_cosp(var) .and. DoCospectrum(var)) then
            call AddDatum(dataline1, 'cov(w_' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
            write(datum, *) Stats%cov(w, var)
            call AddDatum(dataline2, datum, separator)
            call AddDatum(dataline3, 'f_nat*cospec(w_' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
        end if
    end do
    write(udf, '(a)') dataline1(1:len_trim(dataline1) - 1)
    write(udf, '(a)') dataline2(1:len_trim(dataline2) - 1)
    write(udf, '(a)') dataline3(1:len_trim(dataline3) - 1)

    !> Now write spectra and cospectra
    do i = 1, N/2
        call clearstr(dataline)
        !> Frequencies
        call WriteDatumFloat(nf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (nf(i) /= error .and. Ambient%WS /= 0d0) then
            call WriteDatumFloat(nf(i) * (SpecCol(u)%Instr%height - Metadata%d) / Ambient%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if
        !> spectra
        do var = 1, GHGNumVar
            if (RPsetup%out_full_sp(var) .and. DoSpectrum(var)) then
                if (nf(i) /= error .and. Stats%cov(var, var) /= 0d0 .and. Stats%cov(var, var) /= error) then
                    call WriteDatumFloat(nf(i) * Spectrum(i+1)%of(var) / Stats%cov(var, var), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        !> cospectra
        do var = 1, GHGNumVar
            if (RPsetup%out_full_cosp(var) .and. DoCospectrum(var)) then
                if (nf(i) /= error .and. Stats%cov(w, var) /= 0d0 .and. Stats%cov(w, var) /= error) then
                    call WriteDatumFloat(nf(i) * Cospectrum(i+1)%of(var) / Stats%cov(w, var), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
    call LogSay('  Done.')
end subroutine WriteOutFullCoSpectra
