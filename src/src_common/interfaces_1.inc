interface
    subroutine fcn(m, npar, x, fvec, fjac, iflag)
        implicit none
        integer, parameter :: dbl   = kind(0.0d0)
        integer, intent(in)            :: m, npar
        real(kind = dbl), intent(in)    :: x(:)
        real(kind = dbl), intent(inout) :: fvec(:)
        real(kind = dbl), intent(out)   :: fjac(:,:)
        integer, intent(inout)         :: iflag
    end subroutine fcn
end interface

interface
    function SwVerFromString(string)
        use m_common_global_var
        implicit none
        !> in/out variables
        character(*), intent(in) :: string
        type(SwVerType) :: SwVerFromString
        !> Local variables
        character(4) :: chunk
    end function
end interface

interface
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
    end subroutine GasFullOutputUnits

    !> Water is carved out of the per-gas unit handling; answered from the
    !> record, so a second water record past the h2o slot is recognised too.
    logical function GasSlotIsWater(gas_slot)
        use m_common_global_var
        implicit none
        integer, intent(in) :: gas_slot
    end function GasSlotIsWater

    !> Full-output column-name stems per gas slot, repeats disambiguated.
    subroutine FullOutputGasTags(tags)
        use m_common_global_var
        implicit none
        character(*), intent(out) :: tags(GHGNumVar)
    end subroutine FullOutputGasTags

    !> (Co)spectra column names per slot. Shared by the writer and by the
    !> reader that imports the file back, which must agree exactly.
    subroutine SpectralVarTags(tags)
        use m_common_global_var
        implicit none
        character(*), intent(out) :: tags(GHGNumVar)
    end subroutine SpectralVarTags

    !> Species names for the spectral reports, one per configured gas slot,
    !> taken from that slot's record. Free text, unlike the labels above.
    subroutine SpectralGasNames(names)
        use m_common_global_var
        implicit none
        character(*), intent(out) :: names(GHGNumVar)
    end subroutine SpectralGasNames

    !> The slot holding the site's primary water, or 0 if there is none.
    !> For one-per-site quantities only; a gas that needs its own water for
    !> a WPL correction reads E2Col(gas)%moist_ref.
    integer function PrimaryWaterSlot()
        use m_common_global_var
        implicit none
    end function PrimaryWaterSlot

    !> Slot the one-per-site water columns are gated on. Header and rows are
    !> written by different routines and must agree, so both ask this.
    integer function PrimaryWaterOutSlot()
        use m_common_global_var
        implicit none
    end function PrimaryWaterOutSlot

    !> The slot holding the site's primary CO2, or 0 if there is none, and
    !> the same falling back to the historical slot. Only conditional eddy
    !> covariance needs these: it is defined on a CO2/water pair and so has
    !> to name two species rather than iterate over all of them.
    integer function PrimaryCarbonSlot()
        use m_common_global_var
        implicit none
    end function PrimaryCarbonSlot

    integer function PrimaryCarbonOutSlot()
        use m_common_global_var
        implicit none
    end function PrimaryCarbonOutSlot

    !> Width and legend of the packed statistical-flag strings. The units row
    !> names the variables and the data row emits their digits, so both sides
    !> ask this or the cell describes a different set than it carries.
    subroutine StatisticalFlagVars(nvars, legend)
        use m_common_global_var
        implicit none
        integer, intent(out) :: nvars
        character(*), intent(out) :: legend
    end subroutine StatisticalFlagVars

    subroutine TimelagFlagLegend(nvars, legend)
        use m_common_global_var
        implicit none
        integer, intent(out) :: nvars
        character(*), intent(out) :: legend
    end subroutine TimelagFlagLegend

    !> Whether a gas has an in-situ spectral correction, or must go through
    !> FCC with BPCF = 1. A species question, not a slot one.
    logical function HasInSituSpectralCorrection(gas_slot)
        use m_common_global_var
        implicit none
        integer, intent(in) :: gas_slot
    end function HasInSituSpectralCorrection

    !> First gas slot the project configures, for the places that need "an"
    !> analyser rather than a particular one.
    integer function FirstConfiguredGasSlot()
        use m_common_global_var
        implicit none
    end function FirstConfiguredGasSlot

    !> Per-slot full-output scales and labels. Header and row writer share
    !> one call so a label can never describe a scale that was not applied.
    subroutine GasFullOutputUnitsAll(flux_sc, dens_sc, &
        flux_label, conc_label, mixr_label, dens_label)
        use m_common_global_var
        implicit none
        real(kind = dbl), intent(out) :: flux_sc(GHGNumVar), dens_sc(GHGNumVar)
        character(*), intent(out) :: flux_label(GHGNumVar), conc_label(GHGNumVar)
        character(*), intent(out) :: mixr_label(GHGNumVar), dens_label(GHGNumVar)
    end subroutine GasFullOutputUnitsAll

    !> FLUXNET column scales, per species rather than per slot.
    real(kind = dbl) function FluxnetGasScale(gas_slot)
        use m_common_global_var
        implicit none
        integer, intent(in) :: gas_slot
    end function FluxnetGasScale

    real(kind = dbl) function FluxnetGasAdvScale(gas_slot)
        use m_common_global_var
        implicit none
        integer, intent(in) :: gas_slot
    end function FluxnetGasAdvScale
end interface

interface
    function CompareSwVer(v1, v2)
        use m_common_global_var
        implicit none
        !> In/out variables
        type(SwVerType), intent(in) :: v1
        type(SwVerType), intent(in) :: v2
        logical :: CompareSwVer
    end function
end interface

interface
    function EqualSwVer(v1, v2)
        use m_common_global_var
        implicit none
        !> In/out variables
        type(SwVerType), intent(in) :: v1
        type(SwVerType), intent(in) :: v2
        logical :: EqualSwVer
    end function
end interface

interface
    function DefaultSonicOutputRate(model)
        use m_common_global_var
        implicit none
        !> In/out variables
        character(*), intent(in) :: model
        integer :: DefaultSonicOutputRate
    end function
end interface

interface
    function sinc(x, N)
        use m_numeric_kinds
        implicit none
        !> In/out variables
        integer, intent(in) :: N
        real(kind = dbl) , intent(in) :: x(N)
        real(kind = dbl)  :: sinc(N)
    end function
end interface

interface
    function asymmetric_linear_transformation(x, N, pgain, poffset, ngain, noffset)
        use m_numeric_kinds
        implicit none
        !> In/out variables
        integer, intent(in) :: N
        real(kind = dbl) , intent(inout) :: x(N)
        real(kind = dbl) , intent(in) :: pgain
        real(kind = dbl) , intent(in) :: poffset
        real(kind = dbl) , intent(in) :: ngain
        real(kind = dbl) , intent(in) :: noffset
        real(kind = dbl)  :: asymmetric_linear_transformation(N)
    end function
end interface

interface
    function CountRecordsAndValues(Set, nrow, ncol, var1, var2)
        use m_common_global_var
        implicit none
        !> In/out variables
        integer, intent(in) :: nrow, ncol
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        integer, optional, intent(in) :: var1, var2
        integer :: CountRecordsAndValues
    end function
end interface

interface
    function PotentialRadiation(latit) result(RP)
        use m_common_global_var
        implicit none
        !> in/out variables
        real(kind = dbl), intent(in) :: latit
        real(kind = dbl) :: RP(17568)
    end function
end interface


interface
    function DateTimeToHalfHourNumber(date, time) result(indx)
        character(*), intent(in) :: date
        character(*), intent(in) :: time
        integer :: indx
    end function
end interface


interface
    function replace(string, what, with, outlen) result(nstring)
        use m_common_global_var
        implicit none
        !> in/out variables
        integer, intent(in) :: outlen
        character(*), intent(in) :: string
        character(*), intent(in) :: what
        character(*), intent(in) :: with
        character(outlen) :: nstring
    end function
end interface

interface
    function replace2(string, what, with) result(nstring)
        use m_common_global_var
        implicit none
        !> in/out variables
        character(*), intent(in) :: string
        character(*), intent(in) :: what
        character(*), intent(in) :: with
        character(:), allocatable :: nstring
    end function
end interface

interface 
    function strCharIndex(s, c, n) result(ix)
        use m_common_global_var
        implicit none
        !> In/out variables
        character(*),intent(in) :: s
        character(*),intent(in) :: c
        integer, intent(in) :: n
        integer :: ix
    end function
end interface

interface
    function strByInt(i) result(res)
        character(:),allocatable :: res
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
    end function
end interface

interface 
    subroutine AddFloatDatumToDataline(float_datum, dataline, err_label, gain, offset)
        use m_common_global_var
        implicit none
        !> in/out variables
        real(kind = dbl), intent(in) :: float_datum
        real(kind = dbl), optional, intent(in) :: gain
        real(kind = dbl), optional, intent(in) :: offset
        character(*), intent(in) :: err_label
        character(*), intent(inout) :: dataline
        character(DatumLen) :: char_datum
    end subroutine
end interface

interface 
    subroutine WriteDatumFloat(float_datum, char_datum, err_label, gain, offset)
        use m_common_global_var
        implicit none
        !> in/out variables
        real(kind = dbl), intent(in) :: float_datum
        real(kind = dbl), optional, intent(in) :: gain
        real(kind = dbl), optional, intent(in) :: offset
        character(*), intent(out) :: char_datum
        character(*), intent(in) :: err_label
        real(kind = dbl) :: val
    end subroutine
end interface
