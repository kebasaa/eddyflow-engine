!***************************************************************************
! m_cec.f90
! ---------
! Shared Conditional Eddy Covariance implementation following:
! Zahn et al. (2022), Agricultural and Forest Meteorology 315, 108790.
! https://doi.org/10.1016/j.agrformet.2021.108790
!***************************************************************************
module m_cec
    use ieee_arithmetic
    use m_common_global_var
    implicit none
    private

    public :: ResetCecDescriptor
    public :: ResetCecFlux
    public :: ExtractCecDescriptor
    public :: ApplyCecDescriptor

contains

subroutine ResetCecDescriptor(descriptor)
    type(CECDescriptorType), intent(out) :: descriptor

    descriptor%r_ET = error
    descriptor%r_Fc = error
    descriptor%frac_O1 = error
    descriptor%frac_O2 = error
    descriptor%n_valid = 0
    descriptor%n_O1 = 0
    descriptor%n_O2 = 0
    descriptor%h2o_status = cec_rejected
    descriptor%co2_status = cec_rejected
    descriptor%h2o_valid = .false.
    descriptor%co2_valid = .false.
end subroutine ResetCecDescriptor

subroutine ResetCecFlux(flux)
    type(CECFluxType), intent(out) :: flux

    flux%E_cec = error
    flux%Tr_cec = error
    flux%E_cec_ET = error
    flux%Tr_cec_ET = error
    flux%Reco_cec = error
    flux%P_cec = error
    flux%GPP_cec = error
    flux%NEE_cec = error
    flux%r_ET_cec = error
    flux%r_Fc_cec = error
    flux%ok = .false.
end subroutine ResetCecFlux

subroutine ExtractCecDescriptor(primes, stationarity_co2, stationarity_h2o, descriptor)
    real(kind = dbl), intent(in) :: primes(:, :)
    integer, intent(in) :: stationarity_co2
    integer, intent(in) :: stationarity_h2o
    type(CECDescriptorType), intent(out) :: descriptor

    integer :: i
    integer :: nrow
    real(kind = dbl), allocatable :: w_prime(:)
    real(kind = dbl), allocatable :: c_prime(:)
    real(kind = dbl), allocatable :: q_prime(:)
    real(kind = dbl) :: sum_fE
    real(kind = dbl) :: sum_fT
    real(kind = dbl) :: sum_fR
    real(kind = dbl) :: sum_fP
    real(kind = dbl) :: f_E
    real(kind = dbl) :: f_T
    real(kind = dbl) :: f_R
    real(kind = dbl) :: f_P

    call ResetCecDescriptor(descriptor)

    if (size(primes, 2) < h2o) return
    nrow = size(primes, 1)
    if (nrow < 2) return

    allocate(w_prime(nrow), c_prime(nrow), q_prime(nrow))
    w_prime = primes(:, w)
    c_prime = primes(:, co2)
    q_prime = primes(:, h2o)

    call InterpolateShortCecGaps(w_prime, 4)
    call InterpolateShortCecGaps(c_prime, 4)
    call InterpolateShortCecGaps(q_prime, 4)

    sum_fE = 0d0
    sum_fT = 0d0
    sum_fR = 0d0
    sum_fP = 0d0

    do i = 1, nrow
        if (.not. CecValueIsValid(w_prime(i))) cycle
        if (.not. CecValueIsValid(c_prime(i))) cycle
        if (.not. CecValueIsValid(q_prime(i))) cycle

        descriptor%n_valid = descriptor%n_valid + 1
        if (w_prime(i) > 0d0 .and. q_prime(i) > 0d0 .and. c_prime(i) > 0d0) then
            descriptor%n_O1 = descriptor%n_O1 + 1
            sum_fE = sum_fE + w_prime(i) * q_prime(i)
            sum_fR = sum_fR + w_prime(i) * c_prime(i)
        else if (w_prime(i) > 0d0 .and. q_prime(i) > 0d0 .and. c_prime(i) < 0d0) then
            descriptor%n_O2 = descriptor%n_O2 + 1
            sum_fT = sum_fT + w_prime(i) * q_prime(i)
            sum_fP = sum_fP + w_prime(i) * c_prime(i)
        end if
    end do

    !> Zahn et al. retained periods with at least 90% instantaneous data.
    if (10 * descriptor%n_valid < 9 * nrow) return
    if (stationarity_co2 == ierror .or. stationarity_h2o == ierror) return
    if (stationarity_co2 > 25 .or. stationarity_h2o > 25) return

    descriptor%frac_O1 = dble(descriptor%n_O1) / dble(descriptor%n_valid)
    descriptor%frac_O2 = dble(descriptor%n_O2) / dble(descriptor%n_valid)
    if (descriptor%frac_O1 + descriptor%frac_O2 < 0.20d0) return

    f_E = sum_fE / dble(descriptor%n_valid)
    f_T = sum_fT / dble(descriptor%n_valid)
    f_R = sum_fR / dble(descriptor%n_valid)
    f_P = sum_fP / dble(descriptor%n_valid)

    if (f_T /= 0d0) descriptor%r_ET = f_E / f_T
    if (f_P /= 0d0) descriptor%r_Fc = f_R / f_P

    if (descriptor%frac_O1 < 0.05d0 .or. descriptor%r_ET == 0d0) then
        descriptor%h2o_status = cec_all_stomatal
        descriptor%h2o_valid = .true.
    else if (descriptor%frac_O2 < 0.05d0) then
        descriptor%h2o_status = cec_all_nonstomatal
        descriptor%h2o_valid = .true.
    else if (descriptor%r_ET /= error) then
        descriptor%h2o_status = cec_normal
        descriptor%h2o_valid = .true.
    end if

    if (descriptor%r_Fc /= error .and. abs(descriptor%r_Fc + 1d0) < 0.05d0) then
        descriptor%co2_status = cec_singular
    else if (descriptor%frac_O1 < 0.05d0 .or. descriptor%r_Fc == 0d0) then
        descriptor%co2_status = cec_all_stomatal
        descriptor%co2_valid = .true.
    else if (descriptor%frac_O2 < 0.05d0) then
        descriptor%co2_status = cec_all_nonstomatal
        descriptor%co2_valid = .true.
    else if (descriptor%r_Fc /= error) then
        descriptor%co2_status = cec_normal
        descriptor%co2_valid = .true.
    end if
end subroutine ExtractCecDescriptor

subroutine ApplyCecDescriptor(descriptor, H2O_total, Fc_total, do_cec, flux)
    type(CECDescriptorType), intent(in) :: descriptor
    real(kind = dbl), intent(in) :: H2O_total
    real(kind = dbl), intent(in) :: Fc_total
    integer, intent(in) :: do_cec
    type(CECFluxType), intent(out) :: flux

    call ResetCecFlux(flux)
    flux%r_ET_cec = descriptor%r_ET
    flux%r_Fc_cec = descriptor%r_Fc

    if ((do_cec == 1 .or. do_cec == 2) .and. descriptor%h2o_valid &
        .and. H2O_total /= error) then
        select case (descriptor%h2o_status)
            case (cec_normal)
                flux%E_cec = H2O_total / (1d0 + 1d0 / descriptor%r_ET)
                flux%Tr_cec = H2O_total / (1d0 + descriptor%r_ET)
            case (cec_all_stomatal)
                flux%E_cec = 0d0
                flux%Tr_cec = H2O_total
            case (cec_all_nonstomatal)
                flux%E_cec = H2O_total
                flux%Tr_cec = 0d0
        end select
        if (flux%E_cec /= error .and. flux%Tr_cec /= error) then
            flux%E_cec_ET = flux%E_cec * h2o_to_ET
            flux%Tr_cec_ET = flux%Tr_cec * h2o_to_ET
        end if
    end if

    if ((do_cec == 1 .or. do_cec == 3) .and. Fc_total /= error) &
        flux%NEE_cec = Fc_total

    if ((do_cec == 1 .or. do_cec == 3) .and. descriptor%co2_valid &
        .and. Fc_total /= error) then
        select case (descriptor%co2_status)
            case (cec_normal)
                flux%Reco_cec = Fc_total / (1d0 + 1d0 / descriptor%r_Fc)
                flux%P_cec = Fc_total / (1d0 + descriptor%r_Fc)
            case (cec_all_stomatal)
                flux%Reco_cec = 0d0
                flux%P_cec = Fc_total
            case (cec_all_nonstomatal)
                flux%Reco_cec = Fc_total
                flux%P_cec = 0d0
        end select
        flux%GPP_cec = flux%P_cec
    end if

    if (do_cec == 1) then
        flux%ok = descriptor%h2o_valid .and. descriptor%co2_valid
    else if (do_cec == 2) then
        flux%ok = descriptor%h2o_valid
    else if (do_cec == 3) then
        flux%ok = descriptor%co2_valid
    end if
end subroutine ApplyCecDescriptor

subroutine InterpolateShortCecGaps(values, max_gap)
    real(kind = dbl), intent(inout) :: values(:)
    integer, intent(in) :: max_gap

    integer :: first_gap
    integer :: gap_length
    integer :: i
    integer :: k
    real(kind = dbl) :: increment

    i = 2
    do while (i < size(values))
        if (CecValueIsValid(values(i))) then
            i = i + 1
            cycle
        end if

        first_gap = i
        do while (i <= size(values))
            if (CecValueIsValid(values(i))) exit
            i = i + 1
        end do
        gap_length = i - first_gap

        if (i <= size(values)) then
            if (gap_length <= max_gap .and. CecValueIsValid(values(first_gap - 1)) &
                .and. CecValueIsValid(values(i))) then
                increment = (values(i) - values(first_gap - 1)) / dble(gap_length + 1)
                do k = 1, gap_length
                    values(first_gap + k - 1) = values(first_gap - 1) + increment * dble(k)
                end do
            end if
        end if
    end do
end subroutine InterpolateShortCecGaps

logical function CecValueIsValid(value)
    real(kind = dbl), intent(in) :: value

    CecValueIsValid = ieee_is_finite(value) .and. value /= error
end function CecValueIsValid

end module m_cec
