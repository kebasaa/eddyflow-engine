!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine Storage(PrevStats, prevAmbient)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(StatsType), intent(in) :: PrevStats
    type(AmbientStateType), intent(in) :: prevAmbient
    !> local variables
    integer :: gas
    integer :: wsl
    include '../src_common/interfaces_1.inc'
    real(kind = dbl)  :: seconds
    character(10) tmp_date
    character(5) tmp_time


    !> Latent-heat and ET storage follow the primary water record; a
    !> second hygrometer gets its own storage term but must not
    !> overwrite them.
    wsl = PrimaryWaterSlot()

    write(*, '(a)', advance = 'no') '  Calculating storage terms..'

    !> Check that time periods are consecutive. If not, set storage to error and exit
    call AddDateStep(PrevStats%date, PrevStats%time, tmp_date, tmp_time, DateStep)
    if (tmp_date /= Stats%date .or. tmp_time /= Stats%time) then
        Stor%H  = error
        Stor%LE = error
        Stor%of(firstGas:lastGas)  = error
        return
    end if

    !> Initialization
    Stor%H = 0d0
    Stor%LE = 0d0
    Stor%of(:) = 0d0
    seconds = RPsetup%avrg_len * 6d1

    !> The profile-storage block that used to sit here, commented out, worked
    !> on five named species - stH, stCO2, stH2O, stCH4, stGAS4 - against a
    !> dz(5, MaxProfNodes) read from the project. It was the only reader of
    !> bSetup's profile heights, so those settings were parsed and consumed by
    !> nothing; both are gone. Reviving profile storage wants a height per gas
    !> record, not five fixed rows.

    !> If Stor = 0, it means no profile was available or selected,
    !> so calculate it with 1-point formula
    !> Storage for sensible heat
    if (Stor%H == 0) then
        if(Ambient%RhoCp /= error .and. Ambient%Ta /= error &
            .and. prevAmbient%Ta /= error) then
            Stor%H = Ambient%RhoCp * (Ambient%Ta - prevAmbient%Ta) &
                * E2Col(u)%Instr%height / seconds
        else
            Stor%H = error
        end if
    end if

    !> for all gases
    do gas = firstGas, lastGas
        if (Stor%of(gas) == 0) then
            !> Branched on the species, not the slot: which gas occupies a slot
            !> past the first four is decided by the project at run time, so a
            !> `select case (gas)` over fixed slot numbers would silently skip
            !> every gas beyond them. The primary water slot keeps its branch
            !> even when absent, so that an absent H2O still clears LE and ET
            !> the way it always did.
            if (GasSlotIsWater(gas)) then
                if (Stats%chi(gas) /= error &
                    .and. PrevStats%chi(gas) /= error) then
                    Stor%of(gas) = (Stats%chi(gas) - PrevStats%chi(gas)) &
                        / Ambient%Va * E2Col(gas)%Instr%height / seconds
                    !> Latent heat and evapotranspiration storage follow the
                    !> project's primary H2O only. A second H2O measurement
                    !> gets its own storage term but must not overwrite these.
                    if (gas == wsl) then
                        Stor%LE = Stor%of(wsl) * MW_H2O * Ambient%lambda * 1d-3
                        Stor%ET = Stor%of(wsl) * h2o_to_ET
                    end if
                else
                    Stor%of(gas) = error
                    if (gas == wsl) then
                        Stor%LE = error
                        Stor%ET = error
                    end if
                end if
            else
                if (Stats%chi(gas) /= error &
                    .and. PrevStats%chi(gas) /= error) then
                    Stor%of(gas) = (Stats%chi(gas) - PrevStats%chi(gas)) &
                        / Ambient%Va * 1d-3 &
                        * E2Col(gas)%Instr%height * 1d3 / seconds
                else
                    Stor%of(gas) = error
                end if
            end if
        end if
    end do

    write(*, '(a)') ' Done.'
end subroutine Storage
