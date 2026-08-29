!***************************************************************************
! biomet_units_conversions.f90
! ----------------------------
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
! \brief       Convert input units into standard units
! \author      Gerardo Fratini
! \note
!              Radiations (Rg, Rn, Rd, Rr, LWin, LWout, Ruva, Ruvb) \n
!              are not expected to need unit conversion
!              Photons flux densities (PPFD, PPFDd, PPFDr, PPFDbc, APAR) \n
!              are not expected to need unit conversion
!              Alb, PRI, SWC, SHF are not expected to need unit conversion
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardEddyFlowUnits()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i


    !> Per-channel linear calibration, ahead of the unit conversion below -
    !> a channel is calibrated in its own raw units before those units are
    !> converted, not after. See BiometApplyCalibration's own header.
    call BiometApplyCalibration()

    do i = 1, nbVars
        select case(trim(bVars(i)%nature))
            case('TEMPERATURE')
                select case(bVars(i)%unit_in)
                    case('C','�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) + 273.15d0
                        end where
                    case('F','�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
                        end where
                    case('CK')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('CC','C�C')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2 + 273.15d0
                        end where
                    case('CF','C�F')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = (bSet(:, i) * 1d-2 - 32d0) * 5d0 / 9d0 &
                                + 273.15d0
                        end where
                    case default
                end select

            case('RELATIVE_HUMIDITY')
                select case(bVars(i)%unit_in)
                    case('NUMBER','#','DIMENSIONLESS')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case default
                        continue
                end select

            case('PRESSURE')
                select case(bVars(i)%unit_in)
                    case('HPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d2
                        end where
                    case('KPA')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case('MMHG', 'TORR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 133.322368d0
                        end where
                    case('PSI')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 6894.757d0
                        end where
                    case('BAR')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d5
                        end where
                    case('ATM')
                        !> The standard atmosphere (101325 Pa exactly, by
                        !> definition), not the technical atmosphere
                        !> (98066.5 Pa, 1 kgf/cm^2) the previous constant
                        !> here actually was - a 3.3% error for any channel
                        !> stating ATM.
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 101325d0
                        end where
                    case default
                        continue
                end select

            !> Lengths
            !> converted to [m]
            case('LENGTH', 'PRECIPITATION')
                select case(bVars(i)%unit_in)
                    case('NM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-9
                        end where
                    case('UM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-6
                        end where
                    case('MM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case('CM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('KM')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d3
                        end where
                    case default
                        continue
                end select

            case('CONCENTRATION')
                select case(trim(adjustl(bVars(i)%label)))
                !> CO2 is converted to PPM
                case ('CO2')
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case default
                            continue
                    end select
                !> H2O is converted to PPT
                case ('H2O')
                    select case(bVars(i)%unit_in)
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-3
                            end where
                        case('PPB', 'NMOL/MOL', 'NMOL+1MOL-1', 'NMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d-6
                            end where
                        case default
                            continue
                    end select
                !> Any other gas is converted to PPB
                case default
                    select case(bVars(i)%unit_in)
                        case('PPT', 'MMOL/MOL', 'MMOL+1MOL-1', 'MMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d6
                            end where
                        case('PPM', 'UMOL/MOL', 'UMOL+1MOL-1', 'UMOLMOL-1')
                            where (bSet(:, i) /= error)
                                bSet(:, i) = bSet(:, i) * 1d3
                            end where
                        case default
                            continue
                    end select
                end select

            case('SPEED')
                select case(bVars(i)%unit_in)
                    case('CM+1S-1','CM/S','CMS^-1','CMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-2
                        end where
                    case('MM+1S-1','MM/S','MMS^-1','MMS-1')
                        where (bSet(:, i) /= error)
                            bSet(:, i) = bSet(:, i) * 1d-3
                        end where
                    case default
                        continue
                end select

            case('ANGULAR_DIRECTION')
            case('FLOW')
        end select
    end do
end subroutine BiometStandardEddyFlowUnits

!***************************************************************************
!
! \brief       Per-channel linear calibration on raw biomet values, ahead
!              of unit conversion: bSet(:,i) <- bSet(:,i)*gain + offset.
! \details     Adapted from RFlux's convert_rawdata(), whose per-channel
!              info_* lists carry the same offset+raw*gain term (plus a
!              quadratic one this port does not have - see
!              read_biomet_meta_file.f90's own note on why). Lets a
!              biomet channel's raw counts or voltages become physical
!              readings without a pre-processing step outside EddyFlow,
!              the same way a raw high-frequency channel is calibrated
!              before EddyFlow ever sees the equivalent gain/offset on
!              that side.
!
!              Reaches biomet imported through ReadBiometMetaFile - biomet
!              embedded in a .ghg archive (read_licor_ghg_archive.f90),
!              and, when biom_use_native_header (SCTags(58)) is set to
!              false, an external CSV file too: InitExternalBiomet then
!              reads a sidecar .metadata file next to the biomet file,
!              the same key=value format embedded biomet already uses,
!              instead of deriving bVars from the file's own two-line
!              header (RetrieveExtBiometVars). With biom_use_native_header
!              left at its default (true), external biomet still has no
!              way to state _gain/_offset, so a channel there is assumed
!              already in physical units, same as before this existed.
!
!              A channel with neither stated (bVars(i)%gain still
!              nullbVar's error sentinel) is left untouched - biomet is
!              assumed already in physical units unless the metadata file
!              says otherwise, so a project that never states
!              _gain/_offset sees no change at all.
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          read_biomet_meta_file.f90, RFlux-master/R/convert_rawdata.R
!***************************************************************************
subroutine BiometApplyCalibration()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    real(kind = dbl) :: gain, offset

    do i = 1, nbVars
        if (bVars(i)%gain == error .and. bVars(i)%offset == error) cycle
        gain = bVars(i)%gain
        if (gain == error) gain = 1d0
        offset = bVars(i)%offset
        if (offset == error) offset = 0d0
        where (bSet(:, i) /= error)
            bSet(:, i) = bSet(:, i) * gain + offset
        end where
    end do
end subroutine BiometApplyCalibration

!***************************************************************************
!
! \brief       Start from EddyFlow-standardized units to create
!              dataset with FLUXNET-standardized units
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardFluxnetUnits()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i


    !> Most variables will have same units..
    bAggrFluxnet = bAggrEddyFlow

    !> Change units as needed
    do i = 1, nbVars
        !> All temperatures converted to [degC]
        if (trim(bVars(i)%nature) == 'TEMPERATURE' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) - 273.15d0
        !> Air pressure is converted to [kPa]
        if (bVars(i)%fluxnet_base_name == 'PA' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d-3
        !> VPD is converted to [hPa]
        if (bVars(i)%fluxnet_base_name == 'VPD' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d-2
        !> All precipitations are converted to [mm]
        if (trim(bVars(i)%nature) == 'PRECIPITATION' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d3
        !> Snow depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'SNOW_D' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d2
        !> Water table depth is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'WATER_TABLE_DEPTH' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d2
        !> SWC is converted to [%]
        if (bVars(i)%fluxnet_base_name == 'SWC' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d2
        !> RUNOFF is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'RUNOFF' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d3
        !> THROUGHFALL is converted to [mm]
        if (bVars(i)%fluxnet_base_name == 'THROUGHFALL' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d3
        !> DBH is converted to [cm]
        if (bVars(i)%fluxnet_base_name == 'DBH' .and. bAggrEddyFlow(i) /= error) &
            bAggrFluxnet(i) = bAggrEddyFlow(i) * 1d2
    end do
end subroutine BiometStandardFluxnetUnits
