!***************************************************************************
! read_and_apply_cec_ratios.f90
! -----------------------------
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
! \brief   Two-phase CEC helper for the FCC executable.
!
!   ReadNextCecRatio — reads one line from the CEC ratios intermediate file
!                      (written by RP).  Called once per Ex record, regardless
!                      of whether the record will be processed, to keep the
!                      file position in sync with the Ex file.
!
!   ApplyCecRatios   — applies pre-read r_ET / r_Fc to FCC's spectrally-
!                      corrected total fluxes and populates CECFlux.
!                      Called only for records that pass FCC's validity checks.
!
!   Intermediate file format (one line per averaging period):
!     timestamp,do_cec,r_ET_cec,r_Fc_cec
!***************************************************************************

!***************************************************************************
! \brief   Read and parse the next line from the CEC ratios file.
!          Returns error for both ratios on EOF or parse failure.
!***************************************************************************
subroutine ReadNextCecRatio(r_ET, r_Fc)
    use m_fx_global_var
    implicit none
    real(kind = dbl), intent(out) :: r_ET
    real(kind = dbl), intent(out) :: r_Fc

    !> local variables
    character(256) :: line
    character(32)  :: r_ET_field, r_Fc_field
    integer :: io_stat, p1, p2, p3
    logical :: is_open

    r_ET = error
    r_Fc = error

    inquire(unit = ucec, opened = is_open)
    if (.not. is_open) return

    read(ucec, '(a)', iostat = io_stat) line
    if (io_stat /= 0) return   ! EOF

    !> Parse: timestamp,do_cec,r_ET,r_Fc  (four comma-separated fields)
    p1 = index(line, ',')
    if (p1 <= 0) return

    p2 = index(line(p1+1:), ',')
    if (p2 <= 0) return

    p3 = index(line(p1+p2+1:), ',')
    if (p3 <= 0) return

    r_ET_field = adjustl(line(p1+p2+1 : p1+p2+p3-1))
    r_Fc_field = adjustl(line(p1+p2+p3+1:))

    read(r_ET_field, *, iostat = io_stat) r_ET
    if (io_stat /= 0) r_ET = error

    read(r_Fc_field, *, iostat = io_stat) r_Fc
    if (io_stat /= 0) r_Fc = error

end subroutine ReadNextCecRatio

!***************************************************************************
! \brief   Apply pre-computed CEC ratios (r_ET, r_Fc) to FCC's total fluxes.
!          Populates the module-level CECFlux fields.
!
! \param   ET_total  WPL-corrected ET from FCC Flux3 [mmol m-2 s-1]
! \param   Fc_total  WPL-corrected NEE from FCC Flux3 [umol m-2 s-1]
! \param   r_ET      Octant ratio fE/fT from RP (or error)
! \param   r_Fc      Octant ratio fR/fP from RP (or error)
! \param   do_cec    1=H2O+CO2, 2=H2O only, 3=CO2 only
!***************************************************************************
subroutine ApplyCecRatios(ET_total, Fc_total, r_ET, r_Fc, do_cec)
    use m_fx_global_var
    implicit none
    real(kind = dbl), intent(in) :: ET_total
    real(kind = dbl), intent(in) :: Fc_total
    real(kind = dbl), intent(in) :: r_ET
    real(kind = dbl), intent(in) :: r_Fc
    integer, intent(in) :: do_cec

    CECFlux%E_cec    = error
    CECFlux%Tr_cec   = error
    CECFlux%Reco_cec = error
    CECFlux%GPP_cec  = error
    CECFlux%NEE_cec  = error
    CECFlux%r_ET_cec = r_ET
    CECFlux%r_Fc_cec = r_Fc
    CECFlux%ok       = .false.

    if (r_ET == error .and. r_Fc == error) return

    !> H2O partitioning (do_cec = 1 or 2)
    if ((do_cec == 1 .or. do_cec == 2) .and. ET_total /= error) then
        if (r_ET /= error .and. r_ET /= 0d0) then
            CECFlux%E_cec  = ET_total / (1d0 + 1d0 / r_ET)
            CECFlux%Tr_cec = ET_total / (1d0 + r_ET)
        end if
    end if

    !> CO2 partitioning (do_cec = 1 or 3)
    if ((do_cec == 1 .or. do_cec == 3) .and. Fc_total /= error) then
        if (r_Fc /= error .and. abs(r_Fc + 1d0) < 0.05d0) then
            !> Near-singularity: leave Reco/GPP as error
        else if (r_Fc /= error .and. r_Fc /= 0d0) then
            CECFlux%Reco_cec = Fc_total / (1d0 + 1d0 / r_Fc)
            CECFlux%GPP_cec  = Fc_total / (1d0 + r_Fc)
        end if
        if (CECFlux%Reco_cec /= error .and. CECFlux%GPP_cec /= error) then
            CECFlux%NEE_cec = CECFlux%Reco_cec + CECFlux%GPP_cec
        end if
    end if

    CECFlux%ok = .true.

end subroutine ApplyCecRatios
