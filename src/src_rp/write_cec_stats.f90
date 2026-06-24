!***************************************************************************
! write_cec_stats.f90
! -------------------
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
! \brief   Write one CEC-ratios record per averaging period to the
!          intermediate CEC ratios file (_cec_ratios.csv) so that FCC can
!          apply the dimensionless partitioning ratios to its own
!          spectrally-corrected fluxes.
!
!          Called by eddyflow-rp_main.f90 for each period when
!          fcc_follows = .true. and do_cec > 0.
!***************************************************************************
subroutine WriteCecStats()
    use m_rp_global_var
    implicit none
    character(32) :: timestamp
    character(32) :: r_ET_str, r_Fc_str

    timestamp = trim(Stats%date) // 'T' // trim(Stats%time)

    if (CECFlux%r_ET_cec /= error) then
        write(r_ET_str, '(f14.6)') CECFlux%r_ET_cec
    else
        r_ET_str = trim(EddyFlowProj%err_label)
    end if

    if (CECFlux%r_Fc_cec /= error) then
        write(r_Fc_str, '(f14.6)') CECFlux%r_Fc_cec
    else
        r_Fc_str = trim(EddyFlowProj%err_label)
    end if

    write(ucec, '(a,a,i1,a,a,a,a)') &
        trim(timestamp), ',', &
        EddyFlowProj%do_cec, ',', &
        trim(adjustl(r_ET_str)), ',', &
        trim(adjustl(r_Fc_str))

end subroutine WriteCecStats
