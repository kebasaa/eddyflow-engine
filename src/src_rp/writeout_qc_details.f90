!***************************************************************************
! writeout_qc_details.f90
! -----------------------
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
! \brief       Write results of stationarity and integral turbulence test on
!              file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutQCDetails(string, StDiff, DtDiff, STFlg, DTFlg)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    integer, intent(in) :: STFlg(GHGNumVar)
    integer, intent(in) :: DTFlg(GHGNumVar)
    !> local variables
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum

    call clearstr(dataline)

    !> File info
    call AddDatum(dataline, string(1:len_trim(string)), separator)
    !> Differences (%) of stationarity test
    call WriteDatumInt(STDiff%u, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%w, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%ts, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (OutVarPresent(co2)) then
        call WriteDatumInt(STDiff%co2, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STDiff%h2o, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STDiff%ch4, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STDiff%gas4, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    call WriteDatumInt(STDiff%w_u, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%w_ts, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)

    if (OutVarPresent(co2)) then
        call WriteDatumInt(STDiff%w_co2, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STDiff%w_h2o, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STDiff%w_ch4, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STDiff%w_gas4, datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    !> Flags for the stationarity test
    call WriteDatumInt(STFlg(w_u), datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_ts), datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (OutVarPresent(co2)) then
        call WriteDatumInt(STFlg(w_co2), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STFlg(w_h2o), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STFlg(w_ch4), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STFlg(w_gas4), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    !> Differences (%) of the well developed turbulence test
    call WriteDatumInt(DtDiff%u, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%w, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%ts, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> Flags for the well developed turbulence test
    call WriteDatumInt(DTFlg(u), datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(w), datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(ts), datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)

    write(uqc, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutQCDetails
