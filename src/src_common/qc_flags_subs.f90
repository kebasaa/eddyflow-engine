!***************************************************************************
! qc_flags_subs.f90
! -----------------
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
! \brief       Quality flagging according to the chosen method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine QualityFlags(lFlux2, StDiff, DtDiff, STFlg, DTFlg, lQCFlag, printout, raw_ok, rf_ok)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(FluxType), intent(in) :: lFlux2
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    logical, intent(in) :: printout
    !> Whether Essentials%KID/mahrt98_NR reflect the period being flagged
    !> right now, and whether Essentials%AL1/DDI/HF5/HF10/HD5/HD10/DIP do
    !> too - both feed 'vitale_20' only. True from RP's own call, which
    !> just computed them; false from FCC's call, which never re-derives
    !> per-sample raw-signal diagnostics from the ex record it reads back
    !> - it re-emits test_rf's own columns verbatim rather than parsing
    !> them (see pfd_handle.f90's sibling note on the same record).
    !> rf_ok additionally needs RP's own Test%rf, a src_rp-only type this
    !> src_common file cannot reference directly, so the caller resolves
    !> it instead.
    !>
    !> This matters more than "RP's own call vs FCC's own call" suggests:
    !> whenever EddyFlowProj%fcc_follows is set - which it is for any
    !> project using a spectral-correction method beyond the two built-in
    !> analytic ones (moncrieff_97, massman_00), i.e. most real projects -
    !> RP's own QualityFlags call is the one that gets skipped (the
    !> `else` branch below sets every flag to `error` and RP does not
    !> even write full_output in that mode), and FCC's call - raw_ok and
    !> rf_ok both false, always - is the one whose output the user
    !> actually sees. So for most projects 'vitale_20' runs as ITC-only
    !> in practice, not as the wider aggregate the tooltip's first line
    !> promises; see VitaleFlag's own header.
    logical, intent(in) :: raw_ok
    logical, intent(in) :: rf_ok
    type(QCType), intent(out)   :: lQCFlag
    integer, intent(out) :: STFlg(GHGNumVar)
    integer, intent(out) :: DTFlg(GHGNumVar)
    !> local variables
    integer :: gas


    if (printout) write(*,'(a)', advance = 'no') '  Calculating quality flags..'
    if (printout) write(ulog,'(a)', advance = 'no') '  Calculating quality flags..'

    !> Stationarity flags. Every gas slot: the w_* covariance enumeration is
    !> numerically the variable enumeration, so a slot indexes both directly.
    do gas = firstGas, lastGas
        call PartialFlagLF(StDiff%w_gas(gas), STFlg(gas))
    end do
    call PartialFlagLF(StDiff%w_ts,  STFlg(w_ts))
    call PartialFlagLF(StDiff%w_u,   STFlg(w_u))
    !> Developed turbulence flags
    call PartialFlagLF(DtDiff%u, DTFlg(u))
    call PartialFlagLF(DtDiff%w, DTFlg(w))
    call PartialFlagLF(DtDiff%ts, DTFlg(ts))
    DTFlg(u)  = max(DTFlg(u),  DTFlg(w))

    if (.not. EddyFlowProj%fcc_follows) then
        select case(Meth%qcflag(1:len_trim(Meth%qcflag)))
            case ('none')
                lQCFlag%tau = nint(error)
                lQCFlag%H = nint(error)
                lQCFlag%gas = nint(error)
            case ('mauder_foken_04')
                !> Combined flags according to Mauder and Foken (2004)
                call GTK2Flag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call GTK2Flag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                do gas = firstGas, lastGas
                    call GTK2Flag(STFlg(gas), DTFlg(w), lQCFlag%gas(gas))
                end do
            case ('foken_03')
                !> Combined flags according to Foken (2003), retrieved from Foken et al. (2004, HoM)
                call FokenFlag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call FokenFlag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                do gas = firstGas, lastGas
                    call FokenFlag(STFlg(gas), DTFlg(w), lQCFlag%gas(gas))
                end do
            case ('goeckede_06')
                !> Combined flags according to Goeckede et al. (2006)
                call GoeckedeFlag(STFlg(w_u),   DTFlg(u), lQCFlag%tau)
                call GoeckedeFlag(STFlg(w_ts),  DTFlg(w), lQCFlag%H)
                do gas = firstGas, lastGas
                    call GoeckedeFlag(STFlg(gas), DTFlg(w), lQCFlag%gas(gas))
                end do
            case ('vitale_20')
                !> Two-tier severity scheme (Vitale et al. 2020, Biogeosciences),
                !> from a wider net of tests than the three above: see
                !> VitaleFlag's own header for what feeds it.
                !> RFlux's own ITC test (cleanFlux.R) is defined once, on
                !> itc_w alone, and that same SevEr/ModEr set is reused for
                !> all four fluxes including TAU - so this passes DTFlg(w),
                !> not the max(u*,w) DTFlg(u) the other three methods above
                !> use for their own, unrelated TAU convention.
                call VitaleFlag(DTFlg(w), (/u, v, w/), raw_ok, rf_ok, lQCFlag%tau)
                call VitaleFlag(DTFlg(w), (/ts, w, 0/), raw_ok, rf_ok, lQCFlag%H)
                do gas = firstGas, lastGas
                    call VitaleFlag(DTFlg(w), (/gas, w, 0/), raw_ok, rf_ok, lQCFlag%gas(gas))
                end do
        end select
        !> If fluxes are set to error, set to error also the quality flags
        if (lFlux2%H    == error) lQCFlag%H    = nint(error)
        do gas = firstGas, lastGas
            if (lFlux2%gas(gas) == error) lQCFlag%gas(gas) = nint(error)
        end do
    else
        lQCFlag%tau = nint(error)
        lQCFlag%H = nint(error)
        lQCFlag%gas = nint(error)
    end if

    if (printout) write(*, '(a)') ' Done.'
    if (printout) write(ulog, '(a)') ' Done.'
end subroutine QualityFlags

!***************************************************************************
!
! \brief       Partial flags, after Foken et al. (2004, Handbook of Microm.)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PartialFlagLF(val, flag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: val
    integer, intent(out) :: flag

    select case (val)
        case (0:15)
            flag = 1
        case (16:30)
            flag = 2
        case (31:50)
            flag = 3
        case (51:75)
            flag = 4
        case (76:100)
            flag = 5
        case (101:250)
            flag = 6
        case (251:500)
            flag = 7
        case (501:1000)
            flag = 8
        case (1001:)
            flag = 9
        case default
            flag = ierror
    end select
end subroutine PartialFlagLF

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
!subroutine GhgEuropeFlagLF(STFlg, DTFlg, OAFlag)
!    use m_common_global_var
!    implicit none
!    !> in/out variables
!    integer, intent(in) :: STFlg
!    integer, intent(in) :: DTFlg
!    integer, intent(out) :: OAFlag
!
!
!    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
!        OAFlag = idint(error)
!        return
!    end if
!    OAFlag = nint(error)
!
!    !> Flags combination: it might need a feedback from Foken's group.
!    if( STFlg >= 6 .and. DTFlg >= 6)                                        OAFlag = 2
!    if((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 1
!    if((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 0
!    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 0
!end subroutine GhgEuropeFlagLF

!***************************************************************************
!
! \brief       Final flags, after Mauder and Foken (2004), TK2 documentation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GTK2Flag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag


    !> stat test < 30  ==> stat flag <= 2
    !> itc test  < 30  ==> itc flag  <= 2
    !> stat test < 100 ==> stat flag <= 5
    !> itc test  < 100 ==> itc  flag <= 5
    if (STFlg == ierror .or. DTFlg == ierror) then
        OAFlag = 2
        return
    end if

    OAFlag = 2
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 2)) then
        OAFlag = 0
    elseif ((STFlg >= 1 .and. STFlg <= 5) .and. (DTFlg >= 1 .and. DTFlg <= 5)) then
        OAFlag = 1
    else
        OAFlag = 2
    end if
end subroutine GTK2Flag

!***************************************************************************
!
! \brief       Partial flags, after Gockede et al. (2004, AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GoeckedeFlag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag

    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
        OAFlag = 5
        return
    end if

    !> From Table 2 Goecked et al. (2004, AFM)
    !> Note that the range (STFlg >= 5 .and. STFlg <= 9) .and. (DTFlg >= 5 .and. DTFlg <= 6) is not
    !> provided in the paper, here it is included in the flag 5 (last).
    !> Other ranges not supported are left with error codes
    OAFlag = 5
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 1
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 2
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 3
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 5 .and. DTFlg <= 6))   OAFlag = 4
    if((STFlg >= 5 .and. STFlg <= 9) .and. (DTFlg >= 5 .and. DTFlg <= 9))   OAFlag = 5
end subroutine GoeckedeFlag

!***************************************************************************
!
! \brief       Partial flags, after Foken et al. (2003) as retrieved from
!              Foken et al. 2004, Handbook of Micrometeorology, Table 9.4
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FokenFlag(STFlg, DTFlg, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: STFlg
    integer, intent(inout) :: DTFlg
    integer, intent(out) :: OAFlag

    if (STFlg == idint(error) .or. DTFlg == idint(error)) then
        OAFlag = 9
        return
    end if

    !> From Table 9.4 Foken et al. 2004, HoM
    !> Note that values 7 and 8 of the overall flag was somewhat interpreted from the text
    !> of the Table, which is ambiguous (ranges of flags 7 and 8 overlap with previous ones)
    OAFlag = 9
    if(STFlg == 1                    .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 1
    if(STFlg == 2                    .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 2
    if((STFlg >= 1 .and. STFlg <= 2) .and. (DTFlg >= 3 .and. DTFlg <= 4))   OAFlag = 3
    if((STFlg >= 3 .and. STFlg <= 4) .and. (DTFlg >= 1 .and. DTFlg <= 2))   OAFlag = 4
    if((STFlg >= 1 .and. STFlg <= 4) .and. (DTFlg >= 3 .and. DTFlg <= 5))   OAFlag = 5
    if(STFlg == 5                    .and. (DTFlg >= 1 .and. DTFlg <= 5))   OAFlag = 6

    if(STFlg == 6                    .and. (DTFlg >= 1 .and. DTFlg <= 6))   OAFlag = 7
    if((STFlg >= 1 .and. STFlg <= 6) .and. DTFlg == 6)                      OAFlag = 7

    if((STFlg >= 7 .and. STFlg <= 8) .and. (DTFlg >= 1 .and. DTFlg <= 8))   OAFlag = 8
    if((STFlg >= 1 .and. STFlg <= 8) .and. (DTFlg >= 7 .and. DTFlg <= 8))   OAFlag = 8

    if(STFlg == 9                    .or.  DTFlg == 9)                      OAFlag = 9
end subroutine FokenFlag

!***************************************************************************
!
! \brief       Two-tier severity flag (Vitale et al. 2020, Biogeosciences):
!              0 (ok), 1 (moderate, RFlux's "ModEr") or 2 (severe, "SevEr").
! \details     Built entirely from diagnostics EddyFlow already computes,
!              regardless of which quality-flagging method is selected:
!              - itc_flg is DTFlg(u) for TAU, DTFlg(w) for H and every gas -
!                the same grade GTK2Flag/FokenFlag/GoeckedeFlag above use.
!                Vitale et al. (2020)'s own 30%/100% cutoffs happen to land
!                exactly on two of PartialFlagLF's own grade boundaries
!                (16-30 and 31-100 split there), so the grade is reused
!                rather than re-deriving the percentage from scratch.
!              - Essentials%mahrt98_NR(vars(1)), Mahrt's (1998)
!                nonstationarity ratio, computed unconditionally by
!                RU_Mahrt_98 alongside random uncertainty. vars(1) is the
!                flux's own slot in that array: u for TAU, ts for H, the
!                gas slot for a gas flux - the same indexing AL1/DDI/KID
!                use below, since RU_Mahrt_98 fills it the same way.
!              - Essentials%KID(:), computed unconditionally.
!              - Essentials%AL1/DDI/HF5/HF10/HD5/HD10/DIP(:), only when
!                RP's own Test%rf is on (otherwise `error` - RFlux's own
!                extra raw-signal diagnostic suite this ports).
!              vars(2:3) are the other raw variables RFlux's own SevEr/
!              ModEr unions also check for this flux (w always; v as well
!              for TAU); 0 marks an unused slot.
!
!              raw_ok and rf_ok say whether Essentials reflects the period
!              being flagged right now at all: true from RP, which just
!              computed it; false from FCC, which never re-derives
!              per-sample raw-signal diagnostics from the ex record it
!              reads back (see this argument's own note on QualityFlags).
!              Without this gate a period run through FCC would read
!              Essentials's default-initialised (not `error`) zeros as
!              genuine measurements - AL1 <= 0.5 would misfire as severe
!              on literally every FCC-computed period.
!
!              Deliberately narrower than RFlux's cleanFlux(): its
!              low-signal-resolution test (LSR) and its physical-range
!              filter (hard-coded CO2/H2O/sensible-heat/momentum limits)
!              have no EddyFlow port, and its wind-sector exclusion is left
!              out until EddyFlow supports a time-varying sector - both
!              genuine gaps against RFlux, not oversights here.
!
!              Essentials%CCF(:) (the CCF-degradation R^2 test,
!              cross_corr_test.f90) is deliberately NOT read here either,
!              even though it is computed under the same Test%rf gate as
!              AL1/DDI/etc.: RFlux's own cleanFlux.R never consumes its
!              lrt_h/lrt_fc/lrt_le/lrt_tau values in SevEr/ModEr either -
!              qcStat.R exports them purely as diagnostic columns. Feeding
!              it into this flag would be inventing a threshold RFlux
!              itself never defines, not porting one.
!
!              Under fcc_follows this grade degrades to the ITC deviation
!              alone, for the reason raw_ok/rf_ok's own comment on
!              QualityFlags gives - and fcc_follows is set for any project
!              using a spectral-correction method beyond the two built-in
!              analytic ones, which is most real projects, not a rare
!              corner case. Properly closing this needs FCC's ex record to
!              carry KID/mahrt98_NR/AL1/DDI/HF5/HF10/HD5/HD10/DIP the way
!              it already carries the ITC deviation (U_ITC/W_ITC/TS_ITC)
!              and the stationarity percentages (TAU_SS/H_SS/F_SS) -
!              a real follow-up, not attempted here.
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          RFlux-master/R/cleanFlux.R (the SevEr_ind/ModEr_ind unions)
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine VitaleFlag(itc_flg, vars, raw_ok, rf_ok, OAFlag)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: itc_flg
    integer, intent(in) :: vars(3)
    logical, intent(in) :: raw_ok
    logical, intent(in) :: rf_ok
    integer, intent(out) :: OAFlag
    !> local variables
    logical :: sev_flag, mod_flag
    real(kind = dbl) :: hz, scth1, scth2, m98
    integer :: islot, vslot

    sev_flag = .false.
    mod_flag = .false.

    if (itc_flg > 0) then
        sev_flag = sev_flag .or. (itc_flg >= 6)
        mod_flag = mod_flag .or. (itc_flg >= 3 .and. itc_flg <= 5)
    end if

    if (raw_ok) then
        m98 = Essentials%mahrt98_NR(vars(1))
        if (m98 /= error) then
            sev_flag = sev_flag .or. (m98 > 3d0)
            mod_flag = mod_flag .or. (m98 > 2d0 .and. m98 <= 3d0)
        end if
    end if

    hz = Metadata%ac_freq
    scth1 = hz * 60d0 * 30d0 * 0.04d0
    scth2 = hz * 60d0 * 30d0 * 0.01d0

    do islot = 1, 3
        vslot = vars(islot)
        if (vslot <= 0) cycle
        if (.not. raw_ok) cycle

        if (Essentials%KID(vslot) /= error) then
            sev_flag = sev_flag .or. (Essentials%KID(vslot) > 50d0)
            mod_flag = mod_flag .or. (Essentials%KID(vslot) > 30d0 .and. Essentials%KID(vslot) <= 50d0)
        end if

        if (.not. rf_ok) cycle

        if (Essentials%AL1(vslot) /= error) then
            sev_flag = sev_flag .or. (Essentials%AL1(vslot) <= 0.5d0)
            mod_flag = mod_flag .or. (Essentials%AL1(vslot) > 0.5d0 .and. Essentials%AL1(vslot) <= 0.75d0)
        end if
        if (Essentials%DDI(vslot) /= error) then
            sev_flag = sev_flag .or. (Essentials%DDI(vslot) >= hz * 300d0)
            mod_flag = mod_flag .or. (Essentials%DDI(vslot) >= hz * 150d0 .and. Essentials%DDI(vslot) < hz * 300d0)
        end if
        if (Essentials%HF5(vslot) /= error .and. Essentials%HF10(vslot) /= error) then
            sev_flag = sev_flag .or. (Essentials%HF5(vslot) > scth1 .or. Essentials%HF10(vslot) > scth2)
            mod_flag = mod_flag .or. (Essentials%HF5(vslot) > scth1 / 2d0 .or. Essentials%HF10(vslot) > scth2 / 2d0)
        end if
        if (Essentials%HD5(vslot) /= error .and. Essentials%HD10(vslot) /= error) then
            sev_flag = sev_flag .or. (Essentials%HD5(vslot) > scth1 .or. Essentials%HD10(vslot) > scth2)
            mod_flag = mod_flag .or. (Essentials%HD5(vslot) > scth1 / 2d0 .or. Essentials%HD10(vslot) > scth2 / 2d0)
        end if
        if (Essentials%DIP(vslot) /= error) then
            !> DIP is Hartigan's dip-test p-value, in [0,1] - RFlux's own
            !> thresholds (inst_prob_test.R: sev_thr[11]=0.01,
            !> mod_thr[11]=0.05) were 0.0 and 0.1 here, so severe could
            !> never fire (a p-value is never < 0) and moderate fired 5x
            !> too wide.
            sev_flag = sev_flag .or. (Essentials%DIP(vslot) < 0.01d0)
            mod_flag = mod_flag .or. (Essentials%DIP(vslot) >= 0.01d0 .and. Essentials%DIP(vslot) <= 0.05d0)
        end if
    end do

    if (sev_flag) then
        OAFlag = 2
    else if (mod_flag) then
        OAFlag = 1
    else
        OAFlag = 0
    end if
end subroutine VitaleFlag
