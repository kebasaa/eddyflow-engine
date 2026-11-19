!***************************************************************************
! m_fx_global_var.f90
! -------------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Module for global variables in EddyFlow_fcc
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_fx_global_var
    use m_common_global_var
    implicit none
    save

    !> global variables
    integer :: nstat
    integer :: nfull
    integer, parameter :: ndkf = 60     !< TO BE INCREASED????!!!!!!!
    integer, parameter :: unstable = 1
    integer, parameter :: stable   = 2

    real(kind = dbl) :: float_doy
    real(kind = dbl) :: dkf(ndkf + 1)
    !> Full-output unit basis per gas slot. These were single scalars for the
    !> fourth slot, which is why the FCC full output carried no columns at all
    !> for gases past it: there was nothing to label them with.
    real(kind = dbl) :: gas_full_flux_sc(GHGNumVar)
    real(kind = dbl) :: gas_full_dens_sc(GHGNumVar)
    real(kind = dbl), allocatable :: custVars(:)

    character(12), parameter :: fcc_app = 'EddyFlow-FCC'
    character(32) :: gas_full_flux_label(GHGNumVar)
    character(32) :: gas_full_conc_label(GHGNumVar)
    character(32) :: gas_full_mixr_label(GHGNumVar)
    character(32) :: gas_full_dens_label(GHGNumVar)
    character(64) :: UserVarHeader(MaxUserVar)
    character(26), parameter :: SubDirSpecAn = 'eddyflow_spectral_analysis'
    character(LongOutstringLen) :: fluxnet_header

    logical :: MeanBinSpecAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanBinCospAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanStabSpecAvailable(MaxGasClasses, GHGNumVar)
    logical :: MeanStabCospAvailable(MaxGasClasses, GHGNumVar)
    logical :: fcc_var_present(GHGNumVar)

    !> Counters used to explain spectral-assessment eligibility and failures.
    integer :: SADiagSelectedFiles
    integer :: SADiagReadableFiles
    integer :: SADiagMatchedRecords
    integer :: SADiagUsableWT
    integer :: SADiagRejectedUstar
    integer :: SADiagRejectedVM(GHGNumVar)
    integer :: SADiagRejectedFoken(GHGNumVar)
    integer :: SADiagRejectedFlux(GHGNumVar)
    integer :: SADiagAccepted(GHGNumVar)
    integer :: SADiagDegradedUnstable
    integer :: SADiagDegradedStable
    integer, parameter :: SADiagUnstable = 1
    integer, parameter :: SADiagStable = 2
    integer :: SADiagFluxCandidateCount(2, GHGNumVar)
    integer :: SADiagFluxCandidateCapacity
    real(kind = dbl), allocatable :: SADiagFluxCandidate(:, :, :)
    integer, allocatable :: SADiagFluxCandidateClass(:, :, :)
    real(kind = dbl) :: SAAutoMin(2, GHGNumVar)
    real(kind = dbl) :: SAAutoMax(GHGNumVar)
    logical :: SAAutoApplyMin(2, GHGNumVar)
    logical :: SAAutoApplyMax(GHGNumVar)
    character(PathLen) :: SADiagFilePath
    !> The tube-mate ReportSpectralSources suggests each gas take its transfer
    !> function from, 0 for none. ApplyAutomaticSpectralConfiguration writes
    !> it into the output project as gas_<i>_sa_source.
    integer :: SAAutoSource(GHGNumVar) = 0
    !> The tube-mate a gas left on automatic was corrected with in every month
    !> of this run, because sa_auto_source found its own fit implausible; 0
    !> for none. SaveAppliedSpectralSources writes it into the output project,
    !> so the project kept with the results says what produced them.
    integer :: SAAppliedSource(GHGNumVar) = 0

    !> Every suggestion the run makes, gathered as it goes and printed together
    !> as its very last message (ReportRunSuggestions). The readiness report
    !> that makes most of them comes before the period loop, so they used to
    !> sit under thousands of period lines.
    integer, parameter :: MaxRunSuggestions = 512
    character(256) :: RunSuggestions(MaxRunSuggestions)
    !> Which heading each one goes under: 'sources' or 'flux'
    character(8) :: RunSuggestionSection(MaxRunSuggestions)
    integer :: RunSuggestionN = 0

    type(FCCsetupType) :: FCCsetup
    type(FileListType), allocatable :: FullFileList(:)
    type(FileListType), allocatable :: BinnedFileList(:)

    type(MeanSpectraType), allocatable :: MeanBinSpec(:, :)
    type(MeanSpectraType), allocatable :: dMeanBinSpec(:, :)
    type(MeanSpectraType), allocatable :: MeanBinCosp(:, :)

    type(MeanSpectraType) :: MeanStabilityCosp(ndkf, 2)
    type(MassParType) :: MassPar(GHGNumVar, 2)
    type(FluxType) :: Flux1
    type(FluxType) :: Flux2
    type(FluxType) :: Flux3
    !> Whole-run cache for the post-flux despiking pass (test_pfd). Grown
    !> one entry per period in ex_loop, guarded by EddyFlowProj%test_pfd so
    !> a run with the feature off never allocates it.
    type(PfdCacheEntryType), allocatable :: PfdCache(:)
    integer :: PfdCacheN = 0
    integer :: nCecPairs
    type(CECResolvedPairType) :: CecPairList(MaxNumCecPairs)
    type(CECFluxType) :: CECFlux(MaxNumCecPairs)
    type(FCCMetadataType) :: FCCMetadata

    !> tags of the setup ".ini" file for eccoce
    integer, parameter :: Nsn = 1005
    integer, parameter :: Nsc = 155
    logical            :: SNTagFound(Nsn)
    logical            :: SCTagFound(Nsc)
    type (Numerical)   :: SNTags(Nsn)
    type (Text)        :: SCTags(Nsc)
    !> BEGIN GENERATED FCC.SNTags - edit gen_project_tags.py, not this block
    data SNTags(1)%Label   / 'sa_nbins' / &
         SNTags(2)%Label   / 'sa_min_smpl' / &
         SNTags(3)%Label   / '' / &
         SNTags(4)%Label   / '' / &
         SNTags(5)%Label   / '' / &
         SNTags(6)%Label   / '' / &
         SNTags(7)%Label   / '' / &
         SNTags(8)%Label   / '' / &
         SNTags(9)%Label   / '' / &
         SNTags(10)%Label  / 'sa_auto_source' / &
         SNTags(11)%Label  / 'sa_min_co2' / &
         SNTags(12)%Label  / 'sa_min_ch4' / &
         SNTags(13)%Label  / 'sa_min_gas4' / &
         SNTags(14)%Label  / 'sa_min_le' / &
         SNTags(15)%Label  / 'sa_min_h' / &
         SNTags(16)%Label  / '' / &
         SNTags(17)%Label  / '' / &
         SNTags(18)%Label  / '' / &
         SNTags(19)%Label  / '' / &
         SNTags(20)%Label  / '' / &
         SNTags(21)%Label  / '' / &
         SNTags(22)%Label  / '' / &
         SNTags(23)%Label  / '' / &
         SNTags(24)%Label  / '' / &
         SNTags(25)%Label  / '' / &
         SNTags(26)%Label  / '' / &
         SNTags(27)%Label  / '' / &
         SNTags(28)%Label  / '' / &
         SNTags(29)%Label  / '' / &
         SNTags(30)%Label  / '' / &
         SNTags(31)%Label  / '' / &
         SNTags(32)%Label  / '' / &
         SNTags(33)%Label  / '' / &
         SNTags(34)%Label  / '' / &
         SNTags(35)%Label  / '' / &
         SNTags(36)%Label  / '' / &
         SNTags(37)%Label  / '' / &
         SNTags(38)%Label  / '' / &
         SNTags(39)%Label  / '' / &
         SNTags(40)%Label  / '' / &
         SNTags(41)%Label  / '' / &
         SNTags(42)%Label  / '' / &
         SNTags(43)%Label  / '' / &
         SNTags(44)%Label  / '' / &
         SNTags(45)%Label  / '' / &
         SNTags(46)%Label  / '' / &
         SNTags(47)%Label  / '' / &
         SNTags(48)%Label  / '' / &
         SNTags(49)%Label  / '' / &
         SNTags(50)%Label  / '' / &
         SNTags(51)%Label  / '' / &
         SNTags(52)%Label  / '' / &
         SNTags(53)%Label  / '' / &
         SNTags(54)%Label  / '' / &
         SNTags(55)%Label  / '' / &
         SNTags(56)%Label  / '' / &
         SNTags(57)%Label  / '' / &
         SNTags(58)%Label  / '' / &
         SNTags(59)%Label  / '' / &
         SNTags(60)%Label  / '' / &
         SNTags(61)%Label  / '' / &
         SNTags(62)%Label  / '' / &
         SNTags(63)%Label  / '' / &
         SNTags(64)%Label  / '' / &
         SNTags(65)%Label  / '' / &
         SNTags(66)%Label  / '' / &
         SNTags(67)%Label  / '' / &
         SNTags(68)%Label  / '' / &
         SNTags(69)%Label  / '' / &
         SNTags(70)%Label  / '' / &
         SNTags(71)%Label  / '' / &
         SNTags(72)%Label  / '' / &
         SNTags(73)%Label  / '' / &
         SNTags(74)%Label  / '' / &
         SNTags(75)%Label  / '' / &
         SNTags(76)%Label  / '' / &
         SNTags(77)%Label  / '' / &
         SNTags(78)%Label  / '' / &
         SNTags(79)%Label  / '' / &
         SNTags(80)%Label  / '' / &
         SNTags(81)%Label  / '' / &
         SNTags(82)%Label  / '' / &
         SNTags(83)%Label  / '' / &
         SNTags(84)%Label  / '' / &
         SNTags(85)%Label  / '' / &
         SNTags(86)%Label  / '' / &
         SNTags(87)%Label  / '' / &
         SNTags(88)%Label  / '' / &
         SNTags(89)%Label  / '' / &
         SNTags(90)%Label  / '' / &
         SNTags(91)%Label  / '' / &
         SNTags(92)%Label  / 'sa_min_un_ustar' / &
         SNTags(93)%Label  / '' / &
         SNTags(94)%Label  / '' / &
         SNTags(95)%Label  / '' / &
         SNTags(96)%Label  / 'sa_min_un_le' / &
         SNTags(97)%Label  / 'sa_min_un_h' / &
         SNTags(98)%Label  / 'sa_min_st_ustar' / &
         SNTags(99)%Label  / '' / &
         SNTags(100)%Label / '' / &
         SNTags(101)%Label / '' / &
         SNTags(102)%Label / 'sa_min_st_le' / &
         SNTags(103)%Label / 'sa_min_st_h' / &
         SNTags(104)%Label / 'sa_max_ustar' / &
         SNTags(105)%Label / '' / &
         SNTags(106)%Label / '' / &
         SNTags(107)%Label / '' / &
         SNTags(108)%Label / 'sa_max_le' / &
         SNTags(109)%Label / 'sa_max_h' / &
         SNTags(110)%Label / 'gas_1_sa_fmin' / &
         SNTags(111)%Label / 'gas_1_sa_fmax' / &
         SNTags(112)%Label / 'gas_1_sa_hfn_fmin' / &
         SNTags(113)%Label / 'gas_1_sa_min_st' / &
         SNTags(114)%Label / 'gas_1_sa_min_un' / &
         SNTags(115)%Label / 'gas_1_sa_max' / &
         SNTags(116)%Label / 'gas_1_sa_source' / &
         SNTags(117)%Label / 'gas_2_sa_fmin' / &
         SNTags(118)%Label / 'gas_2_sa_fmax' / &
         SNTags(119)%Label / 'gas_2_sa_hfn_fmin' / &
         SNTags(120)%Label / 'gas_2_sa_min_st' / &
         SNTags(121)%Label / 'gas_2_sa_min_un' / &
         SNTags(122)%Label / 'gas_2_sa_max' / &
         SNTags(123)%Label / 'gas_2_sa_source' / &
         SNTags(124)%Label / 'gas_3_sa_fmin' / &
         SNTags(125)%Label / 'gas_3_sa_fmax' / &
         SNTags(126)%Label / 'gas_3_sa_hfn_fmin' / &
         SNTags(127)%Label / 'gas_3_sa_min_st' / &
         SNTags(128)%Label / 'gas_3_sa_min_un' / &
         SNTags(129)%Label / 'gas_3_sa_max' / &
         SNTags(130)%Label / 'gas_3_sa_source' / &
         SNTags(131)%Label / 'gas_4_sa_fmin' / &
         SNTags(132)%Label / 'gas_4_sa_fmax' / &
         SNTags(133)%Label / 'gas_4_sa_hfn_fmin' / &
         SNTags(134)%Label / 'gas_4_sa_min_st' / &
         SNTags(135)%Label / 'gas_4_sa_min_un' / &
         SNTags(136)%Label / 'gas_4_sa_max' / &
         SNTags(137)%Label / 'gas_4_sa_source' / &
         SNTags(138)%Label / 'gas_5_sa_fmin' / &
         SNTags(139)%Label / 'gas_5_sa_fmax' / &
         SNTags(140)%Label / 'gas_5_sa_hfn_fmin' / &
         SNTags(141)%Label / 'gas_5_sa_min_st' / &
         SNTags(142)%Label / 'gas_5_sa_min_un' / &
         SNTags(143)%Label / 'gas_5_sa_max' / &
         SNTags(144)%Label / 'gas_5_sa_source' / &
         SNTags(145)%Label / 'gas_6_sa_fmin' / &
         SNTags(146)%Label / 'gas_6_sa_fmax' / &
         SNTags(147)%Label / 'gas_6_sa_hfn_fmin' / &
         SNTags(148)%Label / 'gas_6_sa_min_st' / &
         SNTags(149)%Label / 'gas_6_sa_min_un' / &
         SNTags(150)%Label / 'gas_6_sa_max' / &
         SNTags(151)%Label / 'gas_6_sa_source' / &
         SNTags(152)%Label / 'gas_7_sa_fmin' / &
         SNTags(153)%Label / 'gas_7_sa_fmax' / &
         SNTags(154)%Label / 'gas_7_sa_hfn_fmin' / &
         SNTags(155)%Label / 'gas_7_sa_min_st' / &
         SNTags(156)%Label / 'gas_7_sa_min_un' / &
         SNTags(157)%Label / 'gas_7_sa_max' / &
         SNTags(158)%Label / 'gas_7_sa_source' / &
         SNTags(159)%Label / 'gas_8_sa_fmin' / &
         SNTags(160)%Label / 'gas_8_sa_fmax' / &
         SNTags(161)%Label / 'gas_8_sa_hfn_fmin' / &
         SNTags(162)%Label / 'gas_8_sa_min_st' / &
         SNTags(163)%Label / 'gas_8_sa_min_un' / &
         SNTags(164)%Label / 'gas_8_sa_max' / &
         SNTags(165)%Label / 'gas_8_sa_source' / &
         SNTags(166)%Label / 'gas_9_sa_fmin' / &
         SNTags(167)%Label / 'gas_9_sa_fmax' / &
         SNTags(168)%Label / 'gas_9_sa_hfn_fmin' / &
         SNTags(169)%Label / 'gas_9_sa_min_st' / &
         SNTags(170)%Label / 'gas_9_sa_min_un' / &
         SNTags(171)%Label / 'gas_9_sa_max' / &
         SNTags(172)%Label / 'gas_9_sa_source' / &
         SNTags(173)%Label / 'gas_10_sa_fmin' / &
         SNTags(174)%Label / 'gas_10_sa_fmax' / &
         SNTags(175)%Label / 'gas_10_sa_hfn_fmin' / &
         SNTags(176)%Label / 'gas_10_sa_min_st' / &
         SNTags(177)%Label / 'gas_10_sa_min_un' / &
         SNTags(178)%Label / 'gas_10_sa_max' / &
         SNTags(179)%Label / 'gas_10_sa_source' / &
         SNTags(180)%Label / 'gas_11_sa_fmin' / &
         SNTags(181)%Label / 'gas_11_sa_fmax' / &
         SNTags(182)%Label / 'gas_11_sa_hfn_fmin' / &
         SNTags(183)%Label / 'gas_11_sa_min_st' / &
         SNTags(184)%Label / 'gas_11_sa_min_un' / &
         SNTags(185)%Label / 'gas_11_sa_max' / &
         SNTags(186)%Label / 'gas_11_sa_source' / &
         SNTags(187)%Label / 'gas_12_sa_fmin' / &
         SNTags(188)%Label / 'gas_12_sa_fmax' / &
         SNTags(189)%Label / 'gas_12_sa_hfn_fmin' / &
         SNTags(190)%Label / 'gas_12_sa_min_st' / &
         SNTags(191)%Label / 'gas_12_sa_min_un' / &
         SNTags(192)%Label / 'gas_12_sa_max' / &
         SNTags(193)%Label / 'gas_12_sa_source' / &
         SNTags(194)%Label / 'gas_13_sa_fmin' / &
         SNTags(195)%Label / 'gas_13_sa_fmax' / &
         SNTags(196)%Label / 'gas_13_sa_hfn_fmin' / &
         SNTags(197)%Label / 'gas_13_sa_min_st' / &
         SNTags(198)%Label / 'gas_13_sa_min_un' / &
         SNTags(199)%Label / 'gas_13_sa_max' / &
         SNTags(200)%Label / 'gas_13_sa_source' /
    data SNTags(201)%Label / 'gas_14_sa_fmin' / &
         SNTags(202)%Label / 'gas_14_sa_fmax' / &
         SNTags(203)%Label / 'gas_14_sa_hfn_fmin' / &
         SNTags(204)%Label / 'gas_14_sa_min_st' / &
         SNTags(205)%Label / 'gas_14_sa_min_un' / &
         SNTags(206)%Label / 'gas_14_sa_max' / &
         SNTags(207)%Label / 'gas_14_sa_source' / &
         SNTags(208)%Label / 'gas_15_sa_fmin' / &
         SNTags(209)%Label / 'gas_15_sa_fmax' / &
         SNTags(210)%Label / 'gas_15_sa_hfn_fmin' / &
         SNTags(211)%Label / 'gas_15_sa_min_st' / &
         SNTags(212)%Label / 'gas_15_sa_min_un' / &
         SNTags(213)%Label / 'gas_15_sa_max' / &
         SNTags(214)%Label / 'gas_15_sa_source' / &
         SNTags(215)%Label / 'gas_16_sa_fmin' / &
         SNTags(216)%Label / 'gas_16_sa_fmax' / &
         SNTags(217)%Label / 'gas_16_sa_hfn_fmin' / &
         SNTags(218)%Label / 'gas_16_sa_min_st' / &
         SNTags(219)%Label / 'gas_16_sa_min_un' / &
         SNTags(220)%Label / 'gas_16_sa_max' / &
         SNTags(221)%Label / 'gas_16_sa_source' / &
         SNTags(222)%Label / 'gas_17_sa_fmin' / &
         SNTags(223)%Label / 'gas_17_sa_fmax' / &
         SNTags(224)%Label / 'gas_17_sa_hfn_fmin' / &
         SNTags(225)%Label / 'gas_17_sa_min_st' / &
         SNTags(226)%Label / 'gas_17_sa_min_un' / &
         SNTags(227)%Label / 'gas_17_sa_max' / &
         SNTags(228)%Label / 'gas_17_sa_source' / &
         SNTags(229)%Label / 'gas_18_sa_fmin' / &
         SNTags(230)%Label / 'gas_18_sa_fmax' / &
         SNTags(231)%Label / 'gas_18_sa_hfn_fmin' / &
         SNTags(232)%Label / 'gas_18_sa_min_st' / &
         SNTags(233)%Label / 'gas_18_sa_min_un' / &
         SNTags(234)%Label / 'gas_18_sa_max' / &
         SNTags(235)%Label / 'gas_18_sa_source' / &
         SNTags(236)%Label / 'gas_19_sa_fmin' / &
         SNTags(237)%Label / 'gas_19_sa_fmax' / &
         SNTags(238)%Label / 'gas_19_sa_hfn_fmin' / &
         SNTags(239)%Label / 'gas_19_sa_min_st' / &
         SNTags(240)%Label / 'gas_19_sa_min_un' / &
         SNTags(241)%Label / 'gas_19_sa_max' / &
         SNTags(242)%Label / 'gas_19_sa_source' / &
         SNTags(243)%Label / 'gas_20_sa_fmin' / &
         SNTags(244)%Label / 'gas_20_sa_fmax' / &
         SNTags(245)%Label / 'gas_20_sa_hfn_fmin' / &
         SNTags(246)%Label / 'gas_20_sa_min_st' / &
         SNTags(247)%Label / 'gas_20_sa_min_un' / &
         SNTags(248)%Label / 'gas_20_sa_max' / &
         SNTags(249)%Label / 'gas_20_sa_source' / &
         SNTags(250)%Label / 'gas_21_sa_fmin' / &
         SNTags(251)%Label / 'gas_21_sa_fmax' / &
         SNTags(252)%Label / 'gas_21_sa_hfn_fmin' / &
         SNTags(253)%Label / 'gas_21_sa_min_st' / &
         SNTags(254)%Label / 'gas_21_sa_min_un' / &
         SNTags(255)%Label / 'gas_21_sa_max' / &
         SNTags(256)%Label / 'gas_21_sa_source' / &
         SNTags(257)%Label / 'gas_22_sa_fmin' / &
         SNTags(258)%Label / 'gas_22_sa_fmax' / &
         SNTags(259)%Label / 'gas_22_sa_hfn_fmin' / &
         SNTags(260)%Label / 'gas_22_sa_min_st' / &
         SNTags(261)%Label / 'gas_22_sa_min_un' / &
         SNTags(262)%Label / 'gas_22_sa_max' / &
         SNTags(263)%Label / 'gas_22_sa_source' / &
         SNTags(264)%Label / 'gas_23_sa_fmin' / &
         SNTags(265)%Label / 'gas_23_sa_fmax' / &
         SNTags(266)%Label / 'gas_23_sa_hfn_fmin' / &
         SNTags(267)%Label / 'gas_23_sa_min_st' / &
         SNTags(268)%Label / 'gas_23_sa_min_un' / &
         SNTags(269)%Label / 'gas_23_sa_max' / &
         SNTags(270)%Label / 'gas_23_sa_source' / &
         SNTags(271)%Label / 'gas_24_sa_fmin' / &
         SNTags(272)%Label / 'gas_24_sa_fmax' / &
         SNTags(273)%Label / 'gas_24_sa_hfn_fmin' / &
         SNTags(274)%Label / 'gas_24_sa_min_st' / &
         SNTags(275)%Label / 'gas_24_sa_min_un' / &
         SNTags(276)%Label / 'gas_24_sa_max' / &
         SNTags(277)%Label / 'gas_24_sa_source' / &
         SNTags(278)%Label / 'gas_25_sa_fmin' / &
         SNTags(279)%Label / 'gas_25_sa_fmax' / &
         SNTags(280)%Label / 'gas_25_sa_hfn_fmin' / &
         SNTags(281)%Label / 'gas_25_sa_min_st' / &
         SNTags(282)%Label / 'gas_25_sa_min_un' / &
         SNTags(283)%Label / 'gas_25_sa_max' / &
         SNTags(284)%Label / 'gas_25_sa_source' / &
         SNTags(285)%Label / 'gas_26_sa_fmin' / &
         SNTags(286)%Label / 'gas_26_sa_fmax' / &
         SNTags(287)%Label / 'gas_26_sa_hfn_fmin' / &
         SNTags(288)%Label / 'gas_26_sa_min_st' / &
         SNTags(289)%Label / 'gas_26_sa_min_un' / &
         SNTags(290)%Label / 'gas_26_sa_max' / &
         SNTags(291)%Label / 'gas_26_sa_source' / &
         SNTags(292)%Label / 'gas_27_sa_fmin' / &
         SNTags(293)%Label / 'gas_27_sa_fmax' / &
         SNTags(294)%Label / 'gas_27_sa_hfn_fmin' / &
         SNTags(295)%Label / 'gas_27_sa_min_st' / &
         SNTags(296)%Label / 'gas_27_sa_min_un' / &
         SNTags(297)%Label / 'gas_27_sa_max' / &
         SNTags(298)%Label / 'gas_27_sa_source' / &
         SNTags(299)%Label / 'gas_28_sa_fmin' / &
         SNTags(300)%Label / 'gas_28_sa_fmax' / &
         SNTags(301)%Label / 'gas_28_sa_hfn_fmin' / &
         SNTags(302)%Label / 'gas_28_sa_min_st' / &
         SNTags(303)%Label / 'gas_28_sa_min_un' / &
         SNTags(304)%Label / 'gas_28_sa_max' / &
         SNTags(305)%Label / 'gas_28_sa_source' / &
         SNTags(306)%Label / 'gas_29_sa_fmin' / &
         SNTags(307)%Label / 'gas_29_sa_fmax' / &
         SNTags(308)%Label / 'gas_29_sa_hfn_fmin' / &
         SNTags(309)%Label / 'gas_29_sa_min_st' / &
         SNTags(310)%Label / 'gas_29_sa_min_un' / &
         SNTags(311)%Label / 'gas_29_sa_max' / &
         SNTags(312)%Label / 'gas_29_sa_source' / &
         SNTags(313)%Label / 'gas_30_sa_fmin' / &
         SNTags(314)%Label / 'gas_30_sa_fmax' / &
         SNTags(315)%Label / 'gas_30_sa_hfn_fmin' / &
         SNTags(316)%Label / 'gas_30_sa_min_st' / &
         SNTags(317)%Label / 'gas_30_sa_min_un' / &
         SNTags(318)%Label / 'gas_30_sa_max' / &
         SNTags(319)%Label / 'gas_30_sa_source' / &
         SNTags(320)%Label / 'gas_31_sa_fmin' / &
         SNTags(321)%Label / 'gas_31_sa_fmax' / &
         SNTags(322)%Label / 'gas_31_sa_hfn_fmin' / &
         SNTags(323)%Label / 'gas_31_sa_min_st' / &
         SNTags(324)%Label / 'gas_31_sa_min_un' / &
         SNTags(325)%Label / 'gas_31_sa_max' / &
         SNTags(326)%Label / 'gas_31_sa_source' / &
         SNTags(327)%Label / 'gas_32_sa_fmin' / &
         SNTags(328)%Label / 'gas_32_sa_fmax' / &
         SNTags(329)%Label / 'gas_32_sa_hfn_fmin' / &
         SNTags(330)%Label / 'gas_32_sa_min_st' / &
         SNTags(331)%Label / 'gas_32_sa_min_un' / &
         SNTags(332)%Label / 'gas_32_sa_max' / &
         SNTags(333)%Label / 'gas_32_sa_source' / &
         SNTags(334)%Label / 'gas_33_sa_fmin' / &
         SNTags(335)%Label / 'gas_33_sa_fmax' / &
         SNTags(336)%Label / 'gas_33_sa_hfn_fmin' / &
         SNTags(337)%Label / 'gas_33_sa_min_st' / &
         SNTags(338)%Label / 'gas_33_sa_min_un' / &
         SNTags(339)%Label / 'gas_33_sa_max' / &
         SNTags(340)%Label / 'gas_33_sa_source' / &
         SNTags(341)%Label / 'gas_34_sa_fmin' / &
         SNTags(342)%Label / 'gas_34_sa_fmax' / &
         SNTags(343)%Label / 'gas_34_sa_hfn_fmin' / &
         SNTags(344)%Label / 'gas_34_sa_min_st' / &
         SNTags(345)%Label / 'gas_34_sa_min_un' / &
         SNTags(346)%Label / 'gas_34_sa_max' / &
         SNTags(347)%Label / 'gas_34_sa_source' / &
         SNTags(348)%Label / 'gas_35_sa_fmin' / &
         SNTags(349)%Label / 'gas_35_sa_fmax' / &
         SNTags(350)%Label / 'gas_35_sa_hfn_fmin' / &
         SNTags(351)%Label / 'gas_35_sa_min_st' / &
         SNTags(352)%Label / 'gas_35_sa_min_un' / &
         SNTags(353)%Label / 'gas_35_sa_max' / &
         SNTags(354)%Label / 'gas_35_sa_source' / &
         SNTags(355)%Label / 'gas_36_sa_fmin' / &
         SNTags(356)%Label / 'gas_36_sa_fmax' / &
         SNTags(357)%Label / 'gas_36_sa_hfn_fmin' / &
         SNTags(358)%Label / 'gas_36_sa_min_st' / &
         SNTags(359)%Label / 'gas_36_sa_min_un' / &
         SNTags(360)%Label / 'gas_36_sa_max' / &
         SNTags(361)%Label / 'gas_36_sa_source' / &
         SNTags(362)%Label / 'gas_37_sa_fmin' / &
         SNTags(363)%Label / 'gas_37_sa_fmax' / &
         SNTags(364)%Label / 'gas_37_sa_hfn_fmin' / &
         SNTags(365)%Label / 'gas_37_sa_min_st' / &
         SNTags(366)%Label / 'gas_37_sa_min_un' / &
         SNTags(367)%Label / 'gas_37_sa_max' / &
         SNTags(368)%Label / 'gas_37_sa_source' / &
         SNTags(369)%Label / 'gas_38_sa_fmin' / &
         SNTags(370)%Label / 'gas_38_sa_fmax' / &
         SNTags(371)%Label / 'gas_38_sa_hfn_fmin' / &
         SNTags(372)%Label / 'gas_38_sa_min_st' / &
         SNTags(373)%Label / 'gas_38_sa_min_un' / &
         SNTags(374)%Label / 'gas_38_sa_max' / &
         SNTags(375)%Label / 'gas_38_sa_source' / &
         SNTags(376)%Label / 'gas_39_sa_fmin' / &
         SNTags(377)%Label / 'gas_39_sa_fmax' / &
         SNTags(378)%Label / 'gas_39_sa_hfn_fmin' / &
         SNTags(379)%Label / 'gas_39_sa_min_st' / &
         SNTags(380)%Label / 'gas_39_sa_min_un' / &
         SNTags(381)%Label / 'gas_39_sa_max' / &
         SNTags(382)%Label / 'gas_39_sa_source' / &
         SNTags(383)%Label / 'gas_40_sa_fmin' / &
         SNTags(384)%Label / 'gas_40_sa_fmax' / &
         SNTags(385)%Label / 'gas_40_sa_hfn_fmin' / &
         SNTags(386)%Label / 'gas_40_sa_min_st' / &
         SNTags(387)%Label / 'gas_40_sa_min_un' / &
         SNTags(388)%Label / 'gas_40_sa_max' / &
         SNTags(389)%Label / 'gas_40_sa_source' / &
         SNTags(390)%Label / 'gas_41_sa_fmin' / &
         SNTags(391)%Label / 'gas_41_sa_fmax' / &
         SNTags(392)%Label / 'gas_41_sa_hfn_fmin' / &
         SNTags(393)%Label / 'gas_41_sa_min_st' / &
         SNTags(394)%Label / 'gas_41_sa_min_un' / &
         SNTags(395)%Label / 'gas_41_sa_max' / &
         SNTags(396)%Label / 'gas_41_sa_source' / &
         SNTags(397)%Label / 'gas_42_sa_fmin' / &
         SNTags(398)%Label / 'gas_42_sa_fmax' / &
         SNTags(399)%Label / 'gas_42_sa_hfn_fmin' / &
         SNTags(400)%Label / 'gas_42_sa_min_st' /
    data SNTags(401)%Label / 'gas_42_sa_min_un' / &
         SNTags(402)%Label / 'gas_42_sa_max' / &
         SNTags(403)%Label / 'gas_42_sa_source' / &
         SNTags(404)%Label / 'gas_43_sa_fmin' / &
         SNTags(405)%Label / 'gas_43_sa_fmax' / &
         SNTags(406)%Label / 'gas_43_sa_hfn_fmin' / &
         SNTags(407)%Label / 'gas_43_sa_min_st' / &
         SNTags(408)%Label / 'gas_43_sa_min_un' / &
         SNTags(409)%Label / 'gas_43_sa_max' / &
         SNTags(410)%Label / 'gas_43_sa_source' / &
         SNTags(411)%Label / 'gas_44_sa_fmin' / &
         SNTags(412)%Label / 'gas_44_sa_fmax' / &
         SNTags(413)%Label / 'gas_44_sa_hfn_fmin' / &
         SNTags(414)%Label / 'gas_44_sa_min_st' / &
         SNTags(415)%Label / 'gas_44_sa_min_un' / &
         SNTags(416)%Label / 'gas_44_sa_max' / &
         SNTags(417)%Label / 'gas_44_sa_source' / &
         SNTags(418)%Label / 'gas_45_sa_fmin' / &
         SNTags(419)%Label / 'gas_45_sa_fmax' / &
         SNTags(420)%Label / 'gas_45_sa_hfn_fmin' / &
         SNTags(421)%Label / 'gas_45_sa_min_st' / &
         SNTags(422)%Label / 'gas_45_sa_min_un' / &
         SNTags(423)%Label / 'gas_45_sa_max' / &
         SNTags(424)%Label / 'gas_45_sa_source' / &
         SNTags(425)%Label / 'gas_46_sa_fmin' / &
         SNTags(426)%Label / 'gas_46_sa_fmax' / &
         SNTags(427)%Label / 'gas_46_sa_hfn_fmin' / &
         SNTags(428)%Label / 'gas_46_sa_min_st' / &
         SNTags(429)%Label / 'gas_46_sa_min_un' / &
         SNTags(430)%Label / 'gas_46_sa_max' / &
         SNTags(431)%Label / 'gas_46_sa_source' / &
         SNTags(432)%Label / 'gas_47_sa_fmin' / &
         SNTags(433)%Label / 'gas_47_sa_fmax' / &
         SNTags(434)%Label / 'gas_47_sa_hfn_fmin' / &
         SNTags(435)%Label / 'gas_47_sa_min_st' / &
         SNTags(436)%Label / 'gas_47_sa_min_un' / &
         SNTags(437)%Label / 'gas_47_sa_max' / &
         SNTags(438)%Label / 'gas_47_sa_source' / &
         SNTags(439)%Label / 'gas_48_sa_fmin' / &
         SNTags(440)%Label / 'gas_48_sa_fmax' / &
         SNTags(441)%Label / 'gas_48_sa_hfn_fmin' / &
         SNTags(442)%Label / 'gas_48_sa_min_st' / &
         SNTags(443)%Label / 'gas_48_sa_min_un' / &
         SNTags(444)%Label / 'gas_48_sa_max' / &
         SNTags(445)%Label / 'gas_48_sa_source' / &
         SNTags(446)%Label / 'gas_49_sa_fmin' / &
         SNTags(447)%Label / 'gas_49_sa_fmax' / &
         SNTags(448)%Label / 'gas_49_sa_hfn_fmin' / &
         SNTags(449)%Label / 'gas_49_sa_min_st' / &
         SNTags(450)%Label / 'gas_49_sa_min_un' / &
         SNTags(451)%Label / 'gas_49_sa_max' / &
         SNTags(452)%Label / 'gas_49_sa_source' / &
         SNTags(453)%Label / 'gas_50_sa_fmin' / &
         SNTags(454)%Label / 'gas_50_sa_fmax' / &
         SNTags(455)%Label / 'gas_50_sa_hfn_fmin' / &
         SNTags(456)%Label / 'gas_50_sa_min_st' / &
         SNTags(457)%Label / 'gas_50_sa_min_un' / &
         SNTags(458)%Label / 'gas_50_sa_max' / &
         SNTags(459)%Label / 'gas_50_sa_source' / &
         SNTags(460)%Label / 'gas_51_sa_fmin' / &
         SNTags(461)%Label / 'gas_51_sa_fmax' / &
         SNTags(462)%Label / 'gas_51_sa_hfn_fmin' / &
         SNTags(463)%Label / 'gas_51_sa_min_st' / &
         SNTags(464)%Label / 'gas_51_sa_min_un' / &
         SNTags(465)%Label / 'gas_51_sa_max' / &
         SNTags(466)%Label / 'gas_51_sa_source' / &
         SNTags(467)%Label / 'gas_52_sa_fmin' / &
         SNTags(468)%Label / 'gas_52_sa_fmax' / &
         SNTags(469)%Label / 'gas_52_sa_hfn_fmin' / &
         SNTags(470)%Label / 'gas_52_sa_min_st' / &
         SNTags(471)%Label / 'gas_52_sa_min_un' / &
         SNTags(472)%Label / 'gas_52_sa_max' / &
         SNTags(473)%Label / 'gas_52_sa_source' / &
         SNTags(474)%Label / 'gas_53_sa_fmin' / &
         SNTags(475)%Label / 'gas_53_sa_fmax' / &
         SNTags(476)%Label / 'gas_53_sa_hfn_fmin' / &
         SNTags(477)%Label / 'gas_53_sa_min_st' / &
         SNTags(478)%Label / 'gas_53_sa_min_un' / &
         SNTags(479)%Label / 'gas_53_sa_max' / &
         SNTags(480)%Label / 'gas_53_sa_source' / &
         SNTags(481)%Label / 'gas_54_sa_fmin' / &
         SNTags(482)%Label / 'gas_54_sa_fmax' / &
         SNTags(483)%Label / 'gas_54_sa_hfn_fmin' / &
         SNTags(484)%Label / 'gas_54_sa_min_st' / &
         SNTags(485)%Label / 'gas_54_sa_min_un' / &
         SNTags(486)%Label / 'gas_54_sa_max' / &
         SNTags(487)%Label / 'gas_54_sa_source' / &
         SNTags(488)%Label / 'gas_55_sa_fmin' / &
         SNTags(489)%Label / 'gas_55_sa_fmax' / &
         SNTags(490)%Label / 'gas_55_sa_hfn_fmin' / &
         SNTags(491)%Label / 'gas_55_sa_min_st' / &
         SNTags(492)%Label / 'gas_55_sa_min_un' / &
         SNTags(493)%Label / 'gas_55_sa_max' / &
         SNTags(494)%Label / 'gas_55_sa_source' / &
         SNTags(495)%Label / 'gas_56_sa_fmin' / &
         SNTags(496)%Label / 'gas_56_sa_fmax' / &
         SNTags(497)%Label / 'gas_56_sa_hfn_fmin' / &
         SNTags(498)%Label / 'gas_56_sa_min_st' / &
         SNTags(499)%Label / 'gas_56_sa_min_un' / &
         SNTags(500)%Label / 'gas_56_sa_max' / &
         SNTags(501)%Label / 'gas_56_sa_source' / &
         SNTags(502)%Label / 'gas_57_sa_fmin' / &
         SNTags(503)%Label / 'gas_57_sa_fmax' / &
         SNTags(504)%Label / 'gas_57_sa_hfn_fmin' / &
         SNTags(505)%Label / 'gas_57_sa_min_st' / &
         SNTags(506)%Label / 'gas_57_sa_min_un' / &
         SNTags(507)%Label / 'gas_57_sa_max' / &
         SNTags(508)%Label / 'gas_57_sa_source' / &
         SNTags(509)%Label / 'gas_58_sa_fmin' / &
         SNTags(510)%Label / 'gas_58_sa_fmax' / &
         SNTags(511)%Label / 'gas_58_sa_hfn_fmin' / &
         SNTags(512)%Label / 'gas_58_sa_min_st' / &
         SNTags(513)%Label / 'gas_58_sa_min_un' / &
         SNTags(514)%Label / 'gas_58_sa_max' / &
         SNTags(515)%Label / 'gas_58_sa_source' / &
         SNTags(516)%Label / 'gas_59_sa_fmin' / &
         SNTags(517)%Label / 'gas_59_sa_fmax' / &
         SNTags(518)%Label / 'gas_59_sa_hfn_fmin' / &
         SNTags(519)%Label / 'gas_59_sa_min_st' / &
         SNTags(520)%Label / 'gas_59_sa_min_un' / &
         SNTags(521)%Label / 'gas_59_sa_max' / &
         SNTags(522)%Label / 'gas_59_sa_source' / &
         SNTags(523)%Label / 'gas_60_sa_fmin' / &
         SNTags(524)%Label / 'gas_60_sa_fmax' / &
         SNTags(525)%Label / 'gas_60_sa_hfn_fmin' / &
         SNTags(526)%Label / 'gas_60_sa_min_st' / &
         SNTags(527)%Label / 'gas_60_sa_min_un' / &
         SNTags(528)%Label / 'gas_60_sa_max' / &
         SNTags(529)%Label / 'gas_60_sa_source' / &
         SNTags(530)%Label / 'gas_61_sa_fmin' / &
         SNTags(531)%Label / 'gas_61_sa_fmax' / &
         SNTags(532)%Label / 'gas_61_sa_hfn_fmin' / &
         SNTags(533)%Label / 'gas_61_sa_min_st' / &
         SNTags(534)%Label / 'gas_61_sa_min_un' / &
         SNTags(535)%Label / 'gas_61_sa_max' / &
         SNTags(536)%Label / 'gas_61_sa_source' / &
         SNTags(537)%Label / 'gas_62_sa_fmin' / &
         SNTags(538)%Label / 'gas_62_sa_fmax' / &
         SNTags(539)%Label / 'gas_62_sa_hfn_fmin' / &
         SNTags(540)%Label / 'gas_62_sa_min_st' / &
         SNTags(541)%Label / 'gas_62_sa_min_un' / &
         SNTags(542)%Label / 'gas_62_sa_max' / &
         SNTags(543)%Label / 'gas_62_sa_source' / &
         SNTags(544)%Label / 'gas_63_sa_fmin' / &
         SNTags(545)%Label / 'gas_63_sa_fmax' / &
         SNTags(546)%Label / 'gas_63_sa_hfn_fmin' / &
         SNTags(547)%Label / 'gas_63_sa_min_st' / &
         SNTags(548)%Label / 'gas_63_sa_min_un' / &
         SNTags(549)%Label / 'gas_63_sa_max' / &
         SNTags(550)%Label / 'gas_63_sa_source' / &
         SNTags(551)%Label / 'gas_64_sa_fmin' / &
         SNTags(552)%Label / 'gas_64_sa_fmax' / &
         SNTags(553)%Label / 'gas_64_sa_hfn_fmin' / &
         SNTags(554)%Label / 'gas_64_sa_min_st' / &
         SNTags(555)%Label / 'gas_64_sa_min_un' / &
         SNTags(556)%Label / 'gas_64_sa_max' / &
         SNTags(557)%Label / 'gas_64_sa_source' / &
         SNTags(558)%Label / 'gas_65_sa_fmin' / &
         SNTags(559)%Label / 'gas_65_sa_fmax' / &
         SNTags(560)%Label / 'gas_65_sa_hfn_fmin' / &
         SNTags(561)%Label / 'gas_65_sa_min_st' / &
         SNTags(562)%Label / 'gas_65_sa_min_un' / &
         SNTags(563)%Label / 'gas_65_sa_max' / &
         SNTags(564)%Label / 'gas_65_sa_source' / &
         SNTags(565)%Label / 'gas_66_sa_fmin' / &
         SNTags(566)%Label / 'gas_66_sa_fmax' / &
         SNTags(567)%Label / 'gas_66_sa_hfn_fmin' / &
         SNTags(568)%Label / 'gas_66_sa_min_st' / &
         SNTags(569)%Label / 'gas_66_sa_min_un' / &
         SNTags(570)%Label / 'gas_66_sa_max' / &
         SNTags(571)%Label / 'gas_66_sa_source' / &
         SNTags(572)%Label / 'gas_67_sa_fmin' / &
         SNTags(573)%Label / 'gas_67_sa_fmax' / &
         SNTags(574)%Label / 'gas_67_sa_hfn_fmin' / &
         SNTags(575)%Label / 'gas_67_sa_min_st' / &
         SNTags(576)%Label / 'gas_67_sa_min_un' / &
         SNTags(577)%Label / 'gas_67_sa_max' / &
         SNTags(578)%Label / 'gas_67_sa_source' / &
         SNTags(579)%Label / 'gas_68_sa_fmin' / &
         SNTags(580)%Label / 'gas_68_sa_fmax' / &
         SNTags(581)%Label / 'gas_68_sa_hfn_fmin' / &
         SNTags(582)%Label / 'gas_68_sa_min_st' / &
         SNTags(583)%Label / 'gas_68_sa_min_un' / &
         SNTags(584)%Label / 'gas_68_sa_max' / &
         SNTags(585)%Label / 'gas_68_sa_source' / &
         SNTags(586)%Label / 'gas_69_sa_fmin' / &
         SNTags(587)%Label / 'gas_69_sa_fmax' / &
         SNTags(588)%Label / 'gas_69_sa_hfn_fmin' / &
         SNTags(589)%Label / 'gas_69_sa_min_st' / &
         SNTags(590)%Label / 'gas_69_sa_min_un' / &
         SNTags(591)%Label / 'gas_69_sa_max' / &
         SNTags(592)%Label / 'gas_69_sa_source' / &
         SNTags(593)%Label / 'gas_70_sa_fmin' / &
         SNTags(594)%Label / 'gas_70_sa_fmax' / &
         SNTags(595)%Label / 'gas_70_sa_hfn_fmin' / &
         SNTags(596)%Label / 'gas_70_sa_min_st' / &
         SNTags(597)%Label / 'gas_70_sa_min_un' / &
         SNTags(598)%Label / 'gas_70_sa_max' / &
         SNTags(599)%Label / 'gas_70_sa_source' / &
         SNTags(600)%Label / 'gas_71_sa_fmin' /
    data SNTags(601)%Label / 'gas_71_sa_fmax' / &
         SNTags(602)%Label / 'gas_71_sa_hfn_fmin' / &
         SNTags(603)%Label / 'gas_71_sa_min_st' / &
         SNTags(604)%Label / 'gas_71_sa_min_un' / &
         SNTags(605)%Label / 'gas_71_sa_max' / &
         SNTags(606)%Label / 'gas_71_sa_source' / &
         SNTags(607)%Label / 'gas_72_sa_fmin' / &
         SNTags(608)%Label / 'gas_72_sa_fmax' / &
         SNTags(609)%Label / 'gas_72_sa_hfn_fmin' / &
         SNTags(610)%Label / 'gas_72_sa_min_st' / &
         SNTags(611)%Label / 'gas_72_sa_min_un' / &
         SNTags(612)%Label / 'gas_72_sa_max' / &
         SNTags(613)%Label / 'gas_72_sa_source' / &
         SNTags(614)%Label / 'gas_73_sa_fmin' / &
         SNTags(615)%Label / 'gas_73_sa_fmax' / &
         SNTags(616)%Label / 'gas_73_sa_hfn_fmin' / &
         SNTags(617)%Label / 'gas_73_sa_min_st' / &
         SNTags(618)%Label / 'gas_73_sa_min_un' / &
         SNTags(619)%Label / 'gas_73_sa_max' / &
         SNTags(620)%Label / 'gas_73_sa_source' / &
         SNTags(621)%Label / 'gas_74_sa_fmin' / &
         SNTags(622)%Label / 'gas_74_sa_fmax' / &
         SNTags(623)%Label / 'gas_74_sa_hfn_fmin' / &
         SNTags(624)%Label / 'gas_74_sa_min_st' / &
         SNTags(625)%Label / 'gas_74_sa_min_un' / &
         SNTags(626)%Label / 'gas_74_sa_max' / &
         SNTags(627)%Label / 'gas_74_sa_source' / &
         SNTags(628)%Label / 'gas_75_sa_fmin' / &
         SNTags(629)%Label / 'gas_75_sa_fmax' / &
         SNTags(630)%Label / 'gas_75_sa_hfn_fmin' / &
         SNTags(631)%Label / 'gas_75_sa_min_st' / &
         SNTags(632)%Label / 'gas_75_sa_min_un' / &
         SNTags(633)%Label / 'gas_75_sa_max' / &
         SNTags(634)%Label / 'gas_75_sa_source' / &
         SNTags(635)%Label / 'gas_76_sa_fmin' / &
         SNTags(636)%Label / 'gas_76_sa_fmax' / &
         SNTags(637)%Label / 'gas_76_sa_hfn_fmin' / &
         SNTags(638)%Label / 'gas_76_sa_min_st' / &
         SNTags(639)%Label / 'gas_76_sa_min_un' / &
         SNTags(640)%Label / 'gas_76_sa_max' / &
         SNTags(641)%Label / 'gas_76_sa_source' / &
         SNTags(642)%Label / 'gas_77_sa_fmin' / &
         SNTags(643)%Label / 'gas_77_sa_fmax' / &
         SNTags(644)%Label / 'gas_77_sa_hfn_fmin' / &
         SNTags(645)%Label / 'gas_77_sa_min_st' / &
         SNTags(646)%Label / 'gas_77_sa_min_un' / &
         SNTags(647)%Label / 'gas_77_sa_max' / &
         SNTags(648)%Label / 'gas_77_sa_source' / &
         SNTags(649)%Label / 'gas_78_sa_fmin' / &
         SNTags(650)%Label / 'gas_78_sa_fmax' / &
         SNTags(651)%Label / 'gas_78_sa_hfn_fmin' / &
         SNTags(652)%Label / 'gas_78_sa_min_st' / &
         SNTags(653)%Label / 'gas_78_sa_min_un' / &
         SNTags(654)%Label / 'gas_78_sa_max' / &
         SNTags(655)%Label / 'gas_78_sa_source' / &
         SNTags(656)%Label / 'gas_79_sa_fmin' / &
         SNTags(657)%Label / 'gas_79_sa_fmax' / &
         SNTags(658)%Label / 'gas_79_sa_hfn_fmin' / &
         SNTags(659)%Label / 'gas_79_sa_min_st' / &
         SNTags(660)%Label / 'gas_79_sa_min_un' / &
         SNTags(661)%Label / 'gas_79_sa_max' / &
         SNTags(662)%Label / 'gas_79_sa_source' / &
         SNTags(663)%Label / 'gas_80_sa_fmin' / &
         SNTags(664)%Label / 'gas_80_sa_fmax' / &
         SNTags(665)%Label / 'gas_80_sa_hfn_fmin' / &
         SNTags(666)%Label / 'gas_80_sa_min_st' / &
         SNTags(667)%Label / 'gas_80_sa_min_un' / &
         SNTags(668)%Label / 'gas_80_sa_max' / &
         SNTags(669)%Label / 'gas_80_sa_source' / &
         SNTags(670)%Label / 'gas_81_sa_fmin' / &
         SNTags(671)%Label / 'gas_81_sa_fmax' / &
         SNTags(672)%Label / 'gas_81_sa_hfn_fmin' / &
         SNTags(673)%Label / 'gas_81_sa_min_st' / &
         SNTags(674)%Label / 'gas_81_sa_min_un' / &
         SNTags(675)%Label / 'gas_81_sa_max' / &
         SNTags(676)%Label / 'gas_81_sa_source' / &
         SNTags(677)%Label / 'gas_82_sa_fmin' / &
         SNTags(678)%Label / 'gas_82_sa_fmax' / &
         SNTags(679)%Label / 'gas_82_sa_hfn_fmin' / &
         SNTags(680)%Label / 'gas_82_sa_min_st' / &
         SNTags(681)%Label / 'gas_82_sa_min_un' / &
         SNTags(682)%Label / 'gas_82_sa_max' / &
         SNTags(683)%Label / 'gas_82_sa_source' / &
         SNTags(684)%Label / 'gas_83_sa_fmin' / &
         SNTags(685)%Label / 'gas_83_sa_fmax' / &
         SNTags(686)%Label / 'gas_83_sa_hfn_fmin' / &
         SNTags(687)%Label / 'gas_83_sa_min_st' / &
         SNTags(688)%Label / 'gas_83_sa_min_un' / &
         SNTags(689)%Label / 'gas_83_sa_max' / &
         SNTags(690)%Label / 'gas_83_sa_source' / &
         SNTags(691)%Label / 'gas_84_sa_fmin' / &
         SNTags(692)%Label / 'gas_84_sa_fmax' / &
         SNTags(693)%Label / 'gas_84_sa_hfn_fmin' / &
         SNTags(694)%Label / 'gas_84_sa_min_st' / &
         SNTags(695)%Label / 'gas_84_sa_min_un' / &
         SNTags(696)%Label / 'gas_84_sa_max' / &
         SNTags(697)%Label / 'gas_84_sa_source' / &
         SNTags(698)%Label / 'gas_85_sa_fmin' / &
         SNTags(699)%Label / 'gas_85_sa_fmax' / &
         SNTags(700)%Label / 'gas_85_sa_hfn_fmin' / &
         SNTags(701)%Label / 'gas_85_sa_min_st' / &
         SNTags(702)%Label / 'gas_85_sa_min_un' / &
         SNTags(703)%Label / 'gas_85_sa_max' / &
         SNTags(704)%Label / 'gas_85_sa_source' / &
         SNTags(705)%Label / 'gas_86_sa_fmin' / &
         SNTags(706)%Label / 'gas_86_sa_fmax' / &
         SNTags(707)%Label / 'gas_86_sa_hfn_fmin' / &
         SNTags(708)%Label / 'gas_86_sa_min_st' / &
         SNTags(709)%Label / 'gas_86_sa_min_un' / &
         SNTags(710)%Label / 'gas_86_sa_max' / &
         SNTags(711)%Label / 'gas_86_sa_source' / &
         SNTags(712)%Label / 'gas_87_sa_fmin' / &
         SNTags(713)%Label / 'gas_87_sa_fmax' / &
         SNTags(714)%Label / 'gas_87_sa_hfn_fmin' / &
         SNTags(715)%Label / 'gas_87_sa_min_st' / &
         SNTags(716)%Label / 'gas_87_sa_min_un' / &
         SNTags(717)%Label / 'gas_87_sa_max' / &
         SNTags(718)%Label / 'gas_87_sa_source' / &
         SNTags(719)%Label / 'gas_88_sa_fmin' / &
         SNTags(720)%Label / 'gas_88_sa_fmax' / &
         SNTags(721)%Label / 'gas_88_sa_hfn_fmin' / &
         SNTags(722)%Label / 'gas_88_sa_min_st' / &
         SNTags(723)%Label / 'gas_88_sa_min_un' / &
         SNTags(724)%Label / 'gas_88_sa_max' / &
         SNTags(725)%Label / 'gas_88_sa_source' / &
         SNTags(726)%Label / 'gas_89_sa_fmin' / &
         SNTags(727)%Label / 'gas_89_sa_fmax' / &
         SNTags(728)%Label / 'gas_89_sa_hfn_fmin' / &
         SNTags(729)%Label / 'gas_89_sa_min_st' / &
         SNTags(730)%Label / 'gas_89_sa_min_un' / &
         SNTags(731)%Label / 'gas_89_sa_max' / &
         SNTags(732)%Label / 'gas_89_sa_source' / &
         SNTags(733)%Label / 'gas_90_sa_fmin' / &
         SNTags(734)%Label / 'gas_90_sa_fmax' / &
         SNTags(735)%Label / 'gas_90_sa_hfn_fmin' / &
         SNTags(736)%Label / 'gas_90_sa_min_st' / &
         SNTags(737)%Label / 'gas_90_sa_min_un' / &
         SNTags(738)%Label / 'gas_90_sa_max' / &
         SNTags(739)%Label / 'gas_90_sa_source' / &
         SNTags(740)%Label / 'gas_91_sa_fmin' / &
         SNTags(741)%Label / 'gas_91_sa_fmax' / &
         SNTags(742)%Label / 'gas_91_sa_hfn_fmin' / &
         SNTags(743)%Label / 'gas_91_sa_min_st' / &
         SNTags(744)%Label / 'gas_91_sa_min_un' / &
         SNTags(745)%Label / 'gas_91_sa_max' / &
         SNTags(746)%Label / 'gas_91_sa_source' / &
         SNTags(747)%Label / 'gas_92_sa_fmin' / &
         SNTags(748)%Label / 'gas_92_sa_fmax' / &
         SNTags(749)%Label / 'gas_92_sa_hfn_fmin' / &
         SNTags(750)%Label / 'gas_92_sa_min_st' / &
         SNTags(751)%Label / 'gas_92_sa_min_un' / &
         SNTags(752)%Label / 'gas_92_sa_max' / &
         SNTags(753)%Label / 'gas_92_sa_source' / &
         SNTags(754)%Label / 'gas_93_sa_fmin' / &
         SNTags(755)%Label / 'gas_93_sa_fmax' / &
         SNTags(756)%Label / 'gas_93_sa_hfn_fmin' / &
         SNTags(757)%Label / 'gas_93_sa_min_st' / &
         SNTags(758)%Label / 'gas_93_sa_min_un' / &
         SNTags(759)%Label / 'gas_93_sa_max' / &
         SNTags(760)%Label / 'gas_93_sa_source' / &
         SNTags(761)%Label / 'gas_94_sa_fmin' / &
         SNTags(762)%Label / 'gas_94_sa_fmax' / &
         SNTags(763)%Label / 'gas_94_sa_hfn_fmin' / &
         SNTags(764)%Label / 'gas_94_sa_min_st' / &
         SNTags(765)%Label / 'gas_94_sa_min_un' / &
         SNTags(766)%Label / 'gas_94_sa_max' / &
         SNTags(767)%Label / 'gas_94_sa_source' / &
         SNTags(768)%Label / 'gas_95_sa_fmin' / &
         SNTags(769)%Label / 'gas_95_sa_fmax' / &
         SNTags(770)%Label / 'gas_95_sa_hfn_fmin' / &
         SNTags(771)%Label / 'gas_95_sa_min_st' / &
         SNTags(772)%Label / 'gas_95_sa_min_un' / &
         SNTags(773)%Label / 'gas_95_sa_max' / &
         SNTags(774)%Label / 'gas_95_sa_source' / &
         SNTags(775)%Label / 'gas_96_sa_fmin' / &
         SNTags(776)%Label / 'gas_96_sa_fmax' / &
         SNTags(777)%Label / 'gas_96_sa_hfn_fmin' / &
         SNTags(778)%Label / 'gas_96_sa_min_st' / &
         SNTags(779)%Label / 'gas_96_sa_min_un' / &
         SNTags(780)%Label / 'gas_96_sa_max' / &
         SNTags(781)%Label / 'gas_96_sa_source' / &
         SNTags(782)%Label / 'gas_97_sa_fmin' / &
         SNTags(783)%Label / 'gas_97_sa_fmax' / &
         SNTags(784)%Label / 'gas_97_sa_hfn_fmin' / &
         SNTags(785)%Label / 'gas_97_sa_min_st' / &
         SNTags(786)%Label / 'gas_97_sa_min_un' / &
         SNTags(787)%Label / 'gas_97_sa_max' / &
         SNTags(788)%Label / 'gas_97_sa_source' / &
         SNTags(789)%Label / 'gas_98_sa_fmin' / &
         SNTags(790)%Label / 'gas_98_sa_fmax' / &
         SNTags(791)%Label / 'gas_98_sa_hfn_fmin' / &
         SNTags(792)%Label / 'gas_98_sa_min_st' / &
         SNTags(793)%Label / 'gas_98_sa_min_un' / &
         SNTags(794)%Label / 'gas_98_sa_max' / &
         SNTags(795)%Label / 'gas_98_sa_source' / &
         SNTags(796)%Label / 'gas_99_sa_fmin' / &
         SNTags(797)%Label / 'gas_99_sa_fmax' / &
         SNTags(798)%Label / 'gas_99_sa_hfn_fmin' / &
         SNTags(799)%Label / 'gas_99_sa_min_st' / &
         SNTags(800)%Label / 'gas_99_sa_min_un' /
    data SNTags(801)%Label / 'gas_99_sa_max' / &
         SNTags(802)%Label / 'gas_99_sa_source' / &
         SNTags(803)%Label / 'gas_100_sa_fmin' / &
         SNTags(804)%Label / 'gas_100_sa_fmax' / &
         SNTags(805)%Label / 'gas_100_sa_hfn_fmin' / &
         SNTags(806)%Label / 'gas_100_sa_min_st' / &
         SNTags(807)%Label / 'gas_100_sa_min_un' / &
         SNTags(808)%Label / 'gas_100_sa_max' / &
         SNTags(809)%Label / 'gas_100_sa_source' / &
         SNTags(810)%Label / 'gas_101_sa_fmin' / &
         SNTags(811)%Label / 'gas_101_sa_fmax' / &
         SNTags(812)%Label / 'gas_101_sa_hfn_fmin' / &
         SNTags(813)%Label / 'gas_101_sa_min_st' / &
         SNTags(814)%Label / 'gas_101_sa_min_un' / &
         SNTags(815)%Label / 'gas_101_sa_max' / &
         SNTags(816)%Label / 'gas_101_sa_source' / &
         SNTags(817)%Label / 'gas_102_sa_fmin' / &
         SNTags(818)%Label / 'gas_102_sa_fmax' / &
         SNTags(819)%Label / 'gas_102_sa_hfn_fmin' / &
         SNTags(820)%Label / 'gas_102_sa_min_st' / &
         SNTags(821)%Label / 'gas_102_sa_min_un' / &
         SNTags(822)%Label / 'gas_102_sa_max' / &
         SNTags(823)%Label / 'gas_102_sa_source' / &
         SNTags(824)%Label / 'gas_103_sa_fmin' / &
         SNTags(825)%Label / 'gas_103_sa_fmax' / &
         SNTags(826)%Label / 'gas_103_sa_hfn_fmin' / &
         SNTags(827)%Label / 'gas_103_sa_min_st' / &
         SNTags(828)%Label / 'gas_103_sa_min_un' / &
         SNTags(829)%Label / 'gas_103_sa_max' / &
         SNTags(830)%Label / 'gas_103_sa_source' / &
         SNTags(831)%Label / 'gas_104_sa_fmin' / &
         SNTags(832)%Label / 'gas_104_sa_fmax' / &
         SNTags(833)%Label / 'gas_104_sa_hfn_fmin' / &
         SNTags(834)%Label / 'gas_104_sa_min_st' / &
         SNTags(835)%Label / 'gas_104_sa_min_un' / &
         SNTags(836)%Label / 'gas_104_sa_max' / &
         SNTags(837)%Label / 'gas_104_sa_source' / &
         SNTags(838)%Label / 'gas_105_sa_fmin' / &
         SNTags(839)%Label / 'gas_105_sa_fmax' / &
         SNTags(840)%Label / 'gas_105_sa_hfn_fmin' / &
         SNTags(841)%Label / 'gas_105_sa_min_st' / &
         SNTags(842)%Label / 'gas_105_sa_min_un' / &
         SNTags(843)%Label / 'gas_105_sa_max' / &
         SNTags(844)%Label / 'gas_105_sa_source' / &
         SNTags(845)%Label / 'gas_106_sa_fmin' / &
         SNTags(846)%Label / 'gas_106_sa_fmax' / &
         SNTags(847)%Label / 'gas_106_sa_hfn_fmin' / &
         SNTags(848)%Label / 'gas_106_sa_min_st' / &
         SNTags(849)%Label / 'gas_106_sa_min_un' / &
         SNTags(850)%Label / 'gas_106_sa_max' / &
         SNTags(851)%Label / 'gas_106_sa_source' / &
         SNTags(852)%Label / 'gas_107_sa_fmin' / &
         SNTags(853)%Label / 'gas_107_sa_fmax' / &
         SNTags(854)%Label / 'gas_107_sa_hfn_fmin' / &
         SNTags(855)%Label / 'gas_107_sa_min_st' / &
         SNTags(856)%Label / 'gas_107_sa_min_un' / &
         SNTags(857)%Label / 'gas_107_sa_max' / &
         SNTags(858)%Label / 'gas_107_sa_source' / &
         SNTags(859)%Label / 'gas_108_sa_fmin' / &
         SNTags(860)%Label / 'gas_108_sa_fmax' / &
         SNTags(861)%Label / 'gas_108_sa_hfn_fmin' / &
         SNTags(862)%Label / 'gas_108_sa_min_st' / &
         SNTags(863)%Label / 'gas_108_sa_min_un' / &
         SNTags(864)%Label / 'gas_108_sa_max' / &
         SNTags(865)%Label / 'gas_108_sa_source' / &
         SNTags(866)%Label / 'gas_109_sa_fmin' / &
         SNTags(867)%Label / 'gas_109_sa_fmax' / &
         SNTags(868)%Label / 'gas_109_sa_hfn_fmin' / &
         SNTags(869)%Label / 'gas_109_sa_min_st' / &
         SNTags(870)%Label / 'gas_109_sa_min_un' / &
         SNTags(871)%Label / 'gas_109_sa_max' / &
         SNTags(872)%Label / 'gas_109_sa_source' / &
         SNTags(873)%Label / 'gas_110_sa_fmin' / &
         SNTags(874)%Label / 'gas_110_sa_fmax' / &
         SNTags(875)%Label / 'gas_110_sa_hfn_fmin' / &
         SNTags(876)%Label / 'gas_110_sa_min_st' / &
         SNTags(877)%Label / 'gas_110_sa_min_un' / &
         SNTags(878)%Label / 'gas_110_sa_max' / &
         SNTags(879)%Label / 'gas_110_sa_source' / &
         SNTags(880)%Label / 'gas_111_sa_fmin' / &
         SNTags(881)%Label / 'gas_111_sa_fmax' / &
         SNTags(882)%Label / 'gas_111_sa_hfn_fmin' / &
         SNTags(883)%Label / 'gas_111_sa_min_st' / &
         SNTags(884)%Label / 'gas_111_sa_min_un' / &
         SNTags(885)%Label / 'gas_111_sa_max' / &
         SNTags(886)%Label / 'gas_111_sa_source' / &
         SNTags(887)%Label / 'gas_112_sa_fmin' / &
         SNTags(888)%Label / 'gas_112_sa_fmax' / &
         SNTags(889)%Label / 'gas_112_sa_hfn_fmin' / &
         SNTags(890)%Label / 'gas_112_sa_min_st' / &
         SNTags(891)%Label / 'gas_112_sa_min_un' / &
         SNTags(892)%Label / 'gas_112_sa_max' / &
         SNTags(893)%Label / 'gas_112_sa_source' / &
         SNTags(894)%Label / 'gas_113_sa_fmin' / &
         SNTags(895)%Label / 'gas_113_sa_fmax' / &
         SNTags(896)%Label / 'gas_113_sa_hfn_fmin' / &
         SNTags(897)%Label / 'gas_113_sa_min_st' / &
         SNTags(898)%Label / 'gas_113_sa_min_un' / &
         SNTags(899)%Label / 'gas_113_sa_max' / &
         SNTags(900)%Label / 'gas_113_sa_source' / &
         SNTags(901)%Label / 'gas_114_sa_fmin' / &
         SNTags(902)%Label / 'gas_114_sa_fmax' / &
         SNTags(903)%Label / 'gas_114_sa_hfn_fmin' / &
         SNTags(904)%Label / 'gas_114_sa_min_st' / &
         SNTags(905)%Label / 'gas_114_sa_min_un' / &
         SNTags(906)%Label / 'gas_114_sa_max' / &
         SNTags(907)%Label / 'gas_114_sa_source' / &
         SNTags(908)%Label / 'gas_115_sa_fmin' / &
         SNTags(909)%Label / 'gas_115_sa_fmax' / &
         SNTags(910)%Label / 'gas_115_sa_hfn_fmin' / &
         SNTags(911)%Label / 'gas_115_sa_min_st' / &
         SNTags(912)%Label / 'gas_115_sa_min_un' / &
         SNTags(913)%Label / 'gas_115_sa_max' / &
         SNTags(914)%Label / 'gas_115_sa_source' / &
         SNTags(915)%Label / 'gas_116_sa_fmin' / &
         SNTags(916)%Label / 'gas_116_sa_fmax' / &
         SNTags(917)%Label / 'gas_116_sa_hfn_fmin' / &
         SNTags(918)%Label / 'gas_116_sa_min_st' / &
         SNTags(919)%Label / 'gas_116_sa_min_un' / &
         SNTags(920)%Label / 'gas_116_sa_max' / &
         SNTags(921)%Label / 'gas_116_sa_source' / &
         SNTags(922)%Label / 'gas_117_sa_fmin' / &
         SNTags(923)%Label / 'gas_117_sa_fmax' / &
         SNTags(924)%Label / 'gas_117_sa_hfn_fmin' / &
         SNTags(925)%Label / 'gas_117_sa_min_st' / &
         SNTags(926)%Label / 'gas_117_sa_min_un' / &
         SNTags(927)%Label / 'gas_117_sa_max' / &
         SNTags(928)%Label / 'gas_117_sa_source' / &
         SNTags(929)%Label / 'gas_118_sa_fmin' / &
         SNTags(930)%Label / 'gas_118_sa_fmax' / &
         SNTags(931)%Label / 'gas_118_sa_hfn_fmin' / &
         SNTags(932)%Label / 'gas_118_sa_min_st' / &
         SNTags(933)%Label / 'gas_118_sa_min_un' / &
         SNTags(934)%Label / 'gas_118_sa_max' / &
         SNTags(935)%Label / 'gas_118_sa_source' / &
         SNTags(936)%Label / 'gas_119_sa_fmin' / &
         SNTags(937)%Label / 'gas_119_sa_fmax' / &
         SNTags(938)%Label / 'gas_119_sa_hfn_fmin' / &
         SNTags(939)%Label / 'gas_119_sa_min_st' / &
         SNTags(940)%Label / 'gas_119_sa_min_un' / &
         SNTags(941)%Label / 'gas_119_sa_max' / &
         SNTags(942)%Label / 'gas_119_sa_source' / &
         SNTags(943)%Label / 'gas_120_sa_fmin' / &
         SNTags(944)%Label / 'gas_120_sa_fmax' / &
         SNTags(945)%Label / 'gas_120_sa_hfn_fmin' / &
         SNTags(946)%Label / 'gas_120_sa_min_st' / &
         SNTags(947)%Label / 'gas_120_sa_min_un' / &
         SNTags(948)%Label / 'gas_120_sa_max' / &
         SNTags(949)%Label / 'gas_120_sa_source' / &
         SNTags(950)%Label / 'gas_121_sa_fmin' / &
         SNTags(951)%Label / 'gas_121_sa_fmax' / &
         SNTags(952)%Label / 'gas_121_sa_hfn_fmin' / &
         SNTags(953)%Label / 'gas_121_sa_min_st' / &
         SNTags(954)%Label / 'gas_121_sa_min_un' / &
         SNTags(955)%Label / 'gas_121_sa_max' / &
         SNTags(956)%Label / 'gas_121_sa_source' / &
         SNTags(957)%Label / 'gas_122_sa_fmin' / &
         SNTags(958)%Label / 'gas_122_sa_fmax' / &
         SNTags(959)%Label / 'gas_122_sa_hfn_fmin' / &
         SNTags(960)%Label / 'gas_122_sa_min_st' / &
         SNTags(961)%Label / 'gas_122_sa_min_un' / &
         SNTags(962)%Label / 'gas_122_sa_max' / &
         SNTags(963)%Label / 'gas_122_sa_source' / &
         SNTags(964)%Label / 'gas_123_sa_fmin' / &
         SNTags(965)%Label / 'gas_123_sa_fmax' / &
         SNTags(966)%Label / 'gas_123_sa_hfn_fmin' / &
         SNTags(967)%Label / 'gas_123_sa_min_st' / &
         SNTags(968)%Label / 'gas_123_sa_min_un' / &
         SNTags(969)%Label / 'gas_123_sa_max' / &
         SNTags(970)%Label / 'gas_123_sa_source' / &
         SNTags(971)%Label / 'gas_124_sa_fmin' / &
         SNTags(972)%Label / 'gas_124_sa_fmax' / &
         SNTags(973)%Label / 'gas_124_sa_hfn_fmin' / &
         SNTags(974)%Label / 'gas_124_sa_min_st' / &
         SNTags(975)%Label / 'gas_124_sa_min_un' / &
         SNTags(976)%Label / 'gas_124_sa_max' / &
         SNTags(977)%Label / 'gas_124_sa_source' / &
         SNTags(978)%Label / 'gas_125_sa_fmin' / &
         SNTags(979)%Label / 'gas_125_sa_fmax' / &
         SNTags(980)%Label / 'gas_125_sa_hfn_fmin' / &
         SNTags(981)%Label / 'gas_125_sa_min_st' / &
         SNTags(982)%Label / 'gas_125_sa_min_un' / &
         SNTags(983)%Label / 'gas_125_sa_max' / &
         SNTags(984)%Label / 'gas_125_sa_source' / &
         SNTags(985)%Label / 'gas_126_sa_fmin' / &
         SNTags(986)%Label / 'gas_126_sa_fmax' / &
         SNTags(987)%Label / 'gas_126_sa_hfn_fmin' / &
         SNTags(988)%Label / 'gas_126_sa_min_st' / &
         SNTags(989)%Label / 'gas_126_sa_min_un' / &
         SNTags(990)%Label / 'gas_126_sa_max' / &
         SNTags(991)%Label / 'gas_126_sa_source' / &
         SNTags(992)%Label / 'gas_127_sa_fmin' / &
         SNTags(993)%Label / 'gas_127_sa_fmax' / &
         SNTags(994)%Label / 'gas_127_sa_hfn_fmin' / &
         SNTags(995)%Label / 'gas_127_sa_min_st' / &
         SNTags(996)%Label / 'gas_127_sa_min_un' / &
         SNTags(997)%Label / 'gas_127_sa_max' / &
         SNTags(998)%Label / 'gas_127_sa_source' / &
         SNTags(999)%Label / 'gas_128_sa_fmin' / &
         SNTags(1000)%Label / 'gas_128_sa_fmax' /
    data SNTags(1001)%Label / 'gas_128_sa_hfn_fmin' / &
         SNTags(1002)%Label / 'gas_128_sa_min_st' / &
         SNTags(1003)%Label / 'gas_128_sa_min_un' / &
         SNTags(1004)%Label / 'gas_128_sa_max' / &
         SNTags(1005)%Label / 'gas_128_sa_source' /
    !> END GENERATED FCC.SNTags

    !> BEGIN GENERATED FCC.SCTags - edit gen_project_tags.py, not this block
    data SCTags(1)%Label   / 'sa_start_date' / &
         SCTags(2)%Label   / 'sa_start_time' / &
         SCTags(3)%Label   / 'sa_end_date' / &
         SCTags(4)%Label   / 'sa_end_time' / &
         SCTags(5)%Label   / 'start_sa_date' / &
         SCTags(6)%Label   / 'end_sa_date' / &
         SCTags(7)%Label   / 'ex_file' / &
         SCTags(8)%Label   / 'sa_bin_spectra' / &
         SCTags(9)%Label   / 'sa_full_spectra' / &
         SCTags(10)%Label  / 'out_path' / &
         SCTags(11)%Label  / 'make_dataset' / &
         SCTags(12)%Label  / 'sa_hfn_elim' / &
         SCTags(13)%Label  / 'h_meth' / &
         SCTags(14)%Label  / 'hf_meth' / &
         SCTags(15)%Label  / 'lptf_model' / &
         SCTags(16)%Label  / 'cosp_model' / &
         SCTags(17)%Label  / 'add_sonic_lptf' / &
         SCTags(18)%Label  / 'lf_meth' / &
         SCTags(19)%Label  / 'sa_mode' / &
         SCTags(20)%Label  / 'sa_file' / &
         SCTags(21)%Label  / 'horst_lens' / &
         SCTags(22)%Label  / 'sa_subset' / &
         SCTags(23)%Label  / 'sa_use_vm_flags' / &
         SCTags(24)%Label  / 'sa_use_foken_low' / &
         SCTags(25)%Label  / 'sa_use_foken_mid' / &
         SCTags(26)%Label  / 'keep_parent_fluxnet_file' / &
         SCTags(27)%Label  / 'automatic_spectra_config' / &
         SCTags(28)%Label  / 'gas_1_sa_months' / &
         SCTags(29)%Label  / 'gas_2_sa_months' / &
         SCTags(30)%Label  / 'gas_3_sa_months' / &
         SCTags(31)%Label  / 'gas_4_sa_months' / &
         SCTags(32)%Label  / 'gas_5_sa_months' / &
         SCTags(33)%Label  / 'gas_6_sa_months' / &
         SCTags(34)%Label  / 'gas_7_sa_months' / &
         SCTags(35)%Label  / 'gas_8_sa_months' / &
         SCTags(36)%Label  / 'gas_9_sa_months' / &
         SCTags(37)%Label  / 'gas_10_sa_months' / &
         SCTags(38)%Label  / 'gas_11_sa_months' / &
         SCTags(39)%Label  / 'gas_12_sa_months' / &
         SCTags(40)%Label  / 'gas_13_sa_months' / &
         SCTags(41)%Label  / 'gas_14_sa_months' / &
         SCTags(42)%Label  / 'gas_15_sa_months' / &
         SCTags(43)%Label  / 'gas_16_sa_months' / &
         SCTags(44)%Label  / 'gas_17_sa_months' / &
         SCTags(45)%Label  / 'gas_18_sa_months' / &
         SCTags(46)%Label  / 'gas_19_sa_months' / &
         SCTags(47)%Label  / 'gas_20_sa_months' / &
         SCTags(48)%Label  / 'gas_21_sa_months' / &
         SCTags(49)%Label  / 'gas_22_sa_months' / &
         SCTags(50)%Label  / 'gas_23_sa_months' / &
         SCTags(51)%Label  / 'gas_24_sa_months' / &
         SCTags(52)%Label  / 'gas_25_sa_months' / &
         SCTags(53)%Label  / 'gas_26_sa_months' / &
         SCTags(54)%Label  / 'gas_27_sa_months' / &
         SCTags(55)%Label  / 'gas_28_sa_months' / &
         SCTags(56)%Label  / 'gas_29_sa_months' / &
         SCTags(57)%Label  / 'gas_30_sa_months' / &
         SCTags(58)%Label  / 'gas_31_sa_months' / &
         SCTags(59)%Label  / 'gas_32_sa_months' / &
         SCTags(60)%Label  / 'gas_33_sa_months' / &
         SCTags(61)%Label  / 'gas_34_sa_months' / &
         SCTags(62)%Label  / 'gas_35_sa_months' / &
         SCTags(63)%Label  / 'gas_36_sa_months' / &
         SCTags(64)%Label  / 'gas_37_sa_months' / &
         SCTags(65)%Label  / 'gas_38_sa_months' / &
         SCTags(66)%Label  / 'gas_39_sa_months' / &
         SCTags(67)%Label  / 'gas_40_sa_months' / &
         SCTags(68)%Label  / 'gas_41_sa_months' / &
         SCTags(69)%Label  / 'gas_42_sa_months' / &
         SCTags(70)%Label  / 'gas_43_sa_months' / &
         SCTags(71)%Label  / 'gas_44_sa_months' / &
         SCTags(72)%Label  / 'gas_45_sa_months' / &
         SCTags(73)%Label  / 'gas_46_sa_months' / &
         SCTags(74)%Label  / 'gas_47_sa_months' / &
         SCTags(75)%Label  / 'gas_48_sa_months' / &
         SCTags(76)%Label  / 'gas_49_sa_months' / &
         SCTags(77)%Label  / 'gas_50_sa_months' / &
         SCTags(78)%Label  / 'gas_51_sa_months' / &
         SCTags(79)%Label  / 'gas_52_sa_months' / &
         SCTags(80)%Label  / 'gas_53_sa_months' / &
         SCTags(81)%Label  / 'gas_54_sa_months' / &
         SCTags(82)%Label  / 'gas_55_sa_months' / &
         SCTags(83)%Label  / 'gas_56_sa_months' / &
         SCTags(84)%Label  / 'gas_57_sa_months' / &
         SCTags(85)%Label  / 'gas_58_sa_months' / &
         SCTags(86)%Label  / 'gas_59_sa_months' / &
         SCTags(87)%Label  / 'gas_60_sa_months' / &
         SCTags(88)%Label  / 'gas_61_sa_months' / &
         SCTags(89)%Label  / 'gas_62_sa_months' / &
         SCTags(90)%Label  / 'gas_63_sa_months' / &
         SCTags(91)%Label  / 'gas_64_sa_months' / &
         SCTags(92)%Label  / 'gas_65_sa_months' / &
         SCTags(93)%Label  / 'gas_66_sa_months' / &
         SCTags(94)%Label  / 'gas_67_sa_months' / &
         SCTags(95)%Label  / 'gas_68_sa_months' / &
         SCTags(96)%Label  / 'gas_69_sa_months' / &
         SCTags(97)%Label  / 'gas_70_sa_months' / &
         SCTags(98)%Label  / 'gas_71_sa_months' / &
         SCTags(99)%Label  / 'gas_72_sa_months' / &
         SCTags(100)%Label / 'gas_73_sa_months' / &
         SCTags(101)%Label / 'gas_74_sa_months' / &
         SCTags(102)%Label / 'gas_75_sa_months' / &
         SCTags(103)%Label / 'gas_76_sa_months' / &
         SCTags(104)%Label / 'gas_77_sa_months' / &
         SCTags(105)%Label / 'gas_78_sa_months' / &
         SCTags(106)%Label / 'gas_79_sa_months' / &
         SCTags(107)%Label / 'gas_80_sa_months' / &
         SCTags(108)%Label / 'gas_81_sa_months' / &
         SCTags(109)%Label / 'gas_82_sa_months' / &
         SCTags(110)%Label / 'gas_83_sa_months' / &
         SCTags(111)%Label / 'gas_84_sa_months' / &
         SCTags(112)%Label / 'gas_85_sa_months' / &
         SCTags(113)%Label / 'gas_86_sa_months' / &
         SCTags(114)%Label / 'gas_87_sa_months' / &
         SCTags(115)%Label / 'gas_88_sa_months' / &
         SCTags(116)%Label / 'gas_89_sa_months' / &
         SCTags(117)%Label / 'gas_90_sa_months' / &
         SCTags(118)%Label / 'gas_91_sa_months' / &
         SCTags(119)%Label / 'gas_92_sa_months' / &
         SCTags(120)%Label / 'gas_93_sa_months' / &
         SCTags(121)%Label / 'gas_94_sa_months' / &
         SCTags(122)%Label / 'gas_95_sa_months' / &
         SCTags(123)%Label / 'gas_96_sa_months' / &
         SCTags(124)%Label / 'gas_97_sa_months' / &
         SCTags(125)%Label / 'gas_98_sa_months' / &
         SCTags(126)%Label / 'gas_99_sa_months' / &
         SCTags(127)%Label / 'gas_100_sa_months' / &
         SCTags(128)%Label / 'gas_101_sa_months' / &
         SCTags(129)%Label / 'gas_102_sa_months' / &
         SCTags(130)%Label / 'gas_103_sa_months' / &
         SCTags(131)%Label / 'gas_104_sa_months' / &
         SCTags(132)%Label / 'gas_105_sa_months' / &
         SCTags(133)%Label / 'gas_106_sa_months' / &
         SCTags(134)%Label / 'gas_107_sa_months' / &
         SCTags(135)%Label / 'gas_108_sa_months' / &
         SCTags(136)%Label / 'gas_109_sa_months' / &
         SCTags(137)%Label / 'gas_110_sa_months' / &
         SCTags(138)%Label / 'gas_111_sa_months' / &
         SCTags(139)%Label / 'gas_112_sa_months' / &
         SCTags(140)%Label / 'gas_113_sa_months' / &
         SCTags(141)%Label / 'gas_114_sa_months' / &
         SCTags(142)%Label / 'gas_115_sa_months' / &
         SCTags(143)%Label / 'gas_116_sa_months' / &
         SCTags(144)%Label / 'gas_117_sa_months' / &
         SCTags(145)%Label / 'gas_118_sa_months' / &
         SCTags(146)%Label / 'gas_119_sa_months' / &
         SCTags(147)%Label / 'gas_120_sa_months' / &
         SCTags(148)%Label / 'gas_121_sa_months' / &
         SCTags(149)%Label / 'gas_122_sa_months' / &
         SCTags(150)%Label / 'gas_123_sa_months' / &
         SCTags(151)%Label / 'gas_124_sa_months' / &
         SCTags(152)%Label / 'gas_125_sa_months' / &
         SCTags(153)%Label / 'gas_126_sa_months' / &
         SCTags(154)%Label / 'gas_127_sa_months' / &
         SCTags(155)%Label / 'gas_128_sa_months' /
    !> END GENERATED FCC.SCTags
end module m_fx_global_var
