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
    integer :: g4l
    integer :: nstat
    integer :: nfull
    integer, parameter :: ndkf = 60     !< TO BE INCREASED????!!!!!!!
    integer, parameter :: unstable = 1
    integer, parameter :: stable   = 2

    real(kind = dbl) :: float_doy
    real(kind = dbl) :: dkf(ndkf + 1)
    real(kind = dbl) :: gas4_full_flux_sc
    real(kind = dbl) :: gas4_full_dens_sc
    real(kind = dbl), allocatable :: custVars(:)

    character(12), parameter :: fcc_app = 'EddyFlow-FCC'
    character(32) :: g4lab
    character(32) :: gas4_full_flux_label
    character(32) :: gas4_full_conc_label
    character(32) :: gas4_full_mixr_label
    character(32) :: gas4_full_dens_label
    character(64) :: UserVarHeader(MaxUserVar)
    character(26), parameter :: SubDirSpecAn = 'eddyflow_spectral_analysis'
    character(16000) :: fluxnet_header

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
    type(CECFluxType) :: CECFlux
    type(FCCMetadataType) :: FCCMetadata

    !> tags of the setup ".ini" file for eccoce
    integer, parameter :: Nsn = 493
    integer, parameter :: Nsc = 31
    logical            :: SNTagFound(Nsn)
    logical            :: SCTagFound(Nsc)
    type (Numerical)   :: SNTags(Nsn)
    type (Text)        :: SCTags(Nsc)
    !> BEGIN GENERATED FCC.SNTags - edit gen_project_tags.py, not this block
    data SNTags(1)%Label   / 'sa_nbins' / &
         SNTags(2)%Label   / 'sa_min_smpl' / &
         SNTags(3)%Label   / 'sa_fmin_co2' / &
         SNTags(4)%Label   / 'sa_fmax_co2' / &
         SNTags(5)%Label   / 'sa_fmin_h2o' / &
         SNTags(6)%Label   / 'sa_fmax_h2o' / &
         SNTags(7)%Label   / 'sa_fmin_ch4' / &
         SNTags(8)%Label   / 'sa_fmax_ch4' / &
         SNTags(9)%Label   / 'sa_fmin_gas4' / &
         SNTags(10)%Label  / 'sa_fmax_gas4' / &
         SNTags(11)%Label  / 'sa_min_co2' / &
         SNTags(12)%Label  / 'sa_min_ch4' / &
         SNTags(13)%Label  / 'sa_min_gas4' / &
         SNTags(14)%Label  / 'sa_min_le' / &
         SNTags(15)%Label  / 'sa_min_h' / &
         SNTags(16)%Label  / 'sa_hfn_co2_fmin' / &
         SNTags(17)%Label  / 'sa_hfn_h2o_fmin' / &
         SNTags(18)%Label  / 'sa_hfn_ch4_fmin' / &
         SNTags(19)%Label  / 'sa_hfn_gas4_fmin' / &
         SNTags(20)%Label  / 'sa_co2_g1_start' / &
         SNTags(21)%Label  / 'sa_co2_g1_stop' / &
         SNTags(22)%Label  / 'sa_co2_g2_start' / &
         SNTags(23)%Label  / 'sa_co2_g2_stop' / &
         SNTags(24)%Label  / 'sa_co2_g3_start' / &
         SNTags(25)%Label  / 'sa_co2_g3_stop' / &
         SNTags(26)%Label  / 'sa_co2_g4_start' / &
         SNTags(27)%Label  / 'sa_co2_g4_stop' / &
         SNTags(28)%Label  / 'sa_co2_g5_start' / &
         SNTags(29)%Label  / 'sa_co2_g5_stop' / &
         SNTags(30)%Label  / 'sa_co2_g6_start' / &
         SNTags(31)%Label  / 'sa_co2_g6_stop' / &
         SNTags(32)%Label  / 'sa_co2_g7_start' / &
         SNTags(33)%Label  / 'sa_co2_g7_stop' / &
         SNTags(34)%Label  / 'sa_co2_g8_start' / &
         SNTags(35)%Label  / 'sa_co2_g8_stop' / &
         SNTags(36)%Label  / 'sa_co2_g9_start' / &
         SNTags(37)%Label  / 'sa_co2_g9_stop' / &
         SNTags(38)%Label  / 'sa_co2_g10_start' / &
         SNTags(39)%Label  / 'sa_co2_g10_stop' / &
         SNTags(40)%Label  / 'sa_co2_g11_start' / &
         SNTags(41)%Label  / 'sa_co2_g11_stop' / &
         SNTags(42)%Label  / 'sa_co2_g12_start' / &
         SNTags(43)%Label  / 'sa_co2_g12_stop' / &
         SNTags(44)%Label  / 'sa_ch4_g1_start' / &
         SNTags(45)%Label  / 'sa_ch4_g1_stop' / &
         SNTags(46)%Label  / 'sa_ch4_g2_start' / &
         SNTags(47)%Label  / 'sa_ch4_g2_stop' / &
         SNTags(48)%Label  / 'sa_ch4_g3_start' / &
         SNTags(49)%Label  / 'sa_ch4_g3_stop' / &
         SNTags(50)%Label  / 'sa_ch4_g4_start' / &
         SNTags(51)%Label  / 'sa_ch4_g4_stop' / &
         SNTags(52)%Label  / 'sa_ch4_g5_start' / &
         SNTags(53)%Label  / 'sa_ch4_g5_stop' / &
         SNTags(54)%Label  / 'sa_ch4_g6_start' / &
         SNTags(55)%Label  / 'sa_ch4_g6_stop' / &
         SNTags(56)%Label  / 'sa_ch4_g7_start' / &
         SNTags(57)%Label  / 'sa_ch4_g7_stop' / &
         SNTags(58)%Label  / 'sa_ch4_g8_start' / &
         SNTags(59)%Label  / 'sa_ch4_g8_stop' / &
         SNTags(60)%Label  / 'sa_ch4_g9_start' / &
         SNTags(61)%Label  / 'sa_ch4_g9_stop' / &
         SNTags(62)%Label  / 'sa_ch4_g10_start' / &
         SNTags(63)%Label  / 'sa_ch4_g10_stop' / &
         SNTags(64)%Label  / 'sa_ch4_g11_start' / &
         SNTags(65)%Label  / 'sa_ch4_g11_stop' / &
         SNTags(66)%Label  / 'sa_ch4_g12_start' / &
         SNTags(67)%Label  / 'sa_ch4_g12_stop' / &
         SNTags(68)%Label  / 'sa_gas4_g1_start' / &
         SNTags(69)%Label  / 'sa_gas4_g1_stop' / &
         SNTags(70)%Label  / 'sa_gas4_g2_start' / &
         SNTags(71)%Label  / 'sa_gas4_g2_stop' / &
         SNTags(72)%Label  / 'sa_gas4_g3_start' / &
         SNTags(73)%Label  / 'sa_gas4_g3_stop' / &
         SNTags(74)%Label  / 'sa_gas4_g4_start' / &
         SNTags(75)%Label  / 'sa_gas4_g4_stop' / &
         SNTags(76)%Label  / 'sa_gas4_g5_start' / &
         SNTags(77)%Label  / 'sa_gas4_g5_stop' / &
         SNTags(78)%Label  / 'sa_gas4_g6_start' / &
         SNTags(79)%Label  / 'sa_gas4_g6_stop' / &
         SNTags(80)%Label  / 'sa_gas4_g7_start' / &
         SNTags(81)%Label  / 'sa_gas4_g7_stop' / &
         SNTags(82)%Label  / 'sa_gas4_g8_start' / &
         SNTags(83)%Label  / 'sa_gas4_g8_stop' / &
         SNTags(84)%Label  / 'sa_gas4_g9_start' / &
         SNTags(85)%Label  / 'sa_gas4_g9_stop' / &
         SNTags(86)%Label  / 'sa_gas4_g10_start' / &
         SNTags(87)%Label  / 'sa_gas4_g10_stop' / &
         SNTags(88)%Label  / 'sa_gas4_g11_start' / &
         SNTags(89)%Label  / 'sa_gas4_g11_stop' / &
         SNTags(90)%Label  / 'sa_gas4_g12_start' / &
         SNTags(91)%Label  / 'sa_gas4_g12_stop' / &
         SNTags(92)%Label  / 'sa_min_un_ustar' / &
         SNTags(93)%Label  / 'sa_min_un_co2' / &
         SNTags(94)%Label  / 'sa_min_un_ch4' / &
         SNTags(95)%Label  / 'sa_min_un_gas4' / &
         SNTags(96)%Label  / 'sa_min_un_le' / &
         SNTags(97)%Label  / 'sa_min_un_h' / &
         SNTags(98)%Label  / 'sa_min_st_ustar' / &
         SNTags(99)%Label  / 'sa_min_st_co2' / &
         SNTags(100)%Label / 'sa_min_st_ch4' / &
         SNTags(101)%Label / 'sa_min_st_gas4' / &
         SNTags(102)%Label / 'sa_min_st_le' / &
         SNTags(103)%Label / 'sa_min_st_h' / &
         SNTags(104)%Label / 'sa_max_ustar' / &
         SNTags(105)%Label / 'sa_max_co2' / &
         SNTags(106)%Label / 'sa_max_ch4' / &
         SNTags(107)%Label / 'sa_max_gas4' / &
         SNTags(108)%Label / 'sa_max_le' / &
         SNTags(109)%Label / 'sa_max_h' / &
         SNTags(110)%Label / 'gas_1_sa_fmin' / &
         SNTags(111)%Label / 'gas_1_sa_fmax' / &
         SNTags(112)%Label / 'gas_1_sa_hfn_fmin' / &
         SNTags(113)%Label / 'gas_1_sa_min_st' / &
         SNTags(114)%Label / 'gas_1_sa_min_un' / &
         SNTags(115)%Label / 'gas_1_sa_max' / &
         SNTags(116)%Label / 'gas_2_sa_fmin' / &
         SNTags(117)%Label / 'gas_2_sa_fmax' / &
         SNTags(118)%Label / 'gas_2_sa_hfn_fmin' / &
         SNTags(119)%Label / 'gas_2_sa_min_st' / &
         SNTags(120)%Label / 'gas_2_sa_min_un' / &
         SNTags(121)%Label / 'gas_2_sa_max' / &
         SNTags(122)%Label / 'gas_3_sa_fmin' / &
         SNTags(123)%Label / 'gas_3_sa_fmax' / &
         SNTags(124)%Label / 'gas_3_sa_hfn_fmin' / &
         SNTags(125)%Label / 'gas_3_sa_min_st' / &
         SNTags(126)%Label / 'gas_3_sa_min_un' / &
         SNTags(127)%Label / 'gas_3_sa_max' / &
         SNTags(128)%Label / 'gas_4_sa_fmin' / &
         SNTags(129)%Label / 'gas_4_sa_fmax' / &
         SNTags(130)%Label / 'gas_4_sa_hfn_fmin' / &
         SNTags(131)%Label / 'gas_4_sa_min_st' / &
         SNTags(132)%Label / 'gas_4_sa_min_un' / &
         SNTags(133)%Label / 'gas_4_sa_max' / &
         SNTags(134)%Label / 'gas_5_sa_fmin' / &
         SNTags(135)%Label / 'gas_5_sa_fmax' / &
         SNTags(136)%Label / 'gas_5_sa_hfn_fmin' / &
         SNTags(137)%Label / 'gas_5_sa_min_st' / &
         SNTags(138)%Label / 'gas_5_sa_min_un' / &
         SNTags(139)%Label / 'gas_5_sa_max' / &
         SNTags(140)%Label / 'gas_6_sa_fmin' / &
         SNTags(141)%Label / 'gas_6_sa_fmax' / &
         SNTags(142)%Label / 'gas_6_sa_hfn_fmin' / &
         SNTags(143)%Label / 'gas_6_sa_min_st' / &
         SNTags(144)%Label / 'gas_6_sa_min_un' / &
         SNTags(145)%Label / 'gas_6_sa_max' / &
         SNTags(146)%Label / 'gas_7_sa_fmin' / &
         SNTags(147)%Label / 'gas_7_sa_fmax' / &
         SNTags(148)%Label / 'gas_7_sa_hfn_fmin' / &
         SNTags(149)%Label / 'gas_7_sa_min_st' / &
         SNTags(150)%Label / 'gas_7_sa_min_un' / &
         SNTags(151)%Label / 'gas_7_sa_max' / &
         SNTags(152)%Label / 'gas_8_sa_fmin' / &
         SNTags(153)%Label / 'gas_8_sa_fmax' / &
         SNTags(154)%Label / 'gas_8_sa_hfn_fmin' / &
         SNTags(155)%Label / 'gas_8_sa_min_st' / &
         SNTags(156)%Label / 'gas_8_sa_min_un' / &
         SNTags(157)%Label / 'gas_8_sa_max' / &
         SNTags(158)%Label / 'gas_9_sa_fmin' / &
         SNTags(159)%Label / 'gas_9_sa_fmax' / &
         SNTags(160)%Label / 'gas_9_sa_hfn_fmin' / &
         SNTags(161)%Label / 'gas_9_sa_min_st' / &
         SNTags(162)%Label / 'gas_9_sa_min_un' / &
         SNTags(163)%Label / 'gas_9_sa_max' / &
         SNTags(164)%Label / 'gas_10_sa_fmin' / &
         SNTags(165)%Label / 'gas_10_sa_fmax' / &
         SNTags(166)%Label / 'gas_10_sa_hfn_fmin' / &
         SNTags(167)%Label / 'gas_10_sa_min_st' / &
         SNTags(168)%Label / 'gas_10_sa_min_un' / &
         SNTags(169)%Label / 'gas_10_sa_max' / &
         SNTags(170)%Label / 'gas_11_sa_fmin' / &
         SNTags(171)%Label / 'gas_11_sa_fmax' / &
         SNTags(172)%Label / 'gas_11_sa_hfn_fmin' / &
         SNTags(173)%Label / 'gas_11_sa_min_st' / &
         SNTags(174)%Label / 'gas_11_sa_min_un' / &
         SNTags(175)%Label / 'gas_11_sa_max' / &
         SNTags(176)%Label / 'gas_12_sa_fmin' / &
         SNTags(177)%Label / 'gas_12_sa_fmax' / &
         SNTags(178)%Label / 'gas_12_sa_hfn_fmin' / &
         SNTags(179)%Label / 'gas_12_sa_min_st' / &
         SNTags(180)%Label / 'gas_12_sa_min_un' / &
         SNTags(181)%Label / 'gas_12_sa_max' / &
         SNTags(182)%Label / 'gas_13_sa_fmin' / &
         SNTags(183)%Label / 'gas_13_sa_fmax' / &
         SNTags(184)%Label / 'gas_13_sa_hfn_fmin' / &
         SNTags(185)%Label / 'gas_13_sa_min_st' / &
         SNTags(186)%Label / 'gas_13_sa_min_un' / &
         SNTags(187)%Label / 'gas_13_sa_max' / &
         SNTags(188)%Label / 'gas_14_sa_fmin' / &
         SNTags(189)%Label / 'gas_14_sa_fmax' / &
         SNTags(190)%Label / 'gas_14_sa_hfn_fmin' / &
         SNTags(191)%Label / 'gas_14_sa_min_st' / &
         SNTags(192)%Label / 'gas_14_sa_min_un' / &
         SNTags(193)%Label / 'gas_14_sa_max' / &
         SNTags(194)%Label / 'gas_15_sa_fmin' / &
         SNTags(195)%Label / 'gas_15_sa_fmax' / &
         SNTags(196)%Label / 'gas_15_sa_hfn_fmin' / &
         SNTags(197)%Label / 'gas_15_sa_min_st' / &
         SNTags(198)%Label / 'gas_15_sa_min_un' / &
         SNTags(199)%Label / 'gas_15_sa_max' / &
         SNTags(200)%Label / 'gas_16_sa_fmin' /
    data SNTags(201)%Label / 'gas_16_sa_fmax' / &
         SNTags(202)%Label / 'gas_16_sa_hfn_fmin' / &
         SNTags(203)%Label / 'gas_16_sa_min_st' / &
         SNTags(204)%Label / 'gas_16_sa_min_un' / &
         SNTags(205)%Label / 'gas_16_sa_max' / &
         SNTags(206)%Label / 'gas_17_sa_fmin' / &
         SNTags(207)%Label / 'gas_17_sa_fmax' / &
         SNTags(208)%Label / 'gas_17_sa_hfn_fmin' / &
         SNTags(209)%Label / 'gas_17_sa_min_st' / &
         SNTags(210)%Label / 'gas_17_sa_min_un' / &
         SNTags(211)%Label / 'gas_17_sa_max' / &
         SNTags(212)%Label / 'gas_18_sa_fmin' / &
         SNTags(213)%Label / 'gas_18_sa_fmax' / &
         SNTags(214)%Label / 'gas_18_sa_hfn_fmin' / &
         SNTags(215)%Label / 'gas_18_sa_min_st' / &
         SNTags(216)%Label / 'gas_18_sa_min_un' / &
         SNTags(217)%Label / 'gas_18_sa_max' / &
         SNTags(218)%Label / 'gas_19_sa_fmin' / &
         SNTags(219)%Label / 'gas_19_sa_fmax' / &
         SNTags(220)%Label / 'gas_19_sa_hfn_fmin' / &
         SNTags(221)%Label / 'gas_19_sa_min_st' / &
         SNTags(222)%Label / 'gas_19_sa_min_un' / &
         SNTags(223)%Label / 'gas_19_sa_max' / &
         SNTags(224)%Label / 'gas_20_sa_fmin' / &
         SNTags(225)%Label / 'gas_20_sa_fmax' / &
         SNTags(226)%Label / 'gas_20_sa_hfn_fmin' / &
         SNTags(227)%Label / 'gas_20_sa_min_st' / &
         SNTags(228)%Label / 'gas_20_sa_min_un' / &
         SNTags(229)%Label / 'gas_20_sa_max' / &
         SNTags(230)%Label / 'gas_21_sa_fmin' / &
         SNTags(231)%Label / 'gas_21_sa_fmax' / &
         SNTags(232)%Label / 'gas_21_sa_hfn_fmin' / &
         SNTags(233)%Label / 'gas_21_sa_min_st' / &
         SNTags(234)%Label / 'gas_21_sa_min_un' / &
         SNTags(235)%Label / 'gas_21_sa_max' / &
         SNTags(236)%Label / 'gas_22_sa_fmin' / &
         SNTags(237)%Label / 'gas_22_sa_fmax' / &
         SNTags(238)%Label / 'gas_22_sa_hfn_fmin' / &
         SNTags(239)%Label / 'gas_22_sa_min_st' / &
         SNTags(240)%Label / 'gas_22_sa_min_un' / &
         SNTags(241)%Label / 'gas_22_sa_max' / &
         SNTags(242)%Label / 'gas_23_sa_fmin' / &
         SNTags(243)%Label / 'gas_23_sa_fmax' / &
         SNTags(244)%Label / 'gas_23_sa_hfn_fmin' / &
         SNTags(245)%Label / 'gas_23_sa_min_st' / &
         SNTags(246)%Label / 'gas_23_sa_min_un' / &
         SNTags(247)%Label / 'gas_23_sa_max' / &
         SNTags(248)%Label / 'gas_24_sa_fmin' / &
         SNTags(249)%Label / 'gas_24_sa_fmax' / &
         SNTags(250)%Label / 'gas_24_sa_hfn_fmin' / &
         SNTags(251)%Label / 'gas_24_sa_min_st' / &
         SNTags(252)%Label / 'gas_24_sa_min_un' / &
         SNTags(253)%Label / 'gas_24_sa_max' / &
         SNTags(254)%Label / 'gas_25_sa_fmin' / &
         SNTags(255)%Label / 'gas_25_sa_fmax' / &
         SNTags(256)%Label / 'gas_25_sa_hfn_fmin' / &
         SNTags(257)%Label / 'gas_25_sa_min_st' / &
         SNTags(258)%Label / 'gas_25_sa_min_un' / &
         SNTags(259)%Label / 'gas_25_sa_max' / &
         SNTags(260)%Label / 'gas_26_sa_fmin' / &
         SNTags(261)%Label / 'gas_26_sa_fmax' / &
         SNTags(262)%Label / 'gas_26_sa_hfn_fmin' / &
         SNTags(263)%Label / 'gas_26_sa_min_st' / &
         SNTags(264)%Label / 'gas_26_sa_min_un' / &
         SNTags(265)%Label / 'gas_26_sa_max' / &
         SNTags(266)%Label / 'gas_27_sa_fmin' / &
         SNTags(267)%Label / 'gas_27_sa_fmax' / &
         SNTags(268)%Label / 'gas_27_sa_hfn_fmin' / &
         SNTags(269)%Label / 'gas_27_sa_min_st' / &
         SNTags(270)%Label / 'gas_27_sa_min_un' / &
         SNTags(271)%Label / 'gas_27_sa_max' / &
         SNTags(272)%Label / 'gas_28_sa_fmin' / &
         SNTags(273)%Label / 'gas_28_sa_fmax' / &
         SNTags(274)%Label / 'gas_28_sa_hfn_fmin' / &
         SNTags(275)%Label / 'gas_28_sa_min_st' / &
         SNTags(276)%Label / 'gas_28_sa_min_un' / &
         SNTags(277)%Label / 'gas_28_sa_max' / &
         SNTags(278)%Label / 'gas_29_sa_fmin' / &
         SNTags(279)%Label / 'gas_29_sa_fmax' / &
         SNTags(280)%Label / 'gas_29_sa_hfn_fmin' / &
         SNTags(281)%Label / 'gas_29_sa_min_st' / &
         SNTags(282)%Label / 'gas_29_sa_min_un' / &
         SNTags(283)%Label / 'gas_29_sa_max' / &
         SNTags(284)%Label / 'gas_30_sa_fmin' / &
         SNTags(285)%Label / 'gas_30_sa_fmax' / &
         SNTags(286)%Label / 'gas_30_sa_hfn_fmin' / &
         SNTags(287)%Label / 'gas_30_sa_min_st' / &
         SNTags(288)%Label / 'gas_30_sa_min_un' / &
         SNTags(289)%Label / 'gas_30_sa_max' / &
         SNTags(290)%Label / 'gas_31_sa_fmin' / &
         SNTags(291)%Label / 'gas_31_sa_fmax' / &
         SNTags(292)%Label / 'gas_31_sa_hfn_fmin' / &
         SNTags(293)%Label / 'gas_31_sa_min_st' / &
         SNTags(294)%Label / 'gas_31_sa_min_un' / &
         SNTags(295)%Label / 'gas_31_sa_max' / &
         SNTags(296)%Label / 'gas_32_sa_fmin' / &
         SNTags(297)%Label / 'gas_32_sa_fmax' / &
         SNTags(298)%Label / 'gas_32_sa_hfn_fmin' / &
         SNTags(299)%Label / 'gas_32_sa_min_st' / &
         SNTags(300)%Label / 'gas_32_sa_min_un' / &
         SNTags(301)%Label / 'gas_32_sa_max' / &
         SNTags(302)%Label / 'gas_33_sa_fmin' / &
         SNTags(303)%Label / 'gas_33_sa_fmax' / &
         SNTags(304)%Label / 'gas_33_sa_hfn_fmin' / &
         SNTags(305)%Label / 'gas_33_sa_min_st' / &
         SNTags(306)%Label / 'gas_33_sa_min_un' / &
         SNTags(307)%Label / 'gas_33_sa_max' / &
         SNTags(308)%Label / 'gas_34_sa_fmin' / &
         SNTags(309)%Label / 'gas_34_sa_fmax' / &
         SNTags(310)%Label / 'gas_34_sa_hfn_fmin' / &
         SNTags(311)%Label / 'gas_34_sa_min_st' / &
         SNTags(312)%Label / 'gas_34_sa_min_un' / &
         SNTags(313)%Label / 'gas_34_sa_max' / &
         SNTags(314)%Label / 'gas_35_sa_fmin' / &
         SNTags(315)%Label / 'gas_35_sa_fmax' / &
         SNTags(316)%Label / 'gas_35_sa_hfn_fmin' / &
         SNTags(317)%Label / 'gas_35_sa_min_st' / &
         SNTags(318)%Label / 'gas_35_sa_min_un' / &
         SNTags(319)%Label / 'gas_35_sa_max' / &
         SNTags(320)%Label / 'gas_36_sa_fmin' / &
         SNTags(321)%Label / 'gas_36_sa_fmax' / &
         SNTags(322)%Label / 'gas_36_sa_hfn_fmin' / &
         SNTags(323)%Label / 'gas_36_sa_min_st' / &
         SNTags(324)%Label / 'gas_36_sa_min_un' / &
         SNTags(325)%Label / 'gas_36_sa_max' / &
         SNTags(326)%Label / 'gas_37_sa_fmin' / &
         SNTags(327)%Label / 'gas_37_sa_fmax' / &
         SNTags(328)%Label / 'gas_37_sa_hfn_fmin' / &
         SNTags(329)%Label / 'gas_37_sa_min_st' / &
         SNTags(330)%Label / 'gas_37_sa_min_un' / &
         SNTags(331)%Label / 'gas_37_sa_max' / &
         SNTags(332)%Label / 'gas_38_sa_fmin' / &
         SNTags(333)%Label / 'gas_38_sa_fmax' / &
         SNTags(334)%Label / 'gas_38_sa_hfn_fmin' / &
         SNTags(335)%Label / 'gas_38_sa_min_st' / &
         SNTags(336)%Label / 'gas_38_sa_min_un' / &
         SNTags(337)%Label / 'gas_38_sa_max' / &
         SNTags(338)%Label / 'gas_39_sa_fmin' / &
         SNTags(339)%Label / 'gas_39_sa_fmax' / &
         SNTags(340)%Label / 'gas_39_sa_hfn_fmin' / &
         SNTags(341)%Label / 'gas_39_sa_min_st' / &
         SNTags(342)%Label / 'gas_39_sa_min_un' / &
         SNTags(343)%Label / 'gas_39_sa_max' / &
         SNTags(344)%Label / 'gas_40_sa_fmin' / &
         SNTags(345)%Label / 'gas_40_sa_fmax' / &
         SNTags(346)%Label / 'gas_40_sa_hfn_fmin' / &
         SNTags(347)%Label / 'gas_40_sa_min_st' / &
         SNTags(348)%Label / 'gas_40_sa_min_un' / &
         SNTags(349)%Label / 'gas_40_sa_max' / &
         SNTags(350)%Label / 'gas_41_sa_fmin' / &
         SNTags(351)%Label / 'gas_41_sa_fmax' / &
         SNTags(352)%Label / 'gas_41_sa_hfn_fmin' / &
         SNTags(353)%Label / 'gas_41_sa_min_st' / &
         SNTags(354)%Label / 'gas_41_sa_min_un' / &
         SNTags(355)%Label / 'gas_41_sa_max' / &
         SNTags(356)%Label / 'gas_42_sa_fmin' / &
         SNTags(357)%Label / 'gas_42_sa_fmax' / &
         SNTags(358)%Label / 'gas_42_sa_hfn_fmin' / &
         SNTags(359)%Label / 'gas_42_sa_min_st' / &
         SNTags(360)%Label / 'gas_42_sa_min_un' / &
         SNTags(361)%Label / 'gas_42_sa_max' / &
         SNTags(362)%Label / 'gas_43_sa_fmin' / &
         SNTags(363)%Label / 'gas_43_sa_fmax' / &
         SNTags(364)%Label / 'gas_43_sa_hfn_fmin' / &
         SNTags(365)%Label / 'gas_43_sa_min_st' / &
         SNTags(366)%Label / 'gas_43_sa_min_un' / &
         SNTags(367)%Label / 'gas_43_sa_max' / &
         SNTags(368)%Label / 'gas_44_sa_fmin' / &
         SNTags(369)%Label / 'gas_44_sa_fmax' / &
         SNTags(370)%Label / 'gas_44_sa_hfn_fmin' / &
         SNTags(371)%Label / 'gas_44_sa_min_st' / &
         SNTags(372)%Label / 'gas_44_sa_min_un' / &
         SNTags(373)%Label / 'gas_44_sa_max' / &
         SNTags(374)%Label / 'gas_45_sa_fmin' / &
         SNTags(375)%Label / 'gas_45_sa_fmax' / &
         SNTags(376)%Label / 'gas_45_sa_hfn_fmin' / &
         SNTags(377)%Label / 'gas_45_sa_min_st' / &
         SNTags(378)%Label / 'gas_45_sa_min_un' / &
         SNTags(379)%Label / 'gas_45_sa_max' / &
         SNTags(380)%Label / 'gas_46_sa_fmin' / &
         SNTags(381)%Label / 'gas_46_sa_fmax' / &
         SNTags(382)%Label / 'gas_46_sa_hfn_fmin' / &
         SNTags(383)%Label / 'gas_46_sa_min_st' / &
         SNTags(384)%Label / 'gas_46_sa_min_un' / &
         SNTags(385)%Label / 'gas_46_sa_max' / &
         SNTags(386)%Label / 'gas_47_sa_fmin' / &
         SNTags(387)%Label / 'gas_47_sa_fmax' / &
         SNTags(388)%Label / 'gas_47_sa_hfn_fmin' / &
         SNTags(389)%Label / 'gas_47_sa_min_st' / &
         SNTags(390)%Label / 'gas_47_sa_min_un' / &
         SNTags(391)%Label / 'gas_47_sa_max' / &
         SNTags(392)%Label / 'gas_48_sa_fmin' / &
         SNTags(393)%Label / 'gas_48_sa_fmax' / &
         SNTags(394)%Label / 'gas_48_sa_hfn_fmin' / &
         SNTags(395)%Label / 'gas_48_sa_min_st' / &
         SNTags(396)%Label / 'gas_48_sa_min_un' / &
         SNTags(397)%Label / 'gas_48_sa_max' / &
         SNTags(398)%Label / 'gas_49_sa_fmin' / &
         SNTags(399)%Label / 'gas_49_sa_fmax' / &
         SNTags(400)%Label / 'gas_49_sa_hfn_fmin' /
    data SNTags(401)%Label / 'gas_49_sa_min_st' / &
         SNTags(402)%Label / 'gas_49_sa_min_un' / &
         SNTags(403)%Label / 'gas_49_sa_max' / &
         SNTags(404)%Label / 'gas_50_sa_fmin' / &
         SNTags(405)%Label / 'gas_50_sa_fmax' / &
         SNTags(406)%Label / 'gas_50_sa_hfn_fmin' / &
         SNTags(407)%Label / 'gas_50_sa_min_st' / &
         SNTags(408)%Label / 'gas_50_sa_min_un' / &
         SNTags(409)%Label / 'gas_50_sa_max' / &
         SNTags(410)%Label / 'gas_51_sa_fmin' / &
         SNTags(411)%Label / 'gas_51_sa_fmax' / &
         SNTags(412)%Label / 'gas_51_sa_hfn_fmin' / &
         SNTags(413)%Label / 'gas_51_sa_min_st' / &
         SNTags(414)%Label / 'gas_51_sa_min_un' / &
         SNTags(415)%Label / 'gas_51_sa_max' / &
         SNTags(416)%Label / 'gas_52_sa_fmin' / &
         SNTags(417)%Label / 'gas_52_sa_fmax' / &
         SNTags(418)%Label / 'gas_52_sa_hfn_fmin' / &
         SNTags(419)%Label / 'gas_52_sa_min_st' / &
         SNTags(420)%Label / 'gas_52_sa_min_un' / &
         SNTags(421)%Label / 'gas_52_sa_max' / &
         SNTags(422)%Label / 'gas_53_sa_fmin' / &
         SNTags(423)%Label / 'gas_53_sa_fmax' / &
         SNTags(424)%Label / 'gas_53_sa_hfn_fmin' / &
         SNTags(425)%Label / 'gas_53_sa_min_st' / &
         SNTags(426)%Label / 'gas_53_sa_min_un' / &
         SNTags(427)%Label / 'gas_53_sa_max' / &
         SNTags(428)%Label / 'gas_54_sa_fmin' / &
         SNTags(429)%Label / 'gas_54_sa_fmax' / &
         SNTags(430)%Label / 'gas_54_sa_hfn_fmin' / &
         SNTags(431)%Label / 'gas_54_sa_min_st' / &
         SNTags(432)%Label / 'gas_54_sa_min_un' / &
         SNTags(433)%Label / 'gas_54_sa_max' / &
         SNTags(434)%Label / 'gas_55_sa_fmin' / &
         SNTags(435)%Label / 'gas_55_sa_fmax' / &
         SNTags(436)%Label / 'gas_55_sa_hfn_fmin' / &
         SNTags(437)%Label / 'gas_55_sa_min_st' / &
         SNTags(438)%Label / 'gas_55_sa_min_un' / &
         SNTags(439)%Label / 'gas_55_sa_max' / &
         SNTags(440)%Label / 'gas_56_sa_fmin' / &
         SNTags(441)%Label / 'gas_56_sa_fmax' / &
         SNTags(442)%Label / 'gas_56_sa_hfn_fmin' / &
         SNTags(443)%Label / 'gas_56_sa_min_st' / &
         SNTags(444)%Label / 'gas_56_sa_min_un' / &
         SNTags(445)%Label / 'gas_56_sa_max' / &
         SNTags(446)%Label / 'gas_57_sa_fmin' / &
         SNTags(447)%Label / 'gas_57_sa_fmax' / &
         SNTags(448)%Label / 'gas_57_sa_hfn_fmin' / &
         SNTags(449)%Label / 'gas_57_sa_min_st' / &
         SNTags(450)%Label / 'gas_57_sa_min_un' / &
         SNTags(451)%Label / 'gas_57_sa_max' / &
         SNTags(452)%Label / 'gas_58_sa_fmin' / &
         SNTags(453)%Label / 'gas_58_sa_fmax' / &
         SNTags(454)%Label / 'gas_58_sa_hfn_fmin' / &
         SNTags(455)%Label / 'gas_58_sa_min_st' / &
         SNTags(456)%Label / 'gas_58_sa_min_un' / &
         SNTags(457)%Label / 'gas_58_sa_max' / &
         SNTags(458)%Label / 'gas_59_sa_fmin' / &
         SNTags(459)%Label / 'gas_59_sa_fmax' / &
         SNTags(460)%Label / 'gas_59_sa_hfn_fmin' / &
         SNTags(461)%Label / 'gas_59_sa_min_st' / &
         SNTags(462)%Label / 'gas_59_sa_min_un' / &
         SNTags(463)%Label / 'gas_59_sa_max' / &
         SNTags(464)%Label / 'gas_60_sa_fmin' / &
         SNTags(465)%Label / 'gas_60_sa_fmax' / &
         SNTags(466)%Label / 'gas_60_sa_hfn_fmin' / &
         SNTags(467)%Label / 'gas_60_sa_min_st' / &
         SNTags(468)%Label / 'gas_60_sa_min_un' / &
         SNTags(469)%Label / 'gas_60_sa_max' / &
         SNTags(470)%Label / 'gas_61_sa_fmin' / &
         SNTags(471)%Label / 'gas_61_sa_fmax' / &
         SNTags(472)%Label / 'gas_61_sa_hfn_fmin' / &
         SNTags(473)%Label / 'gas_61_sa_min_st' / &
         SNTags(474)%Label / 'gas_61_sa_min_un' / &
         SNTags(475)%Label / 'gas_61_sa_max' / &
         SNTags(476)%Label / 'gas_62_sa_fmin' / &
         SNTags(477)%Label / 'gas_62_sa_fmax' / &
         SNTags(478)%Label / 'gas_62_sa_hfn_fmin' / &
         SNTags(479)%Label / 'gas_62_sa_min_st' / &
         SNTags(480)%Label / 'gas_62_sa_min_un' / &
         SNTags(481)%Label / 'gas_62_sa_max' / &
         SNTags(482)%Label / 'gas_63_sa_fmin' / &
         SNTags(483)%Label / 'gas_63_sa_fmax' / &
         SNTags(484)%Label / 'gas_63_sa_hfn_fmin' / &
         SNTags(485)%Label / 'gas_63_sa_min_st' / &
         SNTags(486)%Label / 'gas_63_sa_min_un' / &
         SNTags(487)%Label / 'gas_63_sa_max' / &
         SNTags(488)%Label / 'gas_64_sa_fmin' / &
         SNTags(489)%Label / 'gas_64_sa_fmax' / &
         SNTags(490)%Label / 'gas_64_sa_hfn_fmin' / &
         SNTags(491)%Label / 'gas_64_sa_min_st' / &
         SNTags(492)%Label / 'gas_64_sa_min_un' / &
         SNTags(493)%Label / 'gas_64_sa_max' /
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
         SCTags(27)%Label  / 'automatic_spectra_config' /
    !> END GENERATED FCC.SCTags
end module m_fx_global_var
