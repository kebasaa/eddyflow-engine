!***************************************************************************
! exception_handler.f90
! ---------------------
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
! \brief       Manages error and Warning instances, possibly aborting execution
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExceptionHandler(error_code)
    use m_index_parameters
    use m_log
    implicit none
    !> in/out variables
    integer :: error_code


    select case (error_code)
        case(0)
            call LogSayList(' Fatal error(0)> Occurred while retrieving user home path.')
            call LogSayList(' Fatal error(0)> EddyFlow is not able to locate project files and must terminate.')
            call LogSayList(' Fatal error(0)> Program execution aborted.')
            stop 1
        case(1)
            call LogSayList(' Fatal error(1)> Temporary file "flist.tmp" not created.')
            call LogSayList(' Fatal error(1)> Most likely, no files matching the raw file name format were found.')
            call LogSayList(' Fatal error(1)> Program execution aborted.')
            stop 1
        case(2)
            call LogSayList(' Error(2)> Occurred while scanning file for biomet data.')
            call LogSayList(' Error(2)> EddyFlow could not retrieve biomet data from this file.')
            call LogSayList(' Error(2)> File will be ignored.')
        case(3)
            call LogSayList(' Error(3)> Occurred while reading metadata file in GHG archive. File not found or empty.')
            call LogSayList(' Error(3)> GHG archive appears to be corrupted or invalid and will thus be skipped.')
        case(4)
            call LogSayList(' Error(4)> Occurred while reading raw data file in GHG archive. File not found or empty.')
            call LogSayList(' Error(4)> GHG archive appears to be corrupted or invalid and will thus be skipped.')
        case(5)
            call LogSayList(' Error(5)> Occurred while reading biomet data file or corresponding metadata file in GHG archive.')
            call LogSayList(' Error(5)> File(s) not found or empty. Biomet data not used for this time period.')
        case(6)
            call LogSayList(' Error(6)> Occurred while opening current raw data file. File is empty.')
            call LogSayList(' Error(6)> Raw data file skipped.')
        case(7)
            write(*,*)
            write(ulog,*)
            call LogSayList(' Error(7)> Occurred while opening INI-format file. Looking for a solution..')
        case(8)
            call LogSayList(' Fatal error(8)> No files matching the specified extension were found in the selected folder.')
            call LogSayList(' Fatal error(8)> Program execution aborted.')
            stop 1
        case(14)
            call LogSayList(' Error(14)> Occurred while unzipping GHG archive. File skipped.')
        case(20)
            call LogSayList('  Fatal error(20)> Encountered while interpreting a timestamp template.')
            call LogSayList('  Fatal error(20)> It may apply to retrieval of timestamps from the "Raw file name format"')
            call LogSayList('  Fatal error(20)> or the biomet or dynamic metadata records."')
            call LogSayList('  Fatal error(20)> Program execution aborted.')
            stop 1
        case(21)
            call LogSayList('  Fatal error(21)> Configuration file (project file with extension *.eddyflow) not found.')
            call LogSayList('  Fatal error(21)> Program execution aborted.')
            stop 1
        case(22)
            call LogSayList('  Fatal error(22)> "Alternative metadata file" not found.')
            call LogSayList('  Fatal error(22)> Program execution aborted.')
            stop 1
        case(23)
            call LogSayList('  Fatal error(23)> Occurred while validating "Alternative metadata file".')
            call LogSayList('  Fatal error(23)> Check entries in the "Metadata file editor" and try again.')
            call LogSayList('  Fatal error(23)> Program execution aborted.')
            stop 1
        case(24)
            call LogSayList('  Error(24)> Occurred while opening GHG file.')
            call LogSayList('  Error(24)> File skipped.')
        case(25)
            call LogSayList('  Error(25)> Occurred while validating embedded metadata file.')
            call LogSayList('  Error(25)> GHG file skipped.')
        case(28)
            call LogSayList('  Error(28)> Occurred while opening raw data file.')
            call LogSayList('  Error(28)> Raw file skipped.')
        case(29)
            call LogSayList(' Alert(29)> Occurred while reading rotation matrices from auxiliary planar-fit file.')
            call LogSayList(' Alert(29)> Axis rotation method switched to "Double rotation".')
        case(30)
            call LogSayList(' Alert(30)> Occurred while opening auxiliary file for planar-fit. File not found or empty.')
            call LogSayList(' Alert(30)> Axis rotation method switched to "Double rotation".')
        case(31)
            call LogSayList(' Fatal error(31)> Unrecognized TOB1 data format. Select either IEEE4 or FP2 data formats and try&
                & again.')
            call LogSayList(' Fatal error(31)> Program execution aborted.')
            stop 1
        case(32)
            call LogSayList(' Fatal error(32)> No valid GHG file was found in the selected "Raw data directory".')
            call LogSayList(' Fatal error(32)> The problem could also be due to corrupted archives or invalid metadata files.')
            call LogSayList(' Fatal error(32)> Program execution aborted.')
            stop 1
        case(33)
            call LogSayList(' Error(33)> Number of wind records for this sector is less than requested.')
            call LogSayList(' Error(33)> Planar-fit rotation matrix not calculated for this sector.')
        case(34)
            write(*,*) ' Error(34)> Occurred while calculating planar-fit rotations &
                                     &for this sector: singular matrix found.'
            write(ulog,*) ' Error(34)> Occurred while calculating planar-fit rotations &
                                     &for this sector: singular matrix found.'
            call LogSayList(' Error(34)> Planar-fit rotation matrix not calculated for this sector.')
        case(35)
            call LogSayList(' Fatal error(35)> Oops! Something went wrong. EddyFlow was not able to process any raw file.')
            call LogSayList(' Fatal error(35)> Output files not created.')
            call LogSayList(' Fatal error(35)> Program execution aborted.')
            stop 1
        case(36)
            call LogSayList(' Fatal error(36)> No "Output directory" was selected. Select an "Output directory" before running&
                & EddyFlow.')
            call LogSayList(' Fatal error(36)> Program execution aborted.')
            stop 1
        case(37)
            call LogSayList(' Alert(37)> No valid planar-fit rotation matrix found for any wind sector.')
            call LogSayList(' Alert(37)> Axis rotation method switched to "Double rotation".')
        case(38)
            call LogSayList(' Alert(38)> No sectors selected for planar fit.')
            call LogSayList(' Alert(38)> Forcing to 1 sector of 360 degrees.')
        case(39)
            call LogSayList(' Alert(39)> Error while opening auxiliary file for time-lag optimization. File not found or empty.')
            call LogSayList(' Alert(39)> Time-lag detection method switched to "Covariance maximization".')
        case(40)
            call LogSayList(' Alert(40)> Occurred while retrieving wind sectors configuration for planar-fit.')
            call LogSayList(' Alert(40)> Forcing to 1 sector of 360 degrees.')
        case(41)
            call LogSayList(' Warning(41)> Wind-sector excluded by user. No planar-fit matrix calculated for this sector.')
        case(42)
            call LogSayList(' Error(42)> Method for random uncertainty estimation not recognized.')
            call LogSayList(' Error(42)> Random uncertainty not calculated.')
        case(43)
            call LogSayList(' Alert(43)> Time-lag optimization failed.')
            call LogSayList(' Alert(43)> Switching to method "Covariance maximization" for time-lag detection.')
        case(44)
            call LogSayList(' Error(44)> Occurred while reading or interpreting biomet file.')
            call LogSayList(' Error(44)> Biomet data not used for this period.')
        case(45)
            call LogSayList(' Error(45)> Not enough valid co-spectra were found for fitting models, or fitting procedure failed.')
            call LogSayList(' Error(45)> Stability-sorted ensemble averaged cospectra outputs not created.')
        case(46)
            write(*,*) '  Fatal error(46)> The dataset does not contain any raw file &
                                           &corresponding to the selected sub-period.'
            write(ulog,*) '  Fatal error(46)> The dataset does not contain any raw file &
                                           &corresponding to the selected sub-period.'
            call LogSayList('  Fatal error(46)> Select another sub-period or a different "Raw data directory".')
            call LogSayList('  Fatal error(46)> Try also un-checking the option "Select sub-period" or')
            call LogSayList('  Fatal error(46)> checking the option "Search in sub-folders", in the "Basic settings page".')
            call LogSayList('  Fatal error(46)> Program execution aborted.')
            stop 1
        case(48)
            write(*,*) '  Fatal error(48)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for planar-fit.'
            write(ulog,*) '  Fatal error(48)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for planar-fit.'
            call LogSayList('  Fatal error(48)> Select another sub-period or another "Raw data directory".')
            call LogSayList('  Fatal error(48)> Try also un-checking the option "Select sub-period" in the planar-fit dialogue.')
            call LogSayList('  Fatal error(48)> Program execution aborted.')
            stop 1
        case(49)
            write(*,*) '  Fatal error(49)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for time-lag optimization.'
            write(ulog,*) '  Fatal error(49)> The dataset does not contain any raw file &
                                           &corresponding to the sub-period selected for time-lag optimization.'
            call LogSayList('  Fatal error(49)> Select another sub-period or another "Raw data directory".')
            write(*,*) '  Fatal error(49)> Try also un-checking the option "Select sub-period" &
                                           &in the time-lag optimization dialogue.'
            write(ulog,*) '  Fatal error(49)> Try also un-checking the option "Select sub-period" &
                                           &in the time-lag optimization dialogue.'
            call LogSayList('  Fatal error(49)> Program execution aborted.')
            stop 1
        case(50)
            write(*,*) '  Fatal error(50)> "Essentials" files does not contain any results &
                                           &corresponding to the selected sub-period.'
            write(ulog,*) '  Fatal error(50)> "Essentials" files does not contain any results &
                                           &corresponding to the selected sub-period.'
            call LogSayList('  Fatal error(50)> Select another sub-period.')
            call LogSayList('  Fatal error(50)> Try also un-checking the option "Select sub-period" in "Basic settings" page.')
            call LogSayList('  Fatal error(50)> Program execution aborted.')
            stop 1
        case(51)
            call LogSayList('  Error(51)> Calculation of tube attenuation parameter (lambda) failed for at least one gas.')
            call LogSayList('  Error(51)> Processing continues but tube attenuation is not included in the spectral correction.')
        case(52)
            call LogSayList('  Fatal error(52)> An unexpected, unrecognized internal problem occurred.')
            call LogSayList('  Fatal error(52)> Program execution aborted.')
            stop 1
        case(53)
            call LogSayList('  Warning(53)> No raw data file relevant to current averaging period was found.')
            call LogSayList('  Warning(53)> Skipping to next averaging period.')
        case(54)
            call LogSayList('  Error(54)> Occurred while opening or reading SLT (binary) file.')
            call LogSayList('  Error(54)> Raw file skipped.')
        case(55)
            call LogSayList('  Error(55)> Occurred while opening or reading generic binary file.')
            call LogSayList('  Error(55)> Raw file skipped.')
        case(56)
            call LogSayList('  Error(56)> Occurred while opening or reading TOB1 (binary) file.')
            call LogSayList('  Error(56)> Raw file skipped.')
        case(57)
            call LogSayList('  Error(57)> Occurred while opening or reading ASCII file.')
            call LogSayList('  Error(57)> Raw file skipped.')
        case(58)
            call LogSayList(' Warning(58)> Available samples not enough for an averaging period.')
            call LogSayList(' Warning(58)> Skipping to next averaging period.')
        case(59)
            call LogSayList('  Error(59)> At least one wind component appears to be corrupted (too many implausible values).')
            call LogSayList('  Error(59)> This may also be the result of data exclusion by the "Absolute limits" test or by a')
            call LogSayList('  Error(59)> custom-designed "Flag" in the "Basic Settings" page.')
            call LogSayList('  Error(59)> If the problem occurs for many or all raw files, check those settings.')
            call LogSayList('  Error(59)> Skipping to next averaging period.')
        case(60)
            call LogSayList(' Fatal error(60)> Occurred while opening or reading "essentials" file.')
            call LogSayList(' Fatal error(60)> Execution will be aborted, but you may be able to avoid re-processing raw data')
            call LogSayList(' Fatal error(60)> by fixing the problem with the "essentials" file and using the option')
            call LogSayList(' Fatal error(60)> "Previous data directory" in the "Basic Settings" page.')
            call LogSayList(' Fatal error(60)> Please consult software documentation.')
            call LogSayList(' Fatal error(60)> Program execution aborted.')
            stop 1
        case(61)
            call LogSayList(' Fatal error(61)> No valid data records found in the "essentials" file.')
            call LogSayList(' Fatal error(61)> Program execution aborted.')
            stop 1
        case(62)
            call LogSayList(' Error(62)> Occurred while reading binned (co)spectra file.')
            call LogSayList(' Error(62)> File skipped.')
        case(63)
            call LogSayList(' Error(63)> Occurred while reading full (co)spectra file.')
            call LogSayList(' Error(63)> File skipped.')
        case(64)
            call LogSayList(' Error(64)> Occurred while creating output file.')
            call LogSayList(' Error(64)> Some spectral assessment results will not be written on output file.')
        case(65)
            call LogSayList(' Alert(65)> Occurred while reading auxiliary "spectral assessment" file.')
            call LogSayList(' Alert(65)> High-frequency spectral correction method switched to Moncrieff et al. (1997).')
        case(66)
            call LogSayList(' Alert(66)> Acquisition frequency appears to be set to a value <= zero.')
            call LogSayList(' Alert(66)> High-frequency spectral corrections cannot be calculated.')
            call LogSayList(' Alert(66)> Proceeding without spectral corrections.')
        case(67)
            call LogSayList(' Error(67)> Occurred while opening output file.')
            call LogSayList(' Error(67)> Output dataset not created.')
        case(68)
            call LogSayList(' Error(68)> Occurred while reading dynamic metadata file. File not found or empty.')
            call LogSayList(' Error(68)> Dynamic metadata not used in this run.')
        case(69)
            call LogSayList(' Error(69)> There is a problem with results of the spectral assessment.')
            call LogSayList(' Error(69)> High-frequency spectral correction method switched to Moncrieff et al. (1997).')
        case(70)
            call LogSayList(' Error(70)> Inconsistent number of variables in biomet files.')
            call LogSayList(' Error(70)> EddyFlow cannot resolve the conflict and will thus proceed without using biomet data.')
        case(71)
            call LogSayList(' Error(71)> No valid biomet record imported.')
            call LogSayList(' Error(71)> EddyFlow will proceed without using biomet data.')
        case(72)
            call LogSayList('  Warning(72)> No valid biomet record found for this period.')
        case(73)
            call LogSayList('  Error(73)> The label of at least one biomet variable misses')
            call LogSayList('  Error(73)> or has incomplete positional qualifier ("_x_y_z" suffix).')
            call LogSayList('  Error(73)> EddyFlow will proceed without using biomet data.')
        case(74)
            call LogSayList('  Error(74)> No valid binned (co)spectra files were found')
            call LogSayList('  Error(74)> EddyFlow cannot perform spectral asssessment, nor ')
            call LogSayList('  Error(74)> create ensemble averaged (co)spectra. If the case, spectral ')
            call LogSayList('  Error(74)> correction method will be switched to Moncrieff et al. (2007).')
        case(75)
            call LogSayList(' Error(75)> Not enough valid co-spectra were found for making ensemble averages.')
            call LogSayList(' Error(75)> Time-sorted ensemble averaged cospectra outputs not created.')
        case(76)
            call LogSayList(' Error(76)> EddyFlow could not calculate ensemble spectra.')
            call LogSayList(' Error(76)> Spectral assessment failed. Spectral assessment file not created.')
        case(77)
            call LogSayList(' Error(77)> EddyFlow could not calculate ensemble spectra.')
            call LogSayList(' Error(77)> Ensemble averaged spectral output not created.')
        case(78)
            call LogSayList(' Fatal error(78)> No files matching the expected template were found in the selected folder.')
            call LogSayList(' Fatal error(78)> Program execution aborted.')
            stop 1
        case(79)
            call LogSayList(' Error(79)> Inconsistent variable labels or units in biomet files.')
            call LogSayList(' Error(79)> EddyFlow cannot resolve the conflict and will thus proceed without using biomet data.')
        case(80)
            call LogSayList(' Warning(80)> Implausible altitude value detected. Altitude defaulted to zero.')
        case(81)
            call LogSayList(' Warning(81)> Implausible latitude value detected. Latitude defaulted to zero.')
        case(82)
            call LogSayList(' Warning(82)> Implausible longitude value detected. Longitude defaulted to zero.')
        case(83)
            call LogSayList(' Warning(83)> Implausible canopy height value detected. Canopy height defaulted to zero.')
        case(84)
            call LogSayList(' Warning(84)> Implausible displacement height value detected.')
            call LogSayList(' Warning(84)> Displacement height defaulted to 0.67 times the canopy height.')
        case(85)
            call LogSayList(' Warning(85)> Implausible roughness length value detected.')
            call LogSayList(' Warning(85)> Roughness length  defaulted to 0.15 times the canopy height,')
            call LogSayList(' Warning(85)> or to 1.0mm if canopy height is zero.')
        case(86)
            call LogSayList(' Fatal error(86)> Could not retrieve files from directory. Either directory does not exist')
            call LogSayList(' Fatal error(86)> or it does not contain files matching the selected requirements.')
            call LogSayList(' Fatal error(86)> Program execution aborted.')
            stop 1
        case(87)
            call LogSayList('  Error(87)> Entered or inferred "Binned co-spectra files directory" does not exist.')
            call LogSayList('  Error(87)> EddyFlow cannot perform spectral assessment, calculate ensemble averaged spectra')
            call LogSayList('  Error(87)> or calculate ensemble averaged co-spectra.')
            call LogSayList('  Error(87)> Continuing by switching to Moncrieff et al. (1997) spectral corrections')
            call LogSayList('  Error(87)> if an in-situ method was selected, and ignoring selection of ensemble')
            call LogSayList('  Error(87)> averaged spectra or co-spectra outputs')
        case(88)
            call LogSayList('  Error(88)> Entered or inferred "Full co-spectra files directory" does not exist.')
            call LogSayList('  Error(88)> EddyFlow cannot use spectral correction method of Fratini et al. (2012)')
            call LogSayList('  Error(88)> Continuing by switching to Moncrieff et al. (1997)')
        case(89)
            call LogSayList('  Error(89)> Entered or inferred "Binned co-spectra files directory" does not contain any valid&
                & files.')
            call LogSayList('  Error(89)> EddyFlow cannot perform spectral assessment, calculate ensemble averaged spectra')
            call LogSayList('  Error(89)> or calculate ensemble averaged co-spectra.')
            call LogSayList('  Error(89)> Continuing by switching to Moncrieff et al. (1997) spectral corrections')
            call LogSayList('  Error(89)> if an in-situ method was selected, and ignoring selection of ensemble')
            call LogSayList('  Error(89)> averaged spectra or co-spectra outputs')
        case(90)
            call LogSayList('  Error(90)> Entered or inferred "Binned co-spectra files directory" does not contain')
            call LogSayList('  Error(90)> any files corresponding to the selected start/end period.')
            call LogSayList('  Error(90)> EddyFlow cannot perform spectral assessment, calculate ensemble averaged spectra')
            call LogSayList('  Error(90)> or calculate ensemble averaged co-spectra.')
            call LogSayList('  Error(90)> EddyFlow will continue, switching to Moncrieff et al. (1997) spectral corrections')
            call LogSayList('  Error(90)> if an in-situ method was selected, and ignoring selection of ensemble')
            call LogSayList('  Error(90)> averaged spectra or co-spectra outputs.')
        case(91)
            call LogSayList('  Alert(91)> Time constant of linear detrending cannot be larger than flux averaging interval.')
            call LogSayList('  Alert(91)> Automatically set time constant equal to flux averaging interval.')
        case(92)
            call LogSayList(' Error(92)> Occurred while writing passive gases ensemble spectra on output file.')
            call LogSayList(' Error(92)> File not created.')
        case(93)
            call LogSayList(' Error(93)> Embedded biomet selected, but not processing GHG files. This is not possible.')
            call LogSayList(' Error(93)> EddyFlow will proceed without using biomet data.')
        case(94)
            call LogSayList(' Warning(94)> The selected or inferred angle-of-attack correction method is not applicable')
            call LogSayList(' Warning(94)> to data collected with selected sonic anemometer.')
            call LogSayList(' Warning(94)> Continuing without applying any angle-of-attack correction.')
        case(95)
            call LogSayList(' Warning(95)> The selected "w-boost" correction is not applicable')
            call LogSayList(' Warning(95)> to data collected with selected sonic anemometer.')
            call LogSayList(' Warning(95)> Continuing without applying "w-boost" correction.')
        case(96)
            call LogSayList('  Fatal error(96)> The project file (*.eddyflow) was written by a newer')
            call LogSayList('  Fatal error(96)> version of the EddyFlow GUI than this engine can read.')
            call LogSayList('  Fatal error(96)> Reading it would silently drop or misinterpret settings.')
            call LogSayList('  Fatal error(96)> Update the EddyFlow engine, or re-create the project')
            call LogSayList('  Fatal error(96)> with a matching version of the GUI.')
            call LogSayList('  Fatal error(96)> Program execution aborted.')
            stop 1
        case(97)
            call LogSayList(' Warning(97)> The metadata file describes more instruments than this engine')
            call LogSayList(' Warning(97)> can hold. The surplus instruments have been ignored, and any')
            call LogSayList(' Warning(97)> variables assigned to them will not be used.')
            call LogSayList(' Warning(97)> Reduce the number of instruments in the "Metadata file editor".')
        case(98)
            call LogSayList(' Warning(98)> An INI-format file holds more entries in one section than')
            call LogSayList(' Warning(98)> EddyFlow can store. The surplus entries were ignored, so some')
            call LogSayList(' Warning(98)> settings may have fallen back to their defaults.')
            call LogSayList(' Warning(98)> This points to an unusually large project or metadata file.')
        case(99)
            call LogSayList('  Fatal error(99)> The project file (*.eddyflow) states no gas count.')
            call LogSayList('  Fatal error(99)> Gases, cell measurements and diagnostics are described')
            call LogSayList('  Fatal error(99)> by indexed records from format version 5.0.0 onward, and')
            call LogSayList('  Fatal error(99)> a file predating them would process no gas fluxes at all.')
            call LogSayList('  Fatal error(99)> This is not the same as a project that states gas_num=0:')
            call LogSayList('  Fatal error(99)> a site with an anemometer and no analyser is valid and is')
            call LogSayList('  Fatal error(99)> processed for sensible heat flux.')
            call LogSayList('  Fatal error(99)> Open the project in the EddyFlow GUI and save it - an')
            call LogSayList('  Fatal error(99)> older file is upgraded automatically on opening.')
            call LogSayList('  Fatal error(99)> Program execution aborted.')
            stop 1
        case(100)
            call LogSayList(' Warning(100)> A gas record names a species EddyFlow has no molecular')
            call LogSayList(' Warning(100)> weight or diffusivity for, and carries neither of its own.')
            call LogSayList(' Warning(100)> Nitrous oxide values were used instead. A molecular weight')
            call LogSayList(' Warning(100)> that is wrong by a factor gives a flux wrong by the same')
            call LogSayList(' Warning(100)> factor, and nothing else in the output will look unusual.')
            call LogSayList(' Warning(100)> Set the molecular weight and diffusivity for that gas in')
            call LogSayList(' Warning(100)> the EddyFlow GUI.')
        !> 101 was the fixed full-output format's "gases past the fourth are
        !> dropped" warning. That format is gone and the code is not reused:
        !> a message a user may have seen must not come back meaning
        !> something else.
        case(102)
            call LogSayList(' Warning(102)> A gas states a spectral-assessment month grouping that')
            call LogSayList(' Warning(102)> cannot be read. It must be a list of month ranges, such as')
            call LogSayList(' Warning(102)> 1-12 for one group over the calendar or 1-6,7-12 for two,')
            call LogSayList(' Warning(102)> with every bound between 1 and 12 and no month in more')
            call LogSayList(' Warning(102)> than one group.')
            call LogSayList(' Warning(102)> That gas was given one group spanning the calendar, which')
            call LogSayList(' Warning(102)> is what it would have had if it stated nothing at all.')
        case(103)
            call LogSayList(' Warning(103)> The spectral assessment file gives different transfer')
            call LogSayList(' Warning(103)> function parameters to months this project pools into one')
            call LogSayList(' Warning(103)> group, so it was fitted under a different month grouping')
            call LogSayList(' Warning(103)> than the project now declares.')
            call LogSayList(' Warning(103)> Each group took the parameters of its first month. Re-run')
            call LogSayList(' Warning(103)> the assessment if the grouping was meant to change.')
        case(104)
            call LogSayList(' Warning(104)> This project measures gases but has no humidity: no gas')
            call LogSayList(' Warning(104)> record names H2O, and no biomet relative humidity column')
            call LogSayList(' Warning(104)> is selected. Two things follow, and neither is visible in')
            call LogSayList(' Warning(104)> the output.')
            call LogSayList(' Warning(104)> Air density and heat capacity are computed for DRY air, so')
            call LogSayList(' Warning(104)> density is overestimated and every density-based correction')
            call LogSayList(' Warning(104)> carries that bias.')
            call LogSayList(' Warning(104)> And the humidity correction to sensible heat flux cannot be')
            call LogSayList(' Warning(104)> applied, so the reported H is the uncorrected buoyancy flux,')
            call LogSayList(' Warning(104)> not a true sensible heat flux. Over a wet surface the two')
            call LogSayList(' Warning(104)> differ by several percent.')
            call LogSayList(' Warning(104)> A biomet relative humidity sensor is enough to remove both.')
        case(105)
            call LogSayList(' Warning(105)> A row of the dynamic metadata or calibration file has more')
            call LogSayList(' Warning(105)> fields than this engine can hold, and the ones past the')
            call LogSayList(' Warning(105)> limit were dropped. Those files carry a block per gas, so')
            call LogSayList(' Warning(105)> their width grows with the number of gases in the project.')
            call LogSayList(' Warning(105)> Reduce the number of columns, or split the file.')
        case(106)
            call LogSayList(' Warning(106)> The gas named on the line above is corrected with a')
            call LogSayList(' Warning(106)> hygrometer on a different analyser, because its own')
            call LogSayList(' Warning(106)> carries none. The water-flux term of its WPL correction')
            call LogSayList(' Warning(106)> is taken at that hygrometer''s own time lag, and the')
            call LogSayList(' Warning(106)> dilution to a mixing ratio uses that cell''s humidity')
            call LogSayList(' Warning(106)> rather than its own - a stand-in, not a measurement of')
            call LogSayList(' Warning(106)> the air this gas was sampled in.')
            call LogSayList(' Warning(106)> Declare an H2O column on the gas''s own analyser to')
            call LogSayList(' Warning(106)> remove the compromise: EddyFlow prefers a hygrometer on')
            call LogSayList(' Warning(106)> the gas''s own instrument whenever the project has one.')

        case(107)
            call LogSayList(' Fatal error(107)> The "essentials" file was written by an earlier')
            call LogSayList(' Fatal error(107)> version of EddyFlow-RP and cannot be read.')
            call LogSayList(' Fatal error(107)> Its per-gas water vapour records are narrower than')
            call LogSayList(' Fatal error(107)> this version expects. Reading them would not fail -')
            call LogSayList(' Fatal error(107)> it would run on into the next gas''s fields and mix')
            call LogSayList(' Fatal error(107)> one gas''s water terms into another''s, silently.')
            call LogSayList(' Fatal error(107)> Re-run EddyFlow-RP to regenerate the file, then run')
            call LogSayList(' Fatal error(107)> EddyFlow-FCC again. There is nothing to migrate.')
            call LogSayList(' Fatal error(107)> Program execution aborted.')
            stop 1

        case(108)
            call LogSayList(' Warning(108)> This project describes no gas analyser, so no gas')
            call LogSayList(' Warning(108)> flux can be computed. EddyFlow will process it as an')
            call LogSayList(' Warning(108)> anemometer-only site.')
            call LogSayList(' Warning(108)> Computed: momentum flux and u*, sensible heat from the')
            call LogSayList(' Warning(108)> sonic temperature, Monin-Obukhov length and stability,')
            call LogSayList(' Warning(108)> wind speed and direction, and the turbulence statistics.')
            call LogSayList(' Warning(108)> Not computed: every gas flux, LE, ET, and the WPL and')
            call LogSayList(' Warning(108)> spectral corrections that depend on them. Their columns')
            call LogSayList(' Warning(108)> are written as the error label throughout - the FLUXNET')
            call LogSayList(' Warning(108)> format requires CO2, H2O and CH4 columns to exist even')
            call LogSayList(' Warning(108)> where nothing measured them.')
            call LogSayList(' Warning(108)> The sensible heat flux carries no humidity correction')
            call LogSayList(' Warning(108)> either, there being no hygrometer to supply one.')
        case(109)
            call LogSayList(' Warning(109)> One or more gases state no absolute limits, so the')
            call LogSayList(' Warning(109)> absolute-limits test was not performed on them and')
            call LogSayList(' Warning(109)> their flag digit is 9. Their data is used as it stands:')
            call LogSayList(' Warning(109)> nothing checks whether it is physically plausible, and a')
            call LogSayList(' Warning(109)> fill value the file does not declare would reach the flux.')
            call LogSayList(' Warning(109)> Set a minimum and a maximum for every gas under')
            call LogSayList(' Warning(109)> Advanced > Statistical Analysis > Absolute limits.')
    end select
end subroutine ExceptionHandler
