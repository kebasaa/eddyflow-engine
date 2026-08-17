!***************************************************************************
! init_run_log.f90
! ----------------
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
! \brief       Open the run log, named like every other output of the run.
! \author      Jonathan Muller, ETH Zurich
! \note        The same composition InitOutFiles_rp uses for the full output -
!              project id, what the file holds, the run timestamp, the
!              extension - so the log sorts beside the results it describes
!              and carries the same run-mode tag.
!
!              Called once each binary knows two things: where its output goes,
!              and what this run's timestamp is. Everything said before then was
!              buffered by m_log and is flushed by LogInit.
!
!              RP and FCC each take their own timestamp at their own start, so
!              the two logs of one processing run never collide.
! \sa          m_log.f90
!***************************************************************************
subroutine InitRunLog()
    use m_common_global_var
    implicit none
    character(PathLen) :: LogPath

    call Clearstr(LogPath)
    LogPath = Dir%main_out(1:len_trim(Dir%main_out)) &
            // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
            // Log_FilePadding &
            // Timestamp_FilePadding(1:len_trim(Timestamp_FilePadding)) &
            // LogExt

    call LogInit(LogPath(1:len_trim(LogPath)))
end subroutine InitRunLog
