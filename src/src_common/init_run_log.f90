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
    if (BatchIndex > 0) then
        !> A worker of a parallel pre-pass logs beside the records it was
        !> asked for, not into the output directory. It shares its parent's
        !> start second, so it would otherwise open the parent's own log and
        !> whichever process closed it last would win; and the output
        !> directory belongs to the run, not to an internal detail of how the
        !> pre-pass was computed. The parent folds this file into the run log.
        LogPath = BatchOutPath(1:len_trim(BatchOutPath)) // LogExt
    else
        LogPath = Dir%main_out(1:len_trim(Dir%main_out)) &
                // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                // Log_FilePadding &
                // Timestamp_FilePadding(1:len_trim(Timestamp_FilePadding)) &
                // LogExt
    end if

    call LogInit(LogPath(1:len_trim(LogPath)))
end subroutine InitRunLog
