!***************************************************************************
! m_process_os.f90  (Linux and macOS)
! -----------------------------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
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
! \brief       What this program needs from the operating system about
!              processes: its own identity, and whether the process that
!              started it is still alive. See the Windows file of the same
!              name, src_os/win/, for why.
!
! \note        Must declare the same module and public names as the Windows
!              file; static_checks/test_worker_outlives_parent_static.py
!              holds the two to it.
!
!              Weaker than Windows in one way: there is no handle to hold, so
!              an ID freed by the parent and reused by an unrelated process
!              would read as the parent still running. That needs the ID to
!              come round again within one period, and costs only a worker
!              running to the end of its slice - the behaviour before this.
!
!              pid_t is a C int on both Linux and macOS.
! \author      Jonathan Muller
! \sa          prepass_parallel.f90, init_env.f90
!***************************************************************************
module m_process_os
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    private

    public :: ProcessSelfId, WatchParent, ParentGone

    integer, save :: parent_pid = 0

    interface
        !> Named apart from GNU Fortran's GETPID and KILL extension
        !> intrinsics, which -fall-intrinsics makes visible.
        function p_getpid() bind(C, name = 'getpid')
            import :: c_int
            integer(c_int) :: p_getpid
        end function p_getpid

        function p_kill(pid, sig) bind(C, name = 'kill')
            import :: c_int
            integer(c_int), value :: pid
            integer(c_int), value :: sig
            integer(c_int) :: p_kill
        end function p_kill
    end interface

contains

    !***************************************************************************
    !> \brief This process's ID, as a worker's parent hands it on.
    !***************************************************************************
    integer function ProcessSelfId()
        ProcessSelfId = int(p_getpid())
    end function ProcessSelfId

    !***************************************************************************
    !> \brief Start watching the process with this ID; .false. if it is gone.
    !***************************************************************************
    logical function WatchParent(pid)
        integer, intent(in) :: pid

        WatchParent = .false.
        if (pid <= 0) return
        parent_pid = pid
        WatchParent = .not. ParentGone()
        if (.not. WatchParent) parent_pid = 0
    end function WatchParent

    !***************************************************************************
    !> \brief Has the watched process exited?
    !>
    !> Signal 0 delivers nothing; it only asks whether the process exists.
    !> .false. when nothing is being watched.
    !***************************************************************************
    logical function ParentGone()
        ParentGone = .false.
        if (parent_pid <= 0) return
        ParentGone = p_kill(int(parent_pid, c_int), 0_c_int) /= 0_c_int
    end function ParentGone

end module m_process_os
