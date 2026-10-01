!***************************************************************************
! m_process_os.f90  (Windows)
! ---------------------------
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
!              started it is still alive.
!
!              A pre-pass worker is a separate process, started through a
!              shell script, and nothing ties its life to its parent's.
!              Killing the parent - the interface's Stop button, Task
!              Manager, an error stop while it waits - used to leave every
!              worker running for hours on a slice nobody would ever merge.
!              Measured on the Yatir run: the parent gone, six workers and
!              their six launchers still computing. So a worker now watches
!              its parent and stops when it goes.
!
! \note        One file per platform - this one and src_os/posix/ - chosen by
!              the Makefile. Both declare the same module with the same public
!              names; static_checks/test_worker_outlives_parent_static.py holds
!              them to it, and gen_makefile_deps.py reads only this one.
!
!              64-bit Windows has a single calling convention, so the WINAPI
!              entry points bind directly. DWORD is 32 bits there (LLP64), and
!              a HANDLE is pointer-sized.
! \author      Jonathan Muller
! \sa          prepass_parallel.f90, init_env.f90
!***************************************************************************
module m_process_os
    use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_intptr_t
    implicit none
    private

    public :: ProcessSelfId, WatchParent, ParentGone

    !> SYNCHRONIZE, the only right WaitForSingleObject needs. Asking for no
    !> more means it cannot be refused on a process the same user started.
    integer(c_int32_t), parameter :: SYNCHRONIZE_RIGHT = 1048576_c_int32_t
    integer(c_int32_t), parameter :: WAIT_OBJECT_0 = 0_c_int32_t

    !> Held for the whole life of the worker. Beyond being what the wait is
    !> done on, an open handle keeps the parent's process object alive after it
    !> exits - and with it its ID, which Windows will not hand to a new process
    !> while a handle remains. So a recycled ID cannot pass for the parent.
    integer(c_intptr_t), save :: parent_handle = 0_c_intptr_t

    interface
        function w_GetCurrentProcessId() bind(C, name = 'GetCurrentProcessId')
            import :: c_int32_t
            integer(c_int32_t) :: w_GetCurrentProcessId
        end function w_GetCurrentProcessId

        function w_OpenProcess(access, inherit, pid) bind(C, name = 'OpenProcess')
            import :: c_int32_t, c_int, c_intptr_t
            integer(c_int32_t), value :: access
            integer(c_int), value :: inherit
            integer(c_int32_t), value :: pid
            integer(c_intptr_t) :: w_OpenProcess
        end function w_OpenProcess

        function w_WaitForSingleObject(h, ms) bind(C, name = 'WaitForSingleObject')
            import :: c_intptr_t, c_int32_t
            integer(c_intptr_t), value :: h
            integer(c_int32_t), value :: ms
            integer(c_int32_t) :: w_WaitForSingleObject
        end function w_WaitForSingleObject
    end interface

contains

    !***************************************************************************
    !> \brief This process's ID, as a worker's parent hands it on.
    !***************************************************************************
    integer function ProcessSelfId()
        ProcessSelfId = int(w_GetCurrentProcessId())
    end function ProcessSelfId

    !***************************************************************************
    !> \brief Start watching the process with this ID.
    !>
    !> .false. when it cannot be opened, which for a parent means it has already
    !> gone. Opening is done once, at start-up: an ID opened later could belong
    !> to someone else by then.
    !***************************************************************************
    logical function WatchParent(pid)
        integer, intent(in) :: pid

        WatchParent = .false.
        if (pid <= 0) return
        parent_handle = w_OpenProcess(SYNCHRONIZE_RIGHT, 0_c_int, int(pid, c_int32_t))
        WatchParent = parent_handle /= 0_c_intptr_t
    end function WatchParent

    !***************************************************************************
    !> \brief Has the watched process exited?
    !>
    !> A zero-timeout wait, so it costs nothing. .false. when nothing is being
    !> watched, which leaves a process that was never told its parent - every
    !> run that is not a worker - exactly as it was.
    !***************************************************************************
    logical function ParentGone()
        ParentGone = .false.
        if (parent_handle == 0_c_intptr_t) return
        ParentGone = w_WaitForSingleObject(parent_handle, 0_c_int32_t) == WAIT_OBJECT_0
    end function ParentGone

end module m_process_os
