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

    public :: ProcessSelfId, WatchParent, ParentGone, RequestFullSpeed
    public :: ProcessAlive, DiskFreeBytes

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

    !***************************************************************************
    !> \brief Is a process with this ID running now? Signal 0, as above.
    !***************************************************************************
    logical function ProcessAlive(pid)
        integer, intent(in) :: pid

        ProcessAlive = .false.
        if (pid <= 0) return
        ProcessAlive = p_kill(int(pid, c_int), 0_c_int) == 0_c_int
    end function ProcessAlive

    !***************************************************************************
    !> \brief Ask not to be run as background work - a no-op here.
    !>
    !> The Windows file explains what this is for. Linux and macOS schedule a
    !> windowless process like any other, so there is nothing to ask for.
    !***************************************************************************
    subroutine RequestFullSpeed()
    end subroutine RequestFullSpeed

    !***************************************************************************
    !> \brief How many bytes this user may still write to the disk holding
    !>        `dir`; -1 when that cannot be found out.
    !>
    !> Asked of df rather than statvfs: struct statvfs is laid out differently
    !> on Linux and macOS, and -P fixes df's columns on both. The available
    !> column is what an unprivileged user may still write, which is the
    !> figure that runs out. Read from the end of the line - the capacity
    !> column is the one ending in %, and available is just before it - so a
    !> file system name with a space in it does not shift the columns.
    !***************************************************************************
    integer(8) function DiskFreeBytes(dir)
        character(*), intent(in) :: dir
        character(64) :: scratch
        character(1024) :: line
        integer :: u
        integer :: io_status
        integer :: cmd_status
        integer :: i
        integer :: pct
        integer :: first
        integer(8) :: kb

        DiskFreeBytes = -1
        if (len_trim(dir) == 0) return
        write(scratch, '(a,i0,a)') '/tmp/eddyflow_df_', int(p_getpid()), '.txt'
        call execute_command_line('df -Pk "' // trim(dir) // '" > ' // trim(scratch) &
            // ' 2>/dev/null', exitstat = cmd_status)
        if (cmd_status /= 0) return
        open(newunit = u, file = trim(scratch), status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status /= 0) return
        read(u, '(a)', iostat = io_status)
        if (io_status == 0) read(u, '(a)', iostat = io_status) line
        close(u, status = 'delete')
        if (io_status /= 0) return

        pct = 0
        do i = 2, len_trim(line)
            if (line(i:i) == '%' .and. line(i - 1:i - 1) /= ' ') then
                pct = i
                exit
            end if
        end do
        if (pct == 0) return
        !> Back over the capacity token, the spaces before it, and then the
        !> available token itself.
        i = pct
        do while (i > 1 .and. line(i:i) /= ' ')
            i = i - 1
        end do
        do while (i > 1 .and. line(i:i) == ' ')
            i = i - 1
        end do
        first = i
        do while (first > 1 .and. line(first - 1:first - 1) /= ' ')
            first = first - 1
        end do
        read(line(first:i), *, iostat = io_status) kb
        if (io_status /= 0 .or. kb < 0) return
        DiskFreeBytes = kb * 1024_8
    end function DiskFreeBytes

end module m_process_os
