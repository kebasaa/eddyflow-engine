!***************************************************************************
! init_env.f90
! ------------
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
! \brief       Initialize environment variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitEnv()
    use m_common_global_var
    use m_eddypro_import
    implicit none
    include 'version_and_date.inc'
    !> local variables
    integer :: i
    integer :: ch
    integer :: make_dir
    integer :: io_status
    integer :: aux
    character(32) :: timestring
    character(256) :: switch
    character(PathLen) :: projPath
    character(256) :: arg
    character(32) :: tmpDirPadding
    character(3), parameter :: OS_default = 'win'
    integer, external :: CreateDir
    logical, external :: SwitchTakesValue
    character(PathLen) :: lowerPath
    character(8) :: batchPadding
    integer :: jobs_status


    !> Store current timestamp information
    call hms_current_string(timestring)
    if(timestring(12:12) == ' ') timestring(12:12) = '0'
    Timestamp_FilePadding = '_' // timestring(1:10) // 'T' &
        // timestring(12:13) // timestring(15:16) // timestring(18:19)

    tmpDirPadding = '_' // timestring(1:10) // 'T' &
        // timestring(12:13) // timestring(15:16) // timestring(18:19) // timestring(21:23)

    !> Get command lines switches and values
    OS = ''
    homedir = ''
    NumJobs = 0
    BatchIndex = 0
    BatchCount = 0
    BatchKind = ''
    BatchOutPath = ''
    EddyFlowProj%run_env = ''
    EddyFlowProj%caller = ''
    projPath = ''
    i = 1
    arg_loop: do
        !> Read switch
        call get_command_argument(i, value=switch, status=io_status)
        if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
        i = i + 1

        !> Only a switch that takes a value consumes the token after it.
        !> This used to be unconditional, so the project path swallowed
        !> whatever followed it: a -e written after the path was read as the
        !> path's own value and silently dropped, and the token after that was
        !> then examined as though it were a switch. Nothing but the path
        !> could usefully go last, and nothing said so.
        arg = ''
        io_status = 0
        if (SwitchTakesValue(switch)) then
            call get_command_argument(i, value=arg, status=io_status)
            i = i + 1
        end if

        if (switch(1:1) == '-') then
            select case(trim(adjustl(switch)))

                !> Switch for "system", the host operating system
                case('-s', '--system')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    OS = trim(arg)
                    if (OS(1:1) == '-') OS = ''

                !> Switch for "environment", the 'home' working directory
                case('-e', '--environment')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    homedir = trim(arg)
                    if (homedir(1:1) == '-') homedir = ''

                !> Switch for "mode", whether "embedded" or "desktop" mode
                case('-m', '--mode')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    EddyFlowProj%run_env = trim(arg)
                    if (EddyFlowProj%run_env(1:1) == '-') EddyFlowProj%run_env = ''

                !> Switch for "caller", whether "gui" or "console"
                case('-c', '--caller')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    EddyFlowProj%caller = trim(arg)
                    if (EddyFlowProj%caller(1:1) == '-') EddyFlowProj%caller = ''

                !> Switch for "jobs", how many worker processes the
                !> assessment pre-passes may be split across. 0, and the
                !> absent switch, mean "as many as there are cores"; 1 is
                !> the serial path. Anything unparseable is treated as
                !> absent rather than as an error, in keeping with every
                !> other switch here.
                case('-j', '--jobs')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    if (arg(1:1) /= '-') then
                        read(arg, *, iostat = jobs_status) NumJobs
                        if (jobs_status /= 0 .or. NumJobs < 0) NumJobs = 0
                    end if

                !> Switches for "batch", which this program passes to copies
                !> of itself and a user never types. --batch is
                !> <kind>:<k>:<n>:<first>:<last>, naming the pre-pass ('pf' or
                !> 'to'), which slice of it this process owns, and the period
                !> indices that slice spans; --batch-out is where to leave the
                !> resulting records. Deliberately absent from
                !> CommandLineHelp: they are an implementation detail of the
                !> parallel pre-pass, not an interface.
                case('--batch')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    call ParseBatchArg(arg)


                case('--batch-out')
                    if (io_status > 0 .or. len_trim(switch) == 0) exit arg_loop
                    BatchOutPath = trim(arg)

                !> Software version
                case('-v', '--version')
                    call InformOfSoftwareVersion(sw_ver, build_date)

                !> Minimal command line help
                case('-h', '--help')
                    call CommandLineHelp(sw_ver, build_date)
            end select
        else
            projPath = trim(switch)
            lowerPath = projPath
            !> Its own counter. This loop used to run on `i`, the argument
            !> loop's counter, and left it at len_trim(projPath) + 1 - so the
            !> next iteration read past the end of the command line and every
            !> switch written AFTER the project path was silently dropped.
            do ch = 1, len_trim(lowerPath)
                if (iachar(lowerPath(ch:ch)) >= 65 .and. iachar(lowerPath(ch:ch)) <= 90) &
                    lowerPath(ch:ch) = achar(iachar(lowerPath(ch:ch)) + 32)
            end do
            !> An EddyPro project is accepted here so that it reaches
            !> ImportEddyProProject below. Without this it is discarded
            !> silently and the run falls back to ini/processing.eddyflow -
            !> which either does not exist, or is somebody else's project.
            if (index(lowerPath, '.eddyflow') == 0 .and. &
                index(lowerPath, '.eddypro') == 0) projPath = ''
        end if
    end do arg_loop

    !> Set OS-dependent parameters
    if (len_trim(OS) == 0) OS = OS_default
    call SetOSEnvironment()

    !> A batch worker shares its parent's start second, so in embedded mode
    !> it would land in the parent's temporary directory - where the raw
    !> file list is built under a fixed name, and where concurrent workers
    !> would overwrite each other's. The batch index makes the directory
    !> its own in both modes.
    !>
    !> Timestamp_FilePadding deliberately does NOT get the same treatment:
    !> it is 22 characters wide and most callers concatenate it whole,
    !> trailing blanks and all, so widening it would put spaces in every
    !> output file name. A worker writes nothing into the output directory
    !> except its run log, and InitRunLog tags that itself.
    if (BatchIndex > 0) then
        write(batchPadding, '(a,i2.2)') '_b', BatchIndex
        tmpDirPadding = trim(tmpDirPadding) // trim(batchPadding)
    end if

    !> Default values if args are not passed
    if (len_trim(homedir) == 0) homedir = '..'
    if (len_trim(EddyFlowProj%run_env) == 0) EddyFlowProj%run_env = 'desktop'
    if (len_trim(EddyFlowProj%caller) == 0)  EddyFlowProj%caller  = 'console'

    !> Define default unit number (udf), run specific
    call hms_current_hms(aux, aux, aux, udf)
    if (udf < 200) udf = udf + 200
    udf2 = udf + 1

    !> Define path of key EddyFlow files/dirs
    call AdjDir(homedir, slash)
    IniDir = trim(homedir) // 'ini' // slash
    if (projPath == '') then
        PrjPath = trim(IniDir) // trim(PrjFile)
    else
        PrjPath = projPath
    end if

    !> An EddyPro project is imported into our own format, once, beside the
    !> file it came from, and PrjPath is repointed at the result. Everything
    !> downstream then reads an ordinary project and knows nothing of this.
    !>
    !> Called on whatever path was resolved rather than only on a .eddypro
    !> one, so a project renamed to our extension is recognised too: it
    !> decides by the file's own first line and returns untouched for
    !> anything that is not an EddyPro project.
    call ImportEddyProProject(PrjPath, sw_ver)

    !> Define TmpDir differently if it's in desktop or embedded mode
    if (EddyFlowProj%run_env == 'desktop') then
        TmpDir = trim(homedir) // 'tmp' // slash // 'tmp' &
        // trim(adjustl(tmpDirPadding)) // slash
    else
        TmpDir = trim(homedir) // 'tmp' // slash
    end if

    !> Create TmpDir in case it doesn't exist (for use from command line)
    make_dir = CreateDir('"' // trim(TmpDir) // '"')
end subroutine InitEnv

!***************************************************************************
!
! rief       Whether a command line switch is followed by its own value.
! \note        -v and -h are not, and neither is anything that is not a
!              switch at all - the project path above all. Every switch that
!              IS listed here must also have a case in InitEnv's select, or
!              its value is consumed and then ignored.
!***************************************************************************
logical function SwitchTakesValue(switch)
    implicit none
    character(*), intent(in) :: switch
    character(64) :: token

    token = adjustl(switch)
    select case (trim(token))
        case ('-s', '--system', '-e', '--environment', '-m', '--mode', &
              '-c', '--caller', '-j', '--jobs', &
              '--batch', '--batch-out')
            SwitchTakesValue = .true.
        case default
            SwitchTakesValue = .false.
    end select
end function SwitchTakesValue

!***************************************************************************
!
! \brief       Split the --batch argument into its five colon-separated
!              fields: <kind>:<index>:<count>:<first>:<last>.
! \note        A malformed value is fatal rather than ignored, unlike the
!              user-facing switches. Nobody types this one - it is written by
!              the parent process - so a value that does not parse means the
!              two halves of this program disagree, and continuing would
!              silently pre-process the wrong periods.
!***************************************************************************
subroutine ParseBatchArg(arg)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: arg
    integer :: i
    integer :: from
    integer :: sepa
    integer :: io_status
    integer :: field(4)
    character(64) :: token

    sepa = index(arg, ':')
    if (sepa <= 1) error stop 'Malformed --batch argument.'
    BatchKind = arg(1:sepa - 1)

    from = sepa + 1
    do i = 1, 4
        if (i < 4) then
            sepa = index(arg(from:len_trim(arg)), ':')
            if (sepa <= 1) error stop 'Malformed --batch argument.'
            token = arg(from: from + sepa - 2)
            from = from + sepa
        else
            token = arg(from: len_trim(arg))
        end if
        read(token, *, iostat = io_status) field(i)
        if (io_status /= 0) error stop 'Malformed --batch argument.'
    end do

    BatchIndex = field(1)
    BatchCount = field(2)
    BatchSliceStart = field(3)
    BatchSliceEnd = field(4)
    if (BatchIndex < 1 .or. BatchCount < 1 .or. BatchIndex > BatchCount &
        .or. BatchSliceStart < 1 .or. BatchSliceEnd < BatchSliceStart) &
        error stop 'Malformed --batch argument.'
end subroutine ParseBatchArg

!***************************************************************************
! \file        src/init_env.f90
! \brief       Provide software version on demand (switch "-v")
! \version     4.2.0
! \date        2013-07-02
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InformOfSoftwareVersion(sw_ver, build_date)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: sw_ver
    character(*), intent(in) :: build_date


    write (*, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    write(ulog, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    stop
end subroutine InformOfSoftwareVersion

!***************************************************************************
!
! \brief       Provide software version on demand (switch "-v")
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CommandLineHelp(sw_ver, build_date)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: sw_ver
    character(*), intent(in) :: build_date
    !> Local variables
    character(12) :: prog


    if (app == 'EddyFlow-RP') then
        prog = 'eddyflow_rp'
    else
        prog = 'eddyflow_fcc'
    end if

    write(*, '(a)') ' Help for ' // trim(adjustl(app))
    write(ulog, '(a)') ' Help for ' // trim(adjustl(app))
    call LogSay(' --------------------')
    write (*, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    write(ulog, '(a)') ' ' // trim(adjustl(app)) // ', version ' // trim(adjustl(sw_ver)) // &
        &', build ' // trim(adjustl(build_date)) // '.'
    write(*,*)
    write(ulog,*)
    write(*, '(a)') ' USAGE: ' // trim(prog) // ' [OPTION [ARG]] [PROJ_FILE]'
    write(ulog, '(a)') ' USAGE: ' // trim(prog) // ' [OPTION [ARG]] [PROJ_FILE]'
    write(*,*)
    write(ulog,*)
    call LogSay(' OPTIONS:')
    call LogSay('   [-s | --system [win | linux | mac]]  Operating system; if not provided assumes "win"')
    call LogSay('   [-m | --mode [embedded | desktop]]   Running mode; if not provided assumes "desktop"')
    call LogSay('   [-c | --caller [gui | console]]      Caller; if not provided assumes "console"')
    write(*, '(a)') '   [-e | --environment [DIRECTORY]]     Working directory, to be provided in embedded mode;&
                                                             & if not provided assumes \.'
    write(ulog, '(a)') '   [-e | --environment [DIRECTORY]]     Working directory, to be provided in embedded mode;&
                                                             & if not provided assumes \.'
    call LogSay('   [-j | --jobs [N]]                    Worker processes for the planar-fit and time-lag&
                                                         & pre-passes; 0 or absent uses every core, 1 is serial')
    call LogSay('   [-h | --help]                        Display this help and exit')
    call LogSay('   [-v | --version]                     Output version information and exit')
    write(*, '(a)')
    write(ulog, '(a)')
    write(*, '(a)') ' PROJ_FILE                              Path of project (*.eddyflow) file;&
                                                             & if not provided, assumes ..\ini\processing.eddyflow'
    write(ulog, '(a)') ' PROJ_FILE                              Path of project (*.eddyflow) file;&
                                                             & if not provided, assumes ..\ini\processing.eddyflow'
    stop
end subroutine CommandLineHelp

