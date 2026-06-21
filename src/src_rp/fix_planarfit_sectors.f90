!***************************************************************************
! fix_planarfit_sectors.f90
! -------------------------
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
! \brief       Replaces error planar fits with closest valid sector rotation matrix
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixPlanarfitSectors(GoPlanarFit, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(inout) :: GoPlanarFit(N)
    !> Local variables
    integer :: sec
    integer :: sec2
    real(kind = dbl) :: loc_pfmat(3, 3, -N + 1: 2 * N)
    logical :: loc_go(-N + 1: 2 * N)
    real(kind = dbl)  :: loc_pfb(3, -N + 1: 2 * N)


    !> First, if there is no valid sector, switches to 2D rotations
    do sec = 1, N
        if (GoPlanarFit(sec)) exit
    end do
    if (sec == N + 1) then
        Meth%rot = 'double_rotation'
        call ExceptionHandler(37)
        return
    end if

    !> Define working arrays, only for convenience, for going
    !> clockwise or counter-clockwise easier
    !>
    !> N-1           0   1            N   N+1          2*N
    !>  |--o--x--*---|   |--o--x--*---|    |--o--x--*---|
    !>
    loc_pfmat(:, :, -N + 1: 0)    = PFMat(:, :, 1:N)
    loc_pfmat(:, :, 1: N)         = PFMat(:, :, 1:N)
    loc_pfmat(:, :, N + 1: 2 * N) = PFMat(:, :, 1:N)
    loc_go(-N + 1: 0)    = GoPlanarFit(1:N)
    loc_go(1: N)         = GoPlanarFit(1:N)
    loc_go(N + 1: 2 * N) = GoPlanarFit(1:N)
    loc_pfb(:, -N + 1: 0)    = PFb(:, 1:N)
    loc_pfb(:, 1: N)         = PFb(:, 1:N)
    loc_pfb(:, N + 1: 2 * N) = PFb(:, 1:N)

    do sec = 1, N
        if (.not. GoPlanarFit(sec)) then
            if(PFSetup%fix == 'clockwise') then
                !> Searches clockwise
                do sec2 = sec + 1, 2*N
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        PFb(:, sec) = loc_pfb(:, sec2)
                        GoPlanarFit(Sec) = .true.
                        exit
                    end if
                end do
            elseif(PFSetup%fix == 'counterclockwise') then
                !> Searches counterclockwise
                do sec2 = sec - 1, - N + 1
                    if (loc_go(sec2)) then
                        PFMat(:, :, sec) = loc_pfmat(:, :, sec2)
                        PFb(:, sec) = loc_pfb(:, sec2)
                        GoPlanarFit(sec) = .true.
                        exit
                    end if
                end do
            end if
        end if
    end do
end subroutine FixPlanarfitSectors
