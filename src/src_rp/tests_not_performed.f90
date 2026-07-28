!***************************************************************************
! tests_not_performed.f90
! -----------------------
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
! \brief       Sets flags to 99999... for tests not performed
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestsNotPerformed()
    use m_rp_global_var
    implicit none

    !> The per-variable flags are strings now, so mark every digit as 9
    !> ("not available") rather than relying on the integer form. This also
    !> clears any value left over from the previous averaging period, which
    !> matters because a skipped test never writes its string.
    if(.not.Test%sr) CharHF%sr = repeat('9', FlagStrLen)
    if(.not.Test%ar) CharHF%ar = repeat('9', FlagStrLen)
    if(.not.Test%do) CharHF%do = repeat('9', FlagStrLen)
    if(.not.Test%al) CharHF%al = repeat('9', FlagStrLen)
    if(.not.Test%sk) CharHF%sk = repeat('9', FlagStrLen)
    if(.not.Test%sk) CharSF%sk = repeat('9', FlagStrLen)
    if(.not.Test%ds) CharHF%ds = repeat('9', FlagStrLen)
    if(.not.Test%ds) CharSF%ds = repeat('9', FlagStrLen)
    if(.not.Test%tl) IntHF%tl = 9999
    if(.not.Test%tl) IntSF%tl = 9999
    if(.not.Test%aa) IntHF%aa = 9
    if(.not.Test%ns) IntHF%ns = 9
end subroutine TestsNotPerformed
