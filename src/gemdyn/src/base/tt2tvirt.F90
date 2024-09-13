!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------
!s/r tt2tvirt - Compute virtual temperature from temperature+wload

      subroutine tt2tvirt ( F_tv, F_tt, Minx, Maxx, Miny, Maxy, &
                            Nktv, Nktt, F_i0,F_in,F_j0,F_nj )
      use dyn_fisl_options
      use glb_ld
      use tr3d
      use mem_tracers
      implicit none

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nktv, Nktt, F_i0,F_in,F_j0,F_nj
      real, dimension(Minx:Maxx, Miny:Maxy, Nktt), intent(in ) :: F_tt
      real, dimension(Minx:Maxx, Miny:Maxy, Nktv), intent(out) :: F_tv
!
!     ________________________________________________________________
!
      call sumhydro ( sumq_8, l_minx,l_maxx,l_miny,l_maxy, &
                      G_nk, Tr3d_ntr, trt1, Schm_wload_L )

      call mfottvh2 (F_tt, F_tv, tracers_P(Tr3d_hu)%pntr, real(sumq_8),&
           minx, maxx, miny, maxy, l_nk, F_i0,F_in,F_j0,F_nj, .true.)
!
!     ________________________________________________________________
!
      return
      end subroutine tt2tvirt
