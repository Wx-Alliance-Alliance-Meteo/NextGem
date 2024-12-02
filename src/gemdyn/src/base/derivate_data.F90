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

!** s/r derivate_data

      subroutine derivate_data ( F_zd, F_w, F_u, F_v, F_t, F_s, F_q , &
                                 F_topo, F_orols, F_zmom_8, F_ztht_8, &
                      Minx, Maxx, Miny, Maxy,Nk, F_zd_L, F_w_L, F_q_L )
      use, intrinsic :: iso_fortran_env
      use cstv
      use dyn_fisl_options
      use dynkernel_options
      use gem_options
      use glb_ld
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      logical, intent(in) :: F_zd_L, F_w_L, F_q_L
      real, dimension(Minx:Maxx,Miny:Maxy   ), intent(in):: F_topo, F_orols
      real, dimension(Minx:Maxx,Miny:Maxy   ), intent(inout ):: F_s
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out)   :: F_zd, F_w,F_q
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v, F_t
      real(kind=REAL64), dimension(*), intent(in) :: F_zmom_8, F_ztht_8

      integer :: i,j,k
      integer :: HLT_start, HLT_end, local_np
!
!     ________________________________________________________________
!
      if (F_q_L) call diag_q ( F_q, F_s, F_t, F_topo, F_orols, &
                                   Minx, Maxx, Miny, Maxy, Nk)
      if ( (F_zd_L.or.F_w_L) ) then
         call HLT_split (1, G_nk, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( F_u(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call gem_xch_halo ( F_v(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call gem_xch_halo ( F_t(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call gem_xch_halo ( F_q(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call diag_zd_w( F_zd, F_w, F_u, F_v, F_t, F_q,&
                         F_zmom_8, F_ztht_8,Minx, Maxx, Miny, Maxy,&
                         Nk, F_zd_L, F_w_L)
      endif
!     
!     ________________________________________________________________
!
      return
      end subroutine derivate_data
