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

!** delQ - Precompute Q hor. and vert. derivatives

      subroutine SW_delQ3rd ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
                              F_Qu,F_Qv,F_k0,F_kn )
      use geomh
      use HORgrid_options
      use glb_ld
      use metric
      use sol_mem
      use ver
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_k0:F_kn),&
                                                    intent(INOUT) :: F_q
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,G_nk)     ,&
                               intent(INOUT) :: F_Qu,F_Qv

      integer :: HLT_np, HLT_start, HLT_end
      integer :: i, j, k, nk, ii
      integer :: km1,km2,km3,kp1,kp2,kp3
      real(kind=REAL64) :: dqx(-1:2), dqy(-1:2), qbz, u, v
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy,0:l_nk)::dqdzx,dqdzy
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      nk= F_kn-F_k0+1
      if ( Grd_yinyang_L) then
         call yyg_xchng_8 (F_q(l_minx,l_miny,F_k0), YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, nk, .false., 'CUBIC', .true.)
      else
         call HLT_split (F_k0, F_kn, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 ( F_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif

      do k=1,G_nk
         do j= 1, l_nj
            do i= 1, l_ni
               F_Qu(i,j,k)= Hderiv8(F_q(i-1,j,k), F_q(i,j,k), &
                                    F_q(i+1,j,k), F_q(i+2,j,k), geomh_invDX_8(j)) 
               F_Qv(i,j,k)= Hderiv8(F_q(i,j-1,k), F_q(i,j,k)  , &
                                    F_q(i,j+1,k), F_q(i,j+2,k), geomh_invDY_8) 
              !F_Qu(i,j,k) = (F_q(i+1,j,k)-F_q(i,j,k))*geomh_invDX_8(j) 
              !F_Qv(i,j,k) = (F_q(i,j+1,k)-F_q(i,j,k))*geomh_invDYMv_8(j) 
            end do
         end do
      end do

      call HLT_split (1, G_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 (F_Qu(l_minx,l_miny,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      call gem_xch_halo_8 (F_Qv(l_minx,l_miny,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine SW_delQ3rd
