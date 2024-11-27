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
!**s/r prerhs2nd - Pre-compute stuff for elliptic_rhs

      subroutine prerhs2nd ( F_t2u, F_v2u, F_t2v, F_u2v, &
               F_dq2u, F_dq2v, F_dq2w, Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use gmm_vt0
      use mem_tstp
      use ver
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_dq2u, F_dq2v, F_dq2w

      integer :: i, j, k, km1, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(-0:1) :: t2qv, t2qu, v2q, u2q, dq2u, dq2v
      real(kind=REAL64) :: duu, dvv 
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) :: dqzv, dqzu
!
!     ---------------------------------------------------------------
!
      do k=1, l_nk
         km1=max(k-1,1)
         do j= 1, l_nj
           do i= 1, l_ni
             do n=0,1
                t2qv(n)= Ver_wp_8%m(k)*tt0(i+n,j,k)+Ver_wm_8%m(k)*tt0(i+n,j,km1)
                 v2q(n)= 0.5*(vt0(i,j-1,k)+vt0(i,j,k))
                t2qu(n)= Ver_wp_8%m(k)*tt0(i,j+n,k)+Ver_wm_8%m(k)*tt0(i,j+n,km1)
                 u2q(n)= 0.5*(ut0(i-1,j,k)+ut0(i,j,k))
                dq2u(n)= (qt0(i+n,j,k+1)-qt0(i+n,j,k))*GVM%mc_iJz_8(i+n,j,k)
                dq2v(n)= (qt0(i,j+n,k+1)-qt0(i,j+n,k))*GVM%mc_iJz_8(i,j+n,k)
             end do
             F_t2u(i,j,k)= 0.5*(t2qv(0)+t2qv(1))/Cstv_Tstr_8-1.d0
             F_v2u(i,j,k)= 0.5*( v2q(0)+ v2q(1))
             F_t2v(i,j,k)= 0.5*(t2qu(0)+t2qu(1))/Cstv_Tstr_8-1.d0
             F_u2v(i,j,k)= 0.5*( u2q(0)+ u2q(1))
             F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0) !for Rw
             dqzu(i,j,k)= 0.5*(dq2u(0)+dq2u(1)) !values in thermo, u-grid
             dqzv(i,j,k)= 0.5*(dq2v(0)+dq2v(1)) !values in thermo, v-grid
         end do
         end do
      end do

      !---now interpolate back to momentum level---
      !Note: dqdzu and dqdzv are already staggered appropriatly from previous loop
      do k=1,l_nk
         km1=max(k-1,1)
         do j= 1, l_nj 
            do i= 1, l_ni
               duu = Ver_wp_8%m(k)*dqzu(i,j,k)+Ver_wm_8%m(k)*dqzu(i,j,km1)
               dvv = Ver_wp_8%m(k)*dqzv(i,j,k)+Ver_wm_8%m(k)*dqzv(i,j,km1)
               F_dq2u(i,j,k)= duu * GVM%mc_Jx_8(i,j,k)
               F_dq2v(i,j,k)= dvv * GVM%mc_Jy_8(i,j,k)
            end do
         end do
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( F_t2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_v2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_t2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_u2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2w(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine prerhs2nd
