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

      subroutine delQ ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
                        F_Qu,F_Qv,F_Qw, F_Qq, F_k0,F_kn )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use metric
      use sol_mem
      use ver
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_k0:F_kn), intent(INOUT) :: F_q
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,G_nk), intent(INOUT) :: F_Qu,F_Qv,F_Qw,F_Qq

      integer :: HLT_np, HLT_start, HLT_end
      integer :: i, j, k, nk, ii
      integer :: km1,km2,km3,kp1,kp2,kp3
      real(kind=REAL64) :: dqx(-2:3), dqy(-2:3), qbz, u, v,wp_8(G_nk), wm_8(G_nk)
      real(kind=REAL64), parameter :: one=1.d0,half=0.5d0
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,0:l_nk) :: dqdzx,dqdzy
!
!     ---------------------------------------------------------------
!
      if (Schm_POSO == 5) then
         call delQ5th ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
                        F_Qu,F_Qv,F_Qw, F_Qq, F_k0,F_kn )
      !   call delQ5 ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
      !                  F_Qu,F_Qv,F_Qw, F_Qq, F_k0,F_kn )
         return
      else if (Schm_POSO == 3) then
         call delQ3rd ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
                        F_Qu,F_Qv,F_Qw, F_Qq, F_k0,F_kn )
         return
      endif
      
      nk= F_kn-F_k0+1
      if ( Grd_yinyang_L) then
         call yyg_xchng_8 (F_q, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, nk, .false., 'CUBIC', .true.)
      else
         call HLT_split (F_k0, F_kn, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 ( F_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif
      do k=1,G_nk
         wm_8(k) = (Ver_dqdz_8(k)-Ver_z_8%m(k))/&
                   (Ver_dqdz_8(k)-Ver_dqdz_8(k-1))
         wp_8(k) = one-wm_8(k)
      end do
      do k=0,G_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         kp2=min(k+2,G_nk)
         kp3=min(k+3,G_nk)
         do j= 1,l_nj
            do i= 1,l_ni
               do ii=0,1
                  dqx(ii) = GVM%mc_iJz_8(i+ii,j,k)*(F_q(i+ii,j,k+1)-F_q(i+ii,j,k))
                  dqy(ii) = GVM%mc_iJz_8(i,j+ii,k)*(F_q(i,j+ii,k+1)-F_q(i,j+ii,k))
               end do
               dqdzx(i,j,k)= half*(dqx(1)+dqx(0))
               dqdzy(i,j,k)= half*(dqy(1)+dqy(0))
            end do
         end do
      end do
      call HLT_split (0, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 (dqdzx(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      call gem_xch_halo_8 (dqdzy(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      do k=1,G_nk
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         do j= 1, l_nj
            do i= 1, l_ni
               u = wp_8(k)*dqdzx(i,j,k)+wm_8(k)*dqdzx(i,j,k-1)
               v = wp_8(k)*dqdzy(i,j,k)+wm_8(k)*dqdzy(i,j,k-1)
               F_Qq(i,j,k) = GVM%mc_iJz_8(i,j,k)*(F_q(i,j,k+1)-F_q(i,j,k))
               qbz         = half*(F_q(i,j,k)+F_q(i,j,k+1))
               F_Qw(i,j,k)= F_Qq(i,j,k) - mu_8*qbz
               F_Qu(i,j,k)= (F_q(i+1,j,k)-F_q(i  ,j,k))*geomh_invDXM_8(j) &
                          - GVM%mc_Jx_8(i,j,k) * u
               F_Qv(i,j,k)= geomh_cyv_8(j)*(F_q(i,j+1,k)-F_q(i  ,j,k))*geomh_invDYM_8(j) &
                          - geomh_cyv_8(j)*GVM%mc_Jy_8(i,j,k) * v
            end do
         end do
      end do
      call HLT_split (1, G_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 (F_Qu(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      call gem_xch_halo_8 (F_Qv(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
!     
!     ---------------------------------------------------------------
!     
      return
!      include 'H5th_ope.inc'
      end subroutine delQ
