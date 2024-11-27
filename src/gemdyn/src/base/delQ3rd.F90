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

      subroutine delQ3rd ( F_q, F_minx,F_maxx,F_miny,F_maxy, &
                           F_Qu,F_Qv,F_Qw, F_Qq, F_k0,F_kn )
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
                               intent(INOUT) :: F_Qu,F_Qv,F_Qw,F_Qq

      integer :: HLT_np, HLT_start, HLT_end
      integer :: i, j, k, ii, km2,kp1,kp2
      real(kind=REAL64) :: dqx(-1:2), dqy(-1:2), qbz
      real(kind=REAL64) :: u(l_ni,l_nj),v(l_ni,l_nj)
      real(kind=REAL64),dimension(l_minx:l_maxx,l_miny:l_maxy,0:l_nk)::dqdzx,dqdzy
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!11111111 1 modifs
      k= 0
      do j= 1,l_nj
         do i= 1,l_ni
            do ii=-1,2
               dqx(ii) = GVM%mc_iJz_8(i+ii,j,k)*(F_q(i+ii,j,k+1)-F_q(i+ii,j,k))
               dqy(ii) = GVM%mc_iJz_8(i,j+ii,k)*(F_q(i,j+ii,k+1)-F_q(i,j+ii,k))
            end do
            dqdzx(i,j,k)= half*(dqx(1)+dqx(0))
            dqdzy(i,j,k)= half*(dqy(1)+dqy(0))
         end do
      end do
      do k=1,G_nk
         kp2=min(k+2,G_nk+1)
         do j= 1,l_nj
            do i= 1,l_ni
            do ii=-1,2
                  dqx(ii) = F_q(i+ii,j,k-1) * CDm2t(1,k) & 
                          + F_q(i+ii,j,k  ) * CDm2t(2,k) & 
                          + F_q(i+ii,j,k+1) * CDm2t(3,k) & 
                          + F_q(i+ii,j,kp2) * CDm2t(4,k)

                  dqy(ii) = F_q(i,j+ii,k-1) * CDm2t(1,k) & 
                          + F_q(i,j+ii,k  ) * CDm2t(2,k) & 
                          + F_q(i,j+ii,k+1) * CDm2t(3,k) & 
                          + F_q(i,j+ii,kp2) * CDm2t(4,k)
               end do
               dqdzx(i,j,k) = Hstag8(dqx(-1), dqx(0), dqx(1), dqx(2))
               dqdzy(i,j,k) = Hstag8(dqy(-1), dqy(0), dqy(1), dqy(2))
            end do
         end do
      end do

      do k=1,G_nk
         km2=max(k-2,0)
         kp1=min(k+1,G_nk)
         do j= 1, l_nj
            do i= 1, l_ni
               u(i,j) = dqdzx(i,j,km2) * CWt2m(1,k) &  
                      + dqdzx(i,j,k-1) * CWt2m(2,k) &  
                      + dqdzx(i,j,k  ) * CWt2m(3,k) &  
                      + dqdzx(i,j,kp1) * CWt2m(4,k)
               v(i,j) = dqdzy(i,j,km2) * CWt2m(1,k) &  
                      + dqdzy(i,j,k-1) * CWt2m(2,k) &  
                      + dqdzy(i,j,k  ) * CWt2m(3,k) &  
                      + dqdzy(i,j,kp1) * CWt2m(4,k)
               !--- remains second order for now ---
               F_Qq(i,j,k)= GVM%mc_iJz_8(i,j,k)*(F_q(i,j,k+1)-F_q(i,j,k))
               qbz        = half*(F_q(i,j,k)+F_q(i,j,k+1))
!!$               qbz = F_q(i,j,km1) * CWm2t(1,k) &
!!$                   + F_q(i,j,k  ) * CWm2t(2,k) &
!!$                   + F_q(i,j,kp1) * CWm2t(3,k) &
!!$                   + F_q(i,j,kp2) * CWm2t(4,k)
               F_Qw(i,j,k)= F_Qq(i,j,k) - mu_8*qbz
               F_Qu(i,j,k)= Hderiv8(F_q(i-1,j,k), F_q(i,j,k), &
                                    F_q(i+1,j,k), F_q(i+2,j,k), geomh_invDX_8(j)) &
                          - GVM%mc_Jx_8(i,j,k) * u(i,j)
       !        F_Qv(i,j,k)= Hderiv8(F_q(i,j-1,k), F_q(i,j,k)  , &
       !                             F_q(i,j+1,k), F_q(i,j+2,k), geomh_invDY_8) &
       !                           - GVM%mc_Jy_8(i,j,k) * v(i,j)
               F_Qv(i,j,k)= geomh_cyv_8(j)*Hderiv8(F_q(i,j-1,k), F_q(i,j,k)  , &
                                                   F_q(i,j+1,k), F_q(i,j+2,k), geomh_invDY_8) &
                                                   - geomh_cyv_8(j)*GVM%mc_Jy_8(i,j,k) * v(i,j)
            end do
         end do
         if ( .not. Grd_yinyang_L) then
         if (l_west) then
            i=pil_w
            do j= 1, l_nj
               F_Qu(  i,j,k)= -GVM%mc_Jx_8(i,j,k) * u(i,j)
               F_Qu(1+i,j,k)= (F_q(i+2,j,k)-F_q(i+1,j,k))*geomh_invDXM_8(j) - GVM%mc_Jx_8(i+1,j,k) * u(i+1,j)
            end do
         endif
         if (l_south) then
            j=pil_s
            do i= 1, l_ni
               F_Qv(i,  j,k)= -geomh_cyv_8(j)*GVM%mc_Jy_8(i,j,k) * v(i,j)
               F_Qv(i,1+j,k)=  geomh_cyv_8(j+1)*(F_q(i,j+2,k)-F_q(i,j+1,k))*geomh_invDYM_8(j+1) - geomh_cyv_8(j+1)*GVM%mc_Jy_8(i,j+1,k) * v(i,j+1)
            end do
         endif
         if (l_east) then
            i=l_ni-pil_e
            do j= 1, l_nj
               F_Qu(  i,j,k)= -GVM%mc_Jx_8(i,j,k) * u(i,j)
               F_Qu(i-1,j,k)= (F_q(i,j,k)-F_q(i-1,j,k))*geomh_invDXM_8(j) - GVM%mc_Jx_8(i-1,j,k) * u(i-1,j)
            end do
         endif
         if (l_north) then
            j=l_nj-pil_n
            do i= 1, l_ni
               F_Qv(i,  j,k)= -geomh_cyv_8(j)*GVM%mc_Jy_8(i,j,k) * v(i,j)
               F_Qv(i,j-1,k)=  geomh_cyv_8(j-1)*(F_q(i,j,k)-F_q(i,j-1,k))*geomh_invDYM_8(j-1) - geomh_cyv_8(j-1)*GVM%mc_Jy_8(i,j-1,k) * v(i,j-1)
            end do
         endif
         endif
      end do
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine delQ3rd

