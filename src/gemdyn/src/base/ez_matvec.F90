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

!** matvec - 3D Matrix-vector product subroutine

      subroutine ez_matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                             F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use metric
      use omp_timing
      use sol_mem
      use tdpack
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(inout) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(out) :: F_prod

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, k0, k0t, km, kp, n
      integer :: i0,in,j0,jn
      real(kind=REAL64)  :: r1, r2, aa1, aa2, aa3, bb1, bb2, bb3, cc1,&
                 S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15
      real(kind=REAL64), parameter :: half=0.5d0
      integer(kind=8) :: a,b
!
!     ---------------------------------------------------------------
!
      i0 = 1+pil_w ; in = l_ni-pil_e
      j0 = 1+pil_s ; jn = l_nj-pil_n

      do k= 1, l_nk
         do j= j0, jn
         do i= i0, in
            ext_q(i,j,k)= F_vector(i,j,k)
         end do
         end do
! Possibly for LAM configurations ??
         ext_q(i0-1,:,k) = ext_q(i0,:,k)
         ext_q(in+1,:,k) = ext_q(in,:,k)
         ext_q(:,j0-1,k) = ext_q(:,j0,k)
         ext_q(:,jn+1,k) = ext_q(:,jn,k)
      end do
      do j= j0, jn
      do i= i0, in
         r1=  (GVM%zmom_8(i,j,l_nk+1)-GVM%zmom_8(i,j,l_nk))&
             /(GVM%zmom_8(i,j,l_nk)-GVM%zmom_8(i,j,l_nk-1))
         r2= (GVM%zmom_8(i,j,1)-ver_z_8%m(0))/(GVM%zmom_8(i,j,2)-GVM%zmom_8(i,j,1))
         ext_q(i,j,l_nk+1)=  (1+r1)*F_vector(i,j,l_nk  ) &
                           -    r1 *F_vector(i,j,l_nk-1)
         ext_q(i,j,0)= (1+r2)*F_vector(i,j,1) &
                         -r2 *F_vector(i,j,2)
      end do
      end do

      if ( Grd_yinyang_L) then
         call yyg_xchng_hlt (ext_q, l_minx,l_maxx,l_miny,l_maxy, &
                    l_ni,l_nj, l_nk+2, .false., 'CUBIC', .true.)
      else
         call HLT_split (0, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( ext_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif

      do k= 0, l_nk+1
         do j= j0-1, jn+1
         do i= i0  , in+1
            Qdqdx(i-1,j,k)= (ext_q(i,j,k)-ext_q(i-1,j,k))*geomh_invDXM_8(j) ! on U points/Moment
         end do
         end do
         do j= j0  , jn+1
         do i= i0-1, in+1
            Qdqdy(i,j-1,k)= (ext_q(i,j,k)-ext_q(i,j-1,k))*geomh_invDY_8 ! on V points/Moment
         end do
         end do
         if (k>0) then
            do j= j0-1, jn+1
            do i= i0-1, in+1
               Qdqdz(i,j,k-1)= GVM%mc_iJz_8(i,j,k-1)*(ext_q(i,j,k)-ext_q(i,j,k-1))!*Ver_idz_8%t(k) ! on Thermo levels
            end do
            end do
            do j= j0-1, jn+1
            do i= i0-1, in+1
               Qqbz(i,j,k-1)= half*(ext_q(i,j,k-1)+ext_q(i,j,k)) ! on Thermo levels
            end do
            end do
         endif
      end do
      do k= 1, l_nk
         do j= j0-1, jn+1
         do i= i0  , in+1
            Qbarxz(i-1,j,k)= half*( (Qdqdz(i,j,k)+Qdqdz(i-1,j,k))*Ver_wp_8%m(k) + (Qdqdz(i,j,k-1)+Qdqdz(i-1,j,k-1))*Ver_wm_8%m(k) ) ! on U points/Moment
         end do
         end do
         do j= j0  , jn+1
         do i= i0-1, in+1
            Qbaryz(i,j-1,k)= half*((Qdqdz(i,j,k)+Qdqdz(i,j-1,k))*Ver_wp_8%m(k) + (Qdqdz(i,j,k-1)+Qdqdz(i,j-1,k-1))*Ver_wm_8%m(k) ) ! on V points/Moment
         end do
         end do            
         do j= j0-1, jn+1
         do i= i0-1, in+1
            Qbarz (i,j,k)= Qdqdz(i,j,k)*Ver_wp_8%m(k) + Qdqdz(i,j,k-1)*Ver_wm_8%m(k) ! on Moment levels
         end do
         end do
      end do

      do k= 1, l_nk
         do j= j0, jn
            do i= i0, in
               F_prod(i,j,k)= -gg_8*ext_q(i,j,k) &
               +geomh_invDX_8 (j)*( (Qdqdx(i,j,k)-Qdqdx(i-1,j,k)) - (GVM%mc_Jx_8(i,j,k)*Qbarxz(i,j,k)-GVM%mc_Jx_8(i-1,j,k)*Qbarxz(i-1,j,k))) &
               +geomh_invDYM_8(j)*( (geomh_cyv_8(j)*Qdqdy(i,j,k)-geomh_cyv_8(j-1)*Qdqdy(i,j-1,k)) &
                                   -(GVM%mc_Jy_8(i,j,k)*geomh_cyv_8(j)*Qbaryz(i,j,k)-GVM%mc_Jy_8(i,j-1,k)*geomh_cyv_8(j-1)*Qbaryz(i,j-1,k))) &
               +gama_8*Ver_idz_8%m(k)*((Qdqdz(i,j,k)-mu_8*Qqbz(i,j,k)) - (Qdqdz(i,j,k-1)-mu_8*Qqbz(i,j,k-1))) &
               -gama_8*epsi_8*(Qbarz(i,j,k)-mu_8*ext_q(i,j,k)) &
               +(half*(Qdqdx(i-1,j,k)+Qdqdx(i,j,k)) - GVM%mc_Jx_8(i,j,k)*Qbarz(i,j,k))*GVM%mc_Ix_8(i,j,k) &
               +(half*(Qdqdy(i,j-1,k)+Qdqdy(i,j,k)) - GVM%mc_Jy_8(i,j,k)*Qbarz(i,j,k))*GVM%mc_Iy_8(i,j,k) &
               +gama_8* (Qbarz(i,j,k)-mu_8*ext_q(i,j,k)) * GVM%mc_Iz_8(i,j,k)
            end do
         end do
      end do
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine ez_matvec
