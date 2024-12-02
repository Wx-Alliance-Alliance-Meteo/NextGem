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

!**s/r vertical_metric - Compute vertical metric coefficients

      subroutine vertical_metric ()
      use, intrinsic :: iso_fortran_env
      use gem_options
      use dyn_fisl_options
      use HORgrid_options
      use geomh
      use dcst
      use gmm_vt1
      use tdpack
      use glb_ld
      use metric
      use gmm_geof
      use ver
      use yyg_param
      implicit none

      integer :: i,j,k,km1,km2,km3,kp1,kp2,kp3
      integer :: HLT_start, HLT_end, local_np
      real(kind=REAL64) :: Jzu, Jzv, Jzq, Jx, Jy
      real, parameter :: one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------
!
      call heights ()

      if (Schm_POSO == 5) then
         call vertical_metric5th ( )
         return
      else if (Schm_POSO == 3) then
         call vertical_metric3rd ( )
         return
      endif
      
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            GVM%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*GVM%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do k=G_nk,1,-1
            do i=1-G_halox,l_ni+G_halox
               GVM%lg_pstar_8(i,j,k)=GVM%lg_pstar_8(i,j,k+1)+grav_8*(GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
         do i=1-G_halox,l_ni+G_halox
            GVM%ztht_8(i,j,G_nk)= GVM%zmom_8(i,j,G_nk+1) !temporary for mc_Ix_8 and mc_Iy_8 below
         end do
      end do
!!$omp enddo

!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               GVM%mc_Jx_8 (i,j,k)=(GVM%zmom_8(i+1,j,k)-GVM%zmom_8(i,j,k))*geomh_invDX_8(j)
               GVM%mc_Jy_8 (i,j,k)=(GVM%zmom_8(i,j+1,k)-GVM%zmom_8(i,j,k))*geomh_invDY_8
               GVM%mc_iJz_8(i,j,k)=one/(GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))
               GVM%mc_Ix_8(i,j,k)=log( (zthtu_8(i,j,k)-zthtu_8(i,j,k-1))/(zthtu_8(i-1,j,k)-zthtu_8(i-1,j,k-1)) )*geomh_invDX_8(j)
               GVM%mc_Iy_8(i,j,k)=log( (zthtv_8(i,j,k)-zthtv_8(i,j,k-1))/(zthtv_8(i,j-1,k)-zthtv_8(i,j-1,k-1)) )*geomh_invDY_8
               GVM%mc_Iz_8(i,j,k)=log( (GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                                  /(GVM%zmom_8(i,j,k)-GVM%zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
               GVM%mc_logJz_8(i,j,k)= 0.0
            end do
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
!DIR$ SIMD
         do i=1-G_halox+1,l_ni+G_halox-1
            GVM%mc_iJz_8(i,j,0)=one/(GVM%zmom_8(i,j,1)-ver_z_8%m(0))
            GVM%mc_css_H_8(i,j) = one/(gama_8*(GVM%mc_iJz_8(i,j,G_nk)-half*mu_8))
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
!DIR$ SIMD
         do i=1-G_halox,l_ni+G_halox
            GVM%ztht_8(i,j,G_nk)= ver_z_8%t(G_nk)+(Ver_b_8%t(G_nk)*fis0(i,j)+Ver_c_8%t(G_nk)*orols(i,j))/grav_8
         end do
      end do
!!$omp enddo

!!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            GVM%mc_alfas_H_8(i,j) = ( GVM%mc_iJz_8(i,j,G_nk) + half*mu_8 + Ver_wmstar_8(G_nk)*(GVM%mc_iJz_8(i,j,G_nk-1) -half*mu_8) ) / (GVM%mc_iJz_8(i,j,G_nk)-half*mu_8)
            GVM%mc_betas_H_8(i,j) =                                             Ver_wmstar_8(G_nk)*(GVM%mc_iJz_8(i,j,G_nk-1) +half*mu_8)   / (GVM%mc_iJz_8(i,j,G_nk)-half*mu_8)
         enddo
      enddo
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine vertical_metric
