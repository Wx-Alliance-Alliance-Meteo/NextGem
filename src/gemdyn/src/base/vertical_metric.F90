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
      implicit none

      integer :: i,j,k
      integer :: HLT_start, HLT_end, local_np
      real, parameter :: one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------
!
      if (Grd_yinyang_L) then
         call yyg_xchng_hlt (orography, l_minx,l_maxx,l_miny,l_maxy, &
                             l_ni,l_nj, 6, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, 6, local_np, HLT_start, HLT_end)
         call gem_xch_halo (orography(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      endif
      
      call lvl_heights ( GVM%zmom_8, GVM%ztht_8, &
                         fis0, orols, l_minx,l_maxx,l_miny,l_maxy)
      call lvl_heights ( zmomu_8, zthtu_8, &
                         fis0u, orolsu, l_minx,l_maxx,l_miny,l_maxy)
      call lvl_heights ( zmomv_8, zthtv_8, &
                         fis0v, orolsv, l_minx,l_maxx,l_miny,l_maxy)

!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            GVM%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*GVM%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
!!$omp enddo
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               dgzm(i,j,k)=GVM%zmom_8(i,j,k)-GVM%zmom_8(i,j,k+1)
               dgzt(i,j,k)=GVM%ztht_8(i,j,k)-GVM%ztht_8(i,j,k+1)
            end do
         end do
      end do
!!$omp enddo

! Impossible code in an OMP parallel region
!!$!$omp single
!!$      err=0
!!$      if (minval(dgzm)<0. .or. minval(dgzt)<0. ) err=-1
!!$      call gem_error (err,'vertical_metric','Heights NOT monotonically decreasing from model top')
!!$      if (Lun_debug_L) then
!!$         call glbstat ( dgzm,'DGZM',"metric",l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
!!$                        1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,l_nk )
!!$         call glbstat ( dgzt,'DGZT',"metric",l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
!!$                        1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,l_nk )
!!$      endif
!!$!$omp end single

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
               GVM%mc_Jx_8 (i,j,k)= Hderiv8(GVM%zmom_8(i-2,j,k),GVM%zmom_8(i-1,j,k),GVM%zmom_8(i  ,j,k),&
                                                 GVM%zmom_8(i+1,j,k),GVM%zmom_8(i+2,j,k),GVM%zmom_8(i+3,j,k),geomh_invDX_8(j))
               GVM%mc_Jy_8 (i,j,k)= Hderiv8(GVM%zmom_8(i,j-2,k),GVM%zmom_8(i,j-1,k),GVM%zmom_8(i  ,j,k),&
                                                 GVM%zmom_8(i,j+1,k),GVM%zmom_8(i,j+2,k),GVM%zmom_8(i,j+3,k),geomh_invDY_8)

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

!     call heights_uv () !revisit that code
!
!     ---------------------------------------------------------------
!
      return
      include 'H5th_ope.inc'
      end subroutine vertical_metric
