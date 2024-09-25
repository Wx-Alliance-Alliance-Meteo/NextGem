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

      subroutine vertical_metric3rd ()
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
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            GVM%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*GVM%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do

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

      do k=1,G_nk
         km1=max(k-1,0)
         km2=max(k-2,0)
         km3=max(k-3,0)
         kp1=min(k+1,G_nk+1)
         kp2=min(k+2,G_nk+1)
         kp3=min(k+3,G_nk+1)
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               Jzu = zthtu_8(i,j,km2) * CDt2m(1,k)&
                    +zthtu_8(i,j,km1) * CDt2m(2,k)&
                    +zthtu_8(i,j,k  ) * CDt2m(3,k)&
                    +zthtu_8(i,j,kp1) * CDt2m(4,k)
               Jzv = zthtv_8(i,j,km2) * CDt2m(1,k)&
                    +zthtv_8(i,j,km1) * CDt2m(2,k)&
                    +zthtv_8(i,j,k  ) * CDt2m(3,k)&
                    +zthtv_8(i,j,kp1) * CDt2m(4,k)
               Jzq = GVM%zmom_8(i,j,km1) * CDm2t(1,k)&
                    +GVM%zmom_8(i,j,k  ) * CDm2t(2,k)&
                    +GVM%zmom_8(i,j,kp1) * CDm2t(3,k)&
                    +GVM%zmom_8(i,j,kp2) * CDm2t(4,k)
               Jx= Hderiv8(GVM%zmom_8(i-1,j,k),GVM%zmom_8(i  ,j,k),&
                           GVM%zmom_8(i+1,j,k),GVM%zmom_8(i+2,j,k),geomh_invDX_8(j))
               Jy= Hderiv8(GVM%zmom_8(i,j-1,k),GVM%zmom_8(i  ,j,k),&
                           GVM%zmom_8(i,j+1,k),GVM%zmom_8(i,j+2,k),geomh_invDY_8)
               M_Jxozu(i,j,k)= Jx  / Jzu ! First  term in matvec portion of eqn 58
               M_Jyozv(i,j,k)= Jy  / Jzv ! Second term in matvec portion of eqn 58
               M_iJzq (i,j,k)= one / Jzq ! Third  term in matvec portion of eqn 58
               M_Jzu(i,j,k) = log(Jzu)
               M_Jzv(i,j,k) = log(Jzv)
               M_Jzq(i,j,k) = log(Jzq)
            end do
         end do
      end do
      if (Grd_yinyang_L) then
         call yyg_xchng_8 (M_Jzu, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk+2, .false., 'CUBIC', .true.)
         call yyg_xchng_8 (M_Jzv, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk+2, .false., 'CUBIC', .true.)
         call yyg_xchng_8 (M_Jzq, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk+2, .false., 'CUBIC', .true.)
      else
         call HLT_split (0, 3*(G_nk+1)+2, local_np, HLT_start, HLT_end)
         call gem_xch_halo_8 (M_Jzu(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      endif
      do k=1,G_nk
         km1=max(k-1,0)
         km2=max(k-2,0)
         km3=max(k-3,0)
         kp1=min(k+1,G_nk+1)
         kp2=min(k+2,G_nk+1)
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               ! 5th  term in elliptic RHS portion of eqn 58
               ! 5th  term in matvec portion of eqn 58
               M_logJzu(i,j,k)= Hderiv8(M_Jzu(i-1,j,k),M_Jzu(i  ,j,k),&
                                        M_Jzu(i+1,j,k),M_Jzu(i+2,j,k),geomh_invDX_8(j))
               ! 6th  term in elliptic RHS portion of eqn 58
               ! 6th  term in matvec portion of eqn 58
               M_logJzv(i,j,k)= Hderiv8(M_Jzv(i,j-1,k),M_Jzv(i  ,j,k),M_Jzv(i,j+1,k),M_Jzv(i,j+2,k),geomh_invDY_8)
               ! 7th and 9th  term in elliptic RHS portion of eqn 58
               ! 7th  term in matvec portion of eqn 58
               M_logJzq(i,j,k)= M_Jzq(i,j,km2) * CDt2m(1,k)&
                               +M_Jzq(i,j,km1) * CDt2m(2,k)&
                               +M_Jzq(i,j,k  ) * CDt2m(3,k)&
                               +M_Jzq(i,j,kp1) * CDt2m(4,k)
            end do
         end do
      end do
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               GVM%mc_Jx_8 (i,j,k)= Hderiv8(GVM%zmom_8(i-1,j,k),GVM%zmom_8(i  ,j,k),&
                                            GVM%zmom_8(i+1,j,k),GVM%zmom_8(i+2,j,k),geomh_invDX_8(j))
               GVM%mc_Jy_8 (i,j,k)= Hderiv8(GVM%zmom_8(i,j-1,k),GVM%zmom_8(i  ,j,k),&
                                            GVM%zmom_8(i,j+1,k),GVM%zmom_8(i,j+2,k),geomh_invDY_8)
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
      include 'H3rd_ope.inc'
      end subroutine vertical_metric3rd
