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

!** matvec - 3D Matrix-vector product

      subroutine matvec3rd ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                             F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use ldnh
      use metric
      use omp_timing
      use sol_mem
      use mem_tstp
      use tdpack
      use ver
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(IN ) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn        ,F_nk), intent(OUT) :: F_prod

      integer :: i, j, k
      integer :: km1,km2,km3,kp1,kp2,kp3
      real(kind=REAL64) :: r1(l_ni), r3(l_ni)
      real(kind=REAL64) :: w1(l_ni),w2(l_ni),w3(l_ni),w4(l_ni),&
                           w5(l_ni),w6(l_ni),w9(l_ni)
      real(kind=REAL64) :: s1(l_ni),s2(l_ni),s3(l_ni)
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )

      do k= 1, l_nk
         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            ext_q(i,j,k)= F_vector(i,j,k)
         end do
         end do
      end do
      do j= ds_j0, ds_jn
      do i= ds_i0, ds_in
         r1(i)=  (GVM%zmom_8(i,j,l_nk+1)-GVM%zmom_8(i,j,l_nk))&
                /(GVM%zmom_8(i,j,l_nk)-GVM%zmom_8(i,j,l_nk-1))
         r3(i)= (GVM%zmom_8(i,j,1)-ver_z_8%m(0))/(GVM%zmom_8(i,j,2)-GVM%zmom_8(i,j,1))
         ext_q(i,j,l_nk+1)=GVM%mc_alfas_H_8(i,j) * F_vector(i,j,l_nk) &
                          -GVM%mc_betas_H_8(i,j) * F_vector(i,j,l_nk-1)
        ! ext_q(i,j,l_nk+1)=  (1+r1)*F_vector(i,j,l_nk  ) &
        !                   -    r1 *F_vector(i,j,l_nk-1)
         ext_q(i,j,0)= (1+r3(i))*F_vector(i,j,1) &
                         -r3(i) *F_vector(i,j,2)
      end do
      end do

      call delQ3rd (ext_q, l_minx,l_maxx,l_miny,l_maxy, Qu,Qv,Qw,Qq,0,l_nk+1)

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
      k=1
      do j= ds_j0, ds_jn
!DIR$ SIMD
         do i= ds_i0, ds_in

            r1(i)= GVM%mc_iJz_8(i,j,k  )*(ext_q(i,j,2)-ext_q(i,j,k  )) !dqdz (k)
            r3(i)= half*(ext_q(i,j,k)+ext_q(i,j,2))                    !qbarz(k)
            
            s1(i)= GVM%mc_iJz_8(i+1,j,k)*(ext_q(i+1,j,k+1)-ext_q(i+1,j,k))*m_east (i) !dqdz(i+1,k  )
            s2(i)= GVM%mc_iJz_8(i+1,j,k)*(ext_q(i+1,j,k  )-ext_q(i+1,j,k))*m_east (i) !dqdz(i+1,1)
            s3(i)= GVM%mc_iJz_8(i-1,j,k)*(ext_q(i-1,j,k+1)-ext_q(i-1,j,k))*m_west (i) !dqdz(i-1,k  )
            
            w1(i)= (ext_q(i+1,j,k)-ext_q(i  ,j,k))*geomh_invDXM_8(j)*m_east (i) !dqdx(i  )
            w2(i)= (ext_q(i  ,j,k)-ext_q(i-1,j,k))*geomh_invDXM_8(j)*m_west (i) !dqdx(i-1)
            w3(i)= (ext_q(i,j+1,k)-ext_q(i  ,j,k))*geomh_invDYM_8(j)*m_north(j) !dqdy(j  )
            w4(i)= (ext_q(i,j  ,k)-ext_q(i,j-1,k))*geomh_invDYM_8(j-1)*m_south(j) !dqdy(j-1)
            w5(i)= half*( (r1(i)+s1(i))*Ver_wp_8%m(k) + s2(i)*Ver_wm_8%m(k) )
            w6(i)= half*(s3(i)+r1(i))*Ver_wp_8%m(k)
            w9(i)= r1(i)*Ver_wp_8%m(k)
            
            F_prod(i,j,k)= -gg_8*ext_q(i,j,k) &
            + geomh_invDX_8 (j)*( (w1(i)-w2(i)) - (GVM%mc_Jx_8(i,j,k)*w5(i)-GVM%mc_Jx_8(i-1,j,k)*w6(i))) &
            + geomh_invDYM_8(j)*( (geomh_cyv_8(j)*w3(i)-geomh_cyv_8(j-1)*w4(i))) &
            + gama_8*Ver_idz_8%m(k)*(r1(i)-mu_8*r3(i)) &
            + (half*(w1(i)+w2(i)) - GVM%mc_Jx_8(i,j,k)*w9(i))*GVM%mc_Ix_8(i,j,k) &
            + (half*(w3(i)+w4(i)) - GVM%mc_Jy_8(i,j,k)*w9(i))*GVM%mc_Iy_8(i,j,k) &
            + gama_8* (w9(i)-mu_8*ext_q(i,j,k)) * GVM%mc_Iz_8(i,j,k) &
            - gama_8*epsi_8*w9(i) &
            + gama_8*epsi_8*mu_8*half*( Ver_wp_8%m(k)*(ext_q(i ,j,2)+ext_q(i,j,k )))
         end do 
      end do
      
      do k= 2, l_nk

         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)

         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
!!$               dxQu = (Qu(i,j,k)-Qu(i-1,j,k))*geomh_invDXM_8(j)
!!$               dyQv = geomh_invDYM_8(j)*(Qv(i,j,k)-Qv(i,j-1,k))
!!$               barxQu = half*(Qu(i,j,k)+Qu(i-1,j,k))
!!$               baryQv = half*(Qv(i,j,k)+Qv(i,j-1,k))
!!$               barzQw = Qw(i,j,k)*Ver_wp_8%m(k) + Qw(i,j,k-1)*Ver_wm_8%m(k)
!!$               dzQw   = Ver_idz_8%m(k)*(Qw(i,j,k)-Qw(i,j,k-1))

               dxQu = Hderiv8(Qu(i-2,j,k), Qu(i-1,j,k), &
                              Qu(i  ,j,k), Qu(i+1,j,k), geomh_invDXM_8(j))

               dyQv = Hderiv8(Qv(i,j-2,k)*geomh_cyM_8(j-2), &
                              Qv(i,j-1,k)*geomh_cyM_8(j-1), &
                              Qv(i,j  ,k)*geomh_cyM_8(j  ), &
                              Qv(i,j+1,k)*geomh_cyM_8(j+1), &
                              geomh_invDYM_8(j) )

               barxQu = Hstag8(Qu(i-2,j,k), Qu(i-1,j,k),&
                               Qu(i  ,j,k), Qu(i+1,j,k) )

               baryQv = Hstag8(Qv(i,j-2,k), Qv(i,j-1,k),&
                               Qv(i,j  ,k), Qv(i,j+1,k) )

               barzQw =  Qw(i,j,km2) * CWt2m(1,k) &
                      +  Qw(i,j,km1) * CWt2m(2,k) &
                      +  Qw(i,j,k  ) * CWt2m(3,k) &
                      +  Qw(i,j,kp1) * CWt2m(4,k)

               !--- remains second order for now ---
               dzQw= Ver_idz_8%m(k)*(Qw(i,j,k)-Qw(i,j,k-1))
!!$               dzQw =  Qw(i,j,km2) * CDt2m(1,k) &
!!$                    +  Qw(i,j,km1) * CDt2m(2,k) &
!!$                    +  Qw(i,j,k  ) * CDt2m(3,k) &
!!$                    +  Qw(i,j,kp1) * CDt2m(4,k)
                    
               F_prod(i,j,k)= -gg_8*ext_q(i,j,k) + dxQu + dyQv + gama_8*dzQw &
                              +gama_8*barzQw*(M_logJzq(i,j,k)-epsi_8) &
                              +barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k)
            end do
         end do
      end do
      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine matvec3rd
