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

      subroutine matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
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

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, k0, k0t, km, kp, n, ii
      integer :: km1,km2,km3,kp1,kp2,kp3
      real(kind=REAL64) :: r1(l_ni), r2(l_ni), r3(l_ni), r4(l_ni)
      real(kind=REAL64) :: w1(l_ni),w2(l_ni),w3(l_ni),w4(l_ni),&
          w5(l_ni),w6(l_ni),w7(l_ni),w8(l_ni),w9(l_ni),w10(l_ni)
      real(kind=REAL64) :: s1(l_ni),s2(l_ni),s3(l_ni),s4(l_ni),&
                           s5(l_ni),s6(l_ni),s7(l_ni),s8(l_ni)
      real(kind=REAL64) :: dqx(-2:3), dqy(-2:3), qbz, u, v
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      if (Schm_POSO == 5) then
         call matvec5th ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                          F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
         return
      else if (Schm_POSO == 3) then
         call matvec3rd ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                          F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
         return
      endif

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
         r2(i)= (GVM%zmom_8(i,j,1)-ver_z_8%m(0))/(GVM%zmom_8(i,j,2)-GVM%zmom_8(i,j,1))
         ext_q(i,j,l_nk+1)=GVM%mc_alfas_H_8(i,j) * F_vector(i,j,l_nk) &
                          -GVM%mc_betas_H_8(i,j) * F_vector(i,j,l_nk-1)
        ! ext_q(i,j,l_nk+1)=  (1+r1)*F_vector(i,j,l_nk  ) &
        !                   -    r1 *F_vector(i,j,l_nk-1)
         ext_q(i,j,0)= (1+r2(i))*F_vector(i,j,1) &
                         -r2(i) *F_vector(i,j,2)
      end do
      end do

      if ( Grd_yinyang_L) then
         call yyg_xchng_8 (ext_q, YYG_HALO_q2q,l_minx,l_maxx,l_miny,l_maxy, &
                    l_ni,l_nj, l_nk+2, .false., 'CUBIC', .true.)
      else
         call HLT_split (0, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 ( ext_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
      k=1; km=1 ; kp=2
      do j= ds_j0, ds_jn
!DIR$ SIMD
         do i= ds_i0, ds_in

            r1(i)= GVM%mc_iJz_8(i,j,k  )*(ext_q(i,j,kp)-ext_q(i,j,k  )) !dqdz (k)
            r3(i)= half*(ext_q(i,j,k)+ext_q(i,j,kp))                    !qbarz(k)
            
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
            + gama_8*epsi_8*mu_8*half*( Ver_wp_8%m(k)*(ext_q(i ,j,kp)+ext_q(i,j,k )))
!!$         + gama_8*epsi_8*mu_8*half*( Ver_wm_8%m(k)*ext_q(i,j,km) + Ver_wp_8%m(k)*ext_q(i,j,kp) + ext_q(i,j,k) ) &
               
         end do 
      end do
      
      do k= 2, l_nk
         do j= ds_j0, ds_jn
!DIR$ SIMD
            do i= ds_i0, ds_in

               r1(i)= GVM%mc_iJz_8(i,j,k  )*(ext_q(i,j,k+1)-ext_q(i,j,k  )) !dqdz(k  )
               r2(i)= GVM%mc_iJz_8(i,j,k-1)*(ext_q(i,j,k  )-ext_q(i,j,k-1)) !dqdz(k-1)
               r3(i)= half*(ext_q(i,j,k  )+ext_q(i,j,k+1)) !qbarz(k  )
               r4(i)= half*(ext_q(i,j,k-1)+ext_q(i,j,k  )) !qbarz(k-1)
               
               s1(i)= GVM%mc_iJz_8(i+1,j,k  )*(ext_q(i+1,j,k+1)-ext_q(i+1,j,k  ))*m_east (i) !dqdz(i+1,k  )
               s2(i)= GVM%mc_iJz_8(i+1,j,k-1)*(ext_q(i+1,j,k  )-ext_q(i+1,j,k-1))*m_east (i) !dqdz(i+1,k-1)
               s3(i)= GVM%mc_iJz_8(i-1,j,k  )*(ext_q(i-1,j,k+1)-ext_q(i-1,j,k  ))*m_west (i) !dqdz(i-1,k  )
               s4(i)= GVM%mc_iJz_8(i-1,j,k-1)*(ext_q(i-1,j,k  )-ext_q(i-1,j,k-1))*m_west (i) !dqdz(i-1,k-1)
               s5(i)= GVM%mc_iJz_8(i,j+1,k  )*(ext_q(i,j+1,k+1)-ext_q(i,j+1,k  ))*m_north(j) !dqdz(j+1,k  )
               s6(i)= GVM%mc_iJz_8(i,j+1,k-1)*(ext_q(i,j+1,k  )-ext_q(i,j+1,k-1))*m_north(j) !dqdz(j+1,k-1)
               s7(i)= GVM%mc_iJz_8(i,j-1,k  )*(ext_q(i,j-1,k+1)-ext_q(i,j-1,k  ))*m_south(j) !dqdz(j-1,k  )
               s8(i)= GVM%mc_iJz_8(i,j-1,k-1)*(ext_q(i,j-1,k  )-ext_q(i,j-1,k-1))*m_south(j) !dqdz(j-1,k-1)
               
               w1(i)= (ext_q(i+1,j,k)-ext_q(i  ,j,k))*geomh_invDXM_8(j)*m_east (i) !dqdx(i  )
               w2(i)= (ext_q(i  ,j,k)-ext_q(i-1,j,k))*geomh_invDXM_8(j)*m_west (i) !dqdx(i-1)
               w3(i)= (ext_q(i,j+1,k)-ext_q(i  ,j,k))*geomh_invDYM_8(j)*m_north(j) !dqdy(j  )
               w4(i)= (ext_q(i,j  ,k)-ext_q(i,j-1,k))*geomh_invDYM_8(j-1)*m_south(j) !dqdy(j-1)
               w5(i)= half*( (r1(i)+s1(i))*Ver_wp_8%m(k) + (r2(i)+s2(i))*Ver_wm_8%m(k) )
               w6(i)= half*( (s3(i)+r1(i))*Ver_wp_8%m(k) + (s4(i)+r2(i))*Ver_wm_8%m(k) )

               w7(i)= half*( (r1(i)+s5(i))*Ver_wp_8%m(k) + (r2(i)+s6(i))*Ver_wm_8%m(k) )
               w8(i)= half*( (r1(i)+s7(i))*Ver_wp_8%m(k) + (r2(i)+s8(i))*Ver_wm_8%m(k) )
               
               w9(i) =r1(i)*Ver_wp_8%m(k  ) + r2(i)*Ver_wm_8%m(k  )
               w10(i)=r2(i)*Ver_wp_8%m(k-1) + r3(i)*Ver_wm_8%m(k-1)
               
               F_prod(i,j,k)= -gg_8*ext_q(i,j,k) &
               +geomh_invDX_8 (j)*( (w1(i)-w2(i)) - (GVM%mc_Jx_8(i,j,k)*w5(i)-GVM%mc_Jx_8(i-1,j,k)*w6(i))) &
               +geomh_invDYM_8(j)*( (geomh_cyv_8(j)*w3(i)-geomh_cyv_8(j-1)*w4(i)) &
               -GVM%mc_Jy_8(i,j,k)*geomh_cyv_8(j)*w7(i) + GVM%mc_Jy_8(i,j-1,k)*geomh_cyv_8(j-1)*w8(i)) &
               +gama_8*Ver_idz_8%m(k)*((r1(i)-mu_8*r3(i)) - (r2(i)-mu_8*r4(i))) &
               
               -gama_8*epsi_8*w9(i) &
!               +gama_8*epsi_8*mu_8*ext_q(i,j,k) &
               +gama_8*epsi_8* mu_8*half* ( Ver_wm_8%m(k)*ext_q(i,j,k-1) + Ver_wp_8%m(k)*ext_q(i,j,k+1) + ext_q(i,j,k) ) &
               
               +(half*(w1(i)+w2(i)) - GVM%mc_Jx_8(i,j,k)*w9(i))*GVM%mc_Ix_8(i,j,k) &
               +(half*(w3(i)+w4(i)) - GVM%mc_Jy_8(i,j,k)*w9(i))*GVM%mc_Iy_8(i,j,k) &
               +gama_8* (w9(i)-mu_8*ext_q(i,j,k)) * GVM%mc_Iz_8(i,j,k)

            end do
         end do
      end do
      
      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine matvec
