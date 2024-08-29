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
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )
!!$!      if (F_minx/=ldnh_minx) then
!!$         i0 = 1+pil_w ; in = l_ni-pil_e
!!$         j0 = 1+pil_s ; jn = l_nj-pil_n
!!$!      else
!!$!         i0 = 1 ; in = l_ni
!!$!         j0 = 1 ; jn = l_nj
!!$!      endif

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
         call yyg_xchng_hlt (ext_q, l_minx,l_maxx,l_miny,l_maxy, &
                    l_ni,l_nj, l_nk+2, .false., 'CUBIC', .true.)
      else
         call HLT_split (0, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( ext_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
!!$      i0 = 1+pil_w ; in = l_ni-pil_e
!!$      j0 = 1+pil_s ; jn = l_nj-pil_n
      
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
            w4(i)= (ext_q(i,j  ,k)-ext_q(i,j-1,k))*geomh_invDYM_8(j)*m_south(j) !dqdy(j-1)
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

! CODE SUGGESTION FOR 5TH ORDER PROJECT

      !QD mimics linear at the surface 
      !QW also defined 1:l_nk
      !loop can be extended to G_nk
      !also linear at the lid 

      do k=1,G_nk !extend to include last level
         km1=max(k-1,1)
         km2=max(k-2,1)
         kp1=min(k+1,G_nk+1)
         kp2=min(k+2,G_nk+1)
         kp3=min(k+3,G_nk+1)
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               do ii=-2,3
               dqx(ii) = ext_q(i+ii,j,km2) * QDm2t(1,k)&
                        +ext_q(i+ii,j,km1) * QDm2t(2,k)&
                        +ext_q(i+ii,j,k  ) * QDm2t(3,k)&
                        +ext_q(i+ii,j,k+1) * QDm2t(4,k)&
                        +ext_q(i+ii,j,kp2) * QDm2t(5,k)&
                        +ext_q(i+ii,j,kp3) * QDm2t(6,k)
               dqy(ii) = ext_q(i,j+ii,km2) * QDm2t(1,k)&
                        +ext_q(i,j+ii,km1) * QDm2t(2,k)&
                        +ext_q(i,j+ii,k  ) * QDm2t(3,k)&
                        +ext_q(i,j+ii,k+1) * QDm2t(4,k)&
                        +ext_q(i,j+ii,kp2) * QDm2t(5,k)&
                        +ext_q(i,j+ii,kp3) * QDm2t(6,k)
               end do
               !do not use Hstag because the spacing for q is not constant
               qbz  = ext_q(i,j,km2) * QWm2t(1,k)&
                    +ext_q(i,j,km1) * QWm2t(2,k)&
                    +ext_q(i,j,k  ) * QWm2t(3,k)&
                    +ext_q(i,j,k+1) * QWm2t(4,k)&
                    +ext_q(i,j,kp2) * QWm2t(5,k)&
                    +ext_q(i,j,kp3) * QWm2t(6,k)
                    
               dqdzu(i,j,k)= Hstag8(dqx(-2),dqx(-1),dqx(0),dqx(1),dqx(2),dqx (3)) !are these values in thermo, u-grid
               dqdzv(i,j,k)= Hstag8(dqy(-2),dqy(-1),dqy(0),dqy(1),dqy(2),dqy (3))
               Qw(i,j,k)= M_iJzq(i,j,k)*dqx(0) - mu_8*qbz !Dq for vertical der
            end do
         end do
      end do
      if (Grd_yinyang_L) then
         call yyg_xchng_8 (dqdzu, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk, .false., 'CUBIC', .true.)
         call yyg_xchng_8 (dqdzv, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, 2*G_nk, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 (dqdzu(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      endif

      !rhs & matvec should not address k=1&nk

      do k=2,G_nk
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               u = dqdzu(i,j,km3) * QWt2m(1,k)&
                  +dqdzu(i,j,km2) * QWt2m(2,k)&
                  +dqdzu(i,j,k-1) * QWt2m(3,k)&
                  +dqdzu(i,j,k  ) * QWt2m(4,k)&
                  +dqdzu(i,j,kp1) * QWt2m(5,k)&
                  +dqdzu(i,j,kp2) * QWt2m(6,k)

               v = dqdzv(i,j,km3) * QWt2m(1,k)&
                  +dqdzv(i,j,km2) * QWt2m(2,k)&
                  +dqdzv(i,j,k-1) * QWt2m(3,k)&
                  +dqdzv(i,j,k  ) * QWt2m(4,k)&
                  +dqdzv(i,j,kp1) * QWt2m(5,k)&
                  +dqdzv(i,j,kp2) * QWt2m(6,k)

               Qu(i,j,k)= Hderiv(ext_q(i-2,j,k),ext_q(i-1,j,k),ext_q(i  ,j,k),&
                                 ext_q(i+1,j,k),ext_q(i+2,j,k),ext_q(i+3,j,k),geomh_invDX_8(j))&
                          - M_Jxozu(i,j,k) * u !can apply to Hder to Qu/Qv; Horizontal can do Hder/Hstag

               Qv(i,j,k)= Hderiv(ext_q(i,j-2,k),ext_q(i,j-1,k),ext_q(i  ,j,k),&
                                 ext_q(i,j+1,k),ext_q(i,j+2,k),ext_q(i,j+3,k),geomh_invDY_8)&
                          - M_Jyozv(i,j,k) * v
            end do
         end do
      end do
      if (Grd_yinyang_L) then
         call yyg_xchng_8 (Qu, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk, .false., 'CUBIC', .true.)
         call yyg_xchng_8 (Qv, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, G_nk, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, 2*G_nk, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 (Qu(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,HLT_np,-1 )
      endif


!-----NEW MATVEC OPERATION-----
      do k=2,G_nk-1

         km1=max(k-1,1)
         km2=max(k-2,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         kp3=min(k+3,G_nk)

         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in

              !---x derivative of Qu---
              dxQu = Hderiv8(Qu(i-2,j,k), Qu(i-1,j,k), Qu(i  ,j,k), &
                             Qu(i+1,j,k), Qu(i+2,j,k), Qu(i+3,j,k), geomh_invDXM_8(j))

              !---y derivative of Qv---
              dyQv = Hderiv8( Qv(i,j-2,k)*geomh_cyM_8(j-2), &
                              Qv(i,j-1,k)*geomh_cyM_8(j-1), &
                              Qv(i,j  ,k)*geomh_cyM_8(j  ), &
                              Qv(i,j+1,k)*geomh_cyM_8(j+1), &
                              Qv(i,j+2,k)*geomh_cyM_8(j+2), &
                              Qv(i,j+3,k)*geomh_cyM_8(j+3), &
                              geomh_invDYM_8(j))

              !---x horizontal staggering of Qu---
              barxQu = Hstag8(Qu(i-2,j,k), Qu(i-1,j,k), Qu(i  ,j,k), &
                              Qu(i+1,j,k), Qu(i+2,j,k), Qu(i+3,j,k))

              !---y horizontal staggering of Qv---
              baryQv = Hstag8(Qv(i,j-2,k), Qv(i,j-1,k), Qv(i,j  ,k), &
                              Qv(i,j+1,k), Qv(i,j+2,k), Qv(i,j+3,k))

              !---vertical staggering of Qw---
              !thermo -> mom lvl
              barzQw = Qw(i,j,km3) * QWt2m(1,k)&
                     + Qw(i,j,km2) * QWt2m(2,k)&
                     + Qw(i,j,k-1) * QWt2m(3,k)&
                     + Qw(i,j,k  ) * QWt2m(4,k)&
                     + Qw(i,j,kp1) * QWt2m(5,k)&
                     + Qw(i,j,kp2) * QWt2m(6,k)

              !---vertical derivative of Qv---
              !thermo -> mom lvl
              dzQw =   Qw(i,j,km3) * QDt2m(1,k)&
                     + Qw(i,j,km2) * QDt2m(2,k)&
                     + Qw(i,j,k-1) * QDt2m(3,k)&
                     + Qw(i,j,k  ) * QDt2m(4,k)&
                     + Qw(i,j,kp1) * QDt2m(5,k)&
                     + Qw(i,j,kp2) * QDt2m(6,k)

              !F_prod(i,j,k) = -gg_8*ext_q(i,j,k)                              &
              !              + dxQu + dyQv + gama_8*dzQw                       &
              !              + barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k) &
              !              + gama_8*barzQw*(M_logJzq(i,j,k) - epsi_8)

            end do
         end do
      end do

!-----------------------------
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
               w4(i)= (ext_q(i,j  ,k)-ext_q(i,j-1,k))*geomh_invDYM_8(j)*m_south(j) !dqdy(j-1)
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
      include 'H5th_ope.inc'
      end subroutine matvec
