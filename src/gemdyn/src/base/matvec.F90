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

!** matvec - 3D Matrix-vector product subroutines

      subroutine matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                          F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use metric
      use omp_timing
      use sol_mem
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(out) :: F_prod

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, k0, k0t, km, kp, n
      integer :: i0,in,j0,jn
      real(kind=REAL64)  :: r1, aa1, aa2, aa3, bb1, bb2, bb3, cc1,&
                 S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15
      real(kind=REAL64), parameter :: half=0.5d0
      integer(kind=8) :: a,b
!
!     ---------------------------------------------------------------
!
      if (EZ_newsol) then
         call ez_matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                          F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
         return
      endif
      
      call gtmg_start (72, 'MATVEC1', 25 )
      i0 = 1+pil_w ; in = l_ni-pil_e
      j0 = 1+pil_s ; jn = l_nj-pil_n
      k0=1+Lam_gbpil_T
      k0t=k0 ; if (Schm_opentop_L) k0t=k0-1

!!$omp do collapse(2)
      do k= k0, l_nk
         do j=j0, jn
            do i=i0, in
               fdg2(i,j,k)= F_vector(i,j,k)
            end do
         end do
      end do
!!$omp enddo

!!$omp do
         do j= j0, jn
            do i= i0, in
            fdg2(i,j,l_nk+1)=GVM%mc_alfas_H_8(i,j) * F_vector(i,j,l_nk) &
                            -GVM%mc_betas_H_8(i,j) * F_vector(i,j,l_nk-1)
         end do
      end do
!!$omp enddo

      if (Schm_opentop_L) then
!!$omp do
         do j= j0, jn
            do i= i0, in
               fdg2(i,j,k0t) = GVM%mc_alfat_8(i,j)* F_vector(i,j,k0)
            end do
         end do
!!$omp enddo
      endif

      if ( Grd_yinyang_L) then
         call yyg_xchng_hlt (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, l_nk+1, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( fdg2(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif
      call gtmg_stop (72)
      call gtmg_start (73, 'MATVEC2', 25 )

! Special treatment at k==1 for S1,S4,S5
      if (k0==1) then
         k= 1 ; km=1 ; kp=2
         do j= j0, jn
         do i= i0, in
            aa1=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            aa2= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            bb1=-geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            bb2= geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            cc1= -gama_8*(GVM%mc_iJz_8(i,j,k )+ mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k)) - gg_8
            S1= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1
            aa1=half*GVM%mc_Jx_8(i  ,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            aa2=half*GVM%mc_Jx_8(i-1,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            bb1=half*GVM%mc_Jy_8(i  ,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            bb2=half*GVM%mc_Jy_8(i  ,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            S4= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2))
            aa1=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            aa2=-half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            bb1=-half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            bb2=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            cc1=gama_8*(GVM%mc_iJz_8(i,j,k ) - mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
            S5= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1
            F_prod(i,j,k)= S1*fdg2(i,j,k ) + S4*fdg2(i,j,km) +  S5*fdg2(i,j,kp)
         end do
         end do
         k0=2
      endif
      
      do k=k0,l_nk
         do j= j0, jn
         km=max(k-1,1)
         kp=k+1
         do i= i0, in
            aa1=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            aa2= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            bb1=-geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            bb2= geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            cc1=-gama_8*(GVM%mc_iJz_8(i,j,k ) + GVM%mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
                +(GVM%mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(GVM%mc_iJz_8(i,j,k-1) -mu_8*half) &
                -Ver_wp_8%m(k)*(GVM%mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
            S1= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1
            aa1=half*GVM%mc_Jx_8(i  ,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            aa2=half*GVM%mc_Jx_8(i-1,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            bb1=half*GVM%mc_Jy_8(i  ,j  ,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            bb2=half*GVM%mc_Jy_8(i  ,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            cc1=gama_8*(GVM%mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
                - (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
            S4= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1

            aa1=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            aa2=-half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            bb1=-half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            bb2=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            if ( .not. Grd_yinyang_L) then
               if (l_east  .and. i==l_ni-pil_e)  aa1= 0.
               if (l_north .and. j==l_nj-pil_n)  bb1= 0.
               if (l_west  .and. i==1+pil_w   )  aa2= 0.
               if (l_south .and. j==1+pil_s   )  bb2= 0.
            endif
            cc1=gama_8*(GVM%mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
               (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
            S5= (aa1-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) &
                    + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1
            
            F_prod(i,j,k)= S1*fdg2(i,j,k ) + S4*fdg2(i,j,km) + S5*fdg2(i,j,kp)
         end do
         end do
      end do
      
! All the other S*
      do k= 1+Lam_gbpil_T, l_nk
         do j= j0, jn
         km=max(k-1,1)
         kp=k+1
         do i= i0, in
            r1= half*GVM%mc_Ix_8(i,j,k)-geomh_invDXM_8(j)
            aa1= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km)
            aa2=-half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k )
            aa3=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km))
            if ((.not. Grd_yinyang_L).and.(l_west .and. i==1+pil_w)) then
               aa1= 0.
               aa2= 0.
               aa3= 0.
            end if
            S6= aa1*r1 ; S7= aa2*r1 ; S2= aa3*r1
            
            r1= geomh_invDXM_8(j) + half*GVM%mc_Ix_8(i,j,k)
            aa1= half*GVM%mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km)
            aa2=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k )
            aa3= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km))
            if (( .not. Grd_yinyang_L).and.(l_east  .and. i==l_ni-pil_e)) then
               aa1= 0.
               aa2= 0.
               aa3= 0.
            endif
            S8= aa1*r1 ; S9= aa2*r1 ; S3= aa3*r1
            
            r1 = half*GVM%mc_Iy_8(i,j,k)-geomh_cyM_8(j-1)*geomh_invDYM_8(j)
            bb1= half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km)
            bb2=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k )
            bb3=-geomh_invDYMv_8(j-1) + half*GVM%mc_Jy_8(i,j-1,k)* &
               (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km))
            if (( .not. Grd_yinyang_L).and.(l_south .and. j==1+pil_s   )) then
               bb1= 0.
               bb2= 0.
               bb3= 0.
            endif
            S10= bb3*r1 ; S12= bb1*r1 ; S13= bb2*r1
            
            r1 =  half*GVM%mc_Iy_8(i,j,k)+geomh_cyM_8(j)*geomh_invDYM_8(j)
            bb1=  half*GVM%mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km)
            bb2= -half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k )
            bb3= geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*  &
                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km))
            if (( .not. Grd_yinyang_L).and.(l_north .and. j==l_nj-pil_n)) then
               bb1= 0.
               bb2= 0.
               bb3= 0.
            endif
            S11= bb3*r1 ; S14= bb1*r1 ; S15= bb2*r1

            F_prod(i,j,k)= F_prod(i,j,k) + &
               S2 *fdg2(i-1,j,k ) + S3 *fdg2(i+1,j,k ) +  S6*fdg2(i-1,j,km) +  S7*fdg2(i-1,j,kp) &
             + S8 *fdg2(i+1,j,km) + S9 *fdg2(i+1,j,kp) + S10*fdg2(i,j-1,k ) + S11*fdg2(i,j+1,k ) &
             + S12*fdg2(i,j-1,km) + S13*fdg2(i,j-1,kp) + S14*fdg2(i,j+1,km) + S15*fdg2(i,j+1,kp) 
            end do
         end do
      end do
      
      call gtmg_stop (73)
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine matvec
