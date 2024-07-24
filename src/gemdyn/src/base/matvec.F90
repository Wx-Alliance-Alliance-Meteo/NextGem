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
      integer :: i, j, k, kn, km, kp, n
      integer :: i0,in,j0,jn
      real(kind=REAL64), parameter :: half=0.5d0, zero_k_8=0.d0
      integer(kind=8) :: a,b
!
!     ---------------------------------------------------------------
!
      call gtmg_start (72, 'MATVEC1', 29 )
      i0 = 1+pil_w ; in = l_ni-pil_e
      j0 = 1+pil_s ; jn = l_nj-pil_n

      kn=l_nk
      if (EZ_newsol) kn=1

      do k= 1, l_nk
         do j=j0, jn
            do i=i0, in
               fdg2(i,j,k)= F_vector(i,j,k)
            end do
         end do
      end do

      do j= j0, jn
         do i= i0, in
            fdg2(i,j,l_nk+1)=GVM%mc_alfas_H_8(i,j) * F_vector(i,j,l_nk) &
                            -GVM%mc_betas_H_8(i,j) * F_vector(i,j,l_nk-1)
         end do
      end do

      if ( Grd_yinyang_L) then
         call yyg_xchng_hlt (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, l_nk+1, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( fdg2(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif
      call gtmg_stop (72)
      call gtmg_start (73, 'MATVEC2', 29 )

      k=1; km=1 ; kp=2
      do j= j0, jn
       do i= i0, in
            F_prod(i,j,k)= &
               !+Dx[Dx[q]]
              +geomh_invDXM_8(j)*(geomh_invDX_8(j)*(fdg2(i+1,j,k)  -fdg2(i  ,j,k))*matbc_e_8(i,j,k)   & 
                                 -geomh_invDX_8(j)*(fdg2(i  ,j,k)  -fdg2(i-1,j,k))*matbc_w_8(i,j,k) ) &
               !-Dx[Jx*<Jz^(-1)*Dz[q]>^(xz)]
              -geomh_invDXM_8(j)*(&
                      GVM%mc_Jx_8(i  ,j,k)*matbc_e_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i+1,j,kp)-fdg2(i+1,j,k))*GVM%mc_iJz_8(i+1,j,k )  &
                                                         +(fdg2(i , j,kp)-fdg2(i  ,j,k))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i+1,j,k)-fdg2(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)  &
                                                         +(fdg2(i  ,j,k)-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                     -GVM%mc_Jx_8(i-1,j,k)*matbc_w_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) &
                                                         +(fdg2(i-1,j,kp)-fdg2(i-1,j,k ))*GVM%mc_iJz_8(i-1,j,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km) &
                                                         +(fdg2(i-1,j,k )-fdg2(i-1,j,km))*GVM%mc_iJz_8(i-1,j,km))&
                                         )&
               )&
               !+cos(theta)^(-1)*Dy[cos(theta)*Dy[q]]
              +geomh_invDYM_8(j)*(geomh_cyM_8(j  )*geomh_invDYMv_8(j  )*(fdg2(i,j+1,k)-fdg2(i,j  ,k))*matbc_n_8(i,j,k)&
                                 -geomh_cyM_8(j-1)*geomh_invDYMv_8(j-1)*(fdg2(i,j,k)  -fdg2(i,j-1,k))*matbc_s_8(i,j,k)&
               )&
               !-cos(theta)^(-1)*Dy[cos(theta)*[Jy*<Jz^(-1)*Dz[q]>^(yz)]]
              -geomh_invDYM_8(j)*(&
                        geomh_cyM_8(j)*GVM%mc_Jy_8(i,j,k)    *matbc_n_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i ,j+1,kp)-fdg2(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )  &
                                                         +(fdg2(i ,j  ,kp)-fdg2(i,j  ,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i ,j+1,k )-fdg2(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)  &
                                                         +(fdg2(i ,j  ,k)-fdg2(i,j  ,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       -geomh_cyM_8(j-1)*GVM%mc_Jy_8(i,j-1,k)*matbc_s_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i,j  ,k ) &
                                                         +(fdg2(i,j-1,kp)-fdg2(i,j-1,k ))*GVM%mc_iJz_8(i,j-1,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i,  j,km))*GVM%mc_iJz_8(i,j  ,km) &
                                                         +(fdg2(i,j-1,k )-fdg2(i,j-1,km))*GVM%mc_iJz_8(i,j-1,km))&
                                         )&
               )&
               !+Dz[gama*Jz^(-1)*Dz[q]>]
              +Ver_idz_8%m(k)*gama_8*((fdg2(i,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k  )&
                            -zero_k_8*(fdg2(i,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,k-1)&
               )&
               !-Dz[mu*gama*<q>^(z)]
              -gama_8*mu_8*half*Ver_idz_8%m(k)*(fdg2(i,j,kp)+fdg2(i,j,k))&
               !-epsi*gama*<Jz^(-1) Dz[q]>^(z)
              -epsi_8*gama_8*(   Ver_wp_8%m(k)*(fdg2(i ,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k ) &
                       +zero_k_8*Ver_wm_8%m(k)*(fdg2(i ,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,km) &
               )&
               !+epsi*gama*mu*<q>^(z)^(z)
              +epsi_8*gama_8*mu_8*half*( Ver_wp_8%m(k)*(fdg2(i ,j,kp)+fdg2(i,j,k )) &
                               +zero_k_8*Ver_wm_8%m(k)*(fdg2(i ,j,k )+fdg2(i,j,km)) &
               )&
               !<Dx[q]>^(x)*Dx(ln(Jz))
              +half*GVM%mc_Ix_8(i,j,k)*geomh_invDX_8(j)*( (fdg2(i+1,j,k)-fdg2(i  ,j,k))*matbc_e_8(i,j,k) &
                                                         +(fdg2(i,j,k)  -fdg2(i-1,j,k))*matbc_w_8(i,j,k) &
               )&
               !-<Jx*<Jz^(-1)*Dz[q]>^(xz)>^(x)*Dx(ln(Jz))
              -GVM%mc_Ix_8(i,j,k)*(&
                        half*GVM%mc_Jx_8(i,j,k)*matbc_e_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i+1,j,kp)-fdg2(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )  &
                                                         +(fdg2(i , j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i+1,j,k )-fdg2(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)  &
                                                         +(fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       +half*GVM%mc_Jx_8(i-1,j,k)*matbc_w_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k )  &
                                                         +(fdg2(i-1,j,kp)-fdg2(i-1,j,k ))*GVM%mc_iJz_8(i-1,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)  &
                                                         +(fdg2(i-1,j,k )-fdg2(i-1,j,km))*GVM%mc_iJz_8(i-1,j,km))&
                                         )&
               )&
               !+<Dy[q]>^(y)*Dy(ln(Jz))
              +half*GVM%mc_Iy_8(i,j,k)*(geomh_invDYMv_8(j  )*(fdg2(i,j+1,k)-fdg2(i,j  ,k)*matbc_n_8(i,j,k)) &
                                       +geomh_invDYMv_8(j-1)*(fdg2(i,j,k)  -fdg2(i,j-1,k)*matbc_s_8(i,j,k)))&
               !-<Jy*<Jz^(-1)*Dz[q]>^(yz)>^(y)*Dy(ln(Jz))
              -GVM%mc_Iy_8(i,j,k)*(&
                        half*GVM%mc_Jy_8(i,j,k)*matbc_n_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i ,j+1,kp)-fdg2(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )  &
                                                         +(fdg2(i ,j  ,kp)-fdg2(i,j  ,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i ,j+1,k )-fdg2(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)  &
                                                         +(fdg2(i ,j  ,k )-fdg2(i,j  ,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       +half*GVM%mc_Jy_8(i,j-1,k)*matbc_s_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i,j  ,k ) &
                                                         +(fdg2(i,j-1,kp)-fdg2(i,j-1,k ))*GVM%mc_iJz_8(i,j-1,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i,  j,km))*GVM%mc_iJz_8(i,j  ,km) &
                                                         +(fdg2(i,j-1,k )-fdg2(i,j-1,km))*GVM%mc_iJz_8(i,j-1,km))&
                                         )&
               )&
               !+gama*<Jz^(-1)*Dz[q]>^(z)*Dz(ln(Jz))
              +gama_8*GVM%mc_Iz_8(i,j,k)*(Ver_wp_8%m(k)*(fdg2(i,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k  )&
                                +zero_k_8*Ver_wm_8%m(k)*(fdg2(i,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,k-1)&
               )&
               !-mu*gama*<q>^(z)^(z)*Dz(ln(Jz))
              -mu_8*gama_8*GVM%mc_Iz_8(i,j,k)*half*(Ver_wp_8%m(k)*(fdg2(i,j,kp)+fdg2(i,j,k ))&
                                          +zero_k_8*Ver_wm_8%m(k)*(fdg2(i,j,k )+fdg2(i,j,km))&
               )&
               !-gamma*q
              -gg_8*fdg2(i,j,k)
       end do
      end do

      do k=2,kn
         do j= j0, jn
         km=k-1
         kp=k+1
         do i= i0, in
            F_prod(i,j,k)= &
               !+Dx[Dx[q]]
              +geomh_invDXM_8(j)*(geomh_invDX_8(j)*(fdg2(i+1,j,k)  -fdg2(i  ,j,k))*matbc_e_8(i,j,k)   & 
                                 -geomh_invDX_8(j)*(fdg2(i  ,j,k)  -fdg2(i-1,j,k))*matbc_w_8(i,j,k) ) &
               !-Dx[Jx*<Jz^(-1)*Dz[q]>^(xz)]
              -geomh_invDXM_8(j)*(&
                      GVM%mc_Jx_8(i  ,j,k)*matbc_e_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i+1,j,kp)-fdg2(i+1,j,k))*GVM%mc_iJz_8(i+1,j,k )  &
                                                         +(fdg2(i , j,kp)-fdg2(i  ,j,k))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i+1,j,k)-fdg2(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)  &
                                                         +(fdg2(i  ,j,k)-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                     -GVM%mc_Jx_8(i-1,j,k)*matbc_w_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) &
                                                         +(fdg2(i-1,j,kp)-fdg2(i-1,j,k ))*GVM%mc_iJz_8(i-1,j,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km) &
                                                         +(fdg2(i-1,j,k )-fdg2(i-1,j,km))*GVM%mc_iJz_8(i-1,j,km))&
                                         )&
               )&
               !+cos(theta)^(-1)*Dy[cos(theta)*Dy[q]]
              +geomh_invDYM_8(j)*(geomh_cyM_8(j  )*geomh_invDYMv_8(j  )*(fdg2(i,j+1,k)-fdg2(i,j  ,k))*matbc_n_8(i,j,k)&
                                 -geomh_cyM_8(j-1)*geomh_invDYMv_8(j-1)*(fdg2(i,j,k)  -fdg2(i,j-1,k))*matbc_s_8(i,j,k)&
               )&
               !-cos(theta)^(-1)*Dy[cos(theta)*[Jy*<Jz^(-1)*Dz[q]>^(yz)]]
              -geomh_invDYM_8(j)*(&
                        geomh_cyM_8(j)*GVM%mc_Jy_8(i,j,k)    *matbc_n_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i ,j+1,kp)-fdg2(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )  &
                                                         +(fdg2(i ,j  ,kp)-fdg2(i,j  ,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i ,j+1,k )-fdg2(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)  &
                                                         +(fdg2(i ,j  ,k)-fdg2(i,j  ,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       -geomh_cyM_8(j-1)*GVM%mc_Jy_8(i,j-1,k)*matbc_s_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i,j  ,k ) &
                                                         +(fdg2(i,j-1,kp)-fdg2(i,j-1,k ))*GVM%mc_iJz_8(i,j-1,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i,  j,km))*GVM%mc_iJz_8(i,j  ,km) &
                                                         +(fdg2(i,j-1,k )-fdg2(i,j-1,km))*GVM%mc_iJz_8(i,j-1,km))&
                                         )&
               )&
               !+Dz[gama*Jz^(-1)*Dz[q]>]
              +Ver_idz_8%m(k)*gama_8*((fdg2(i,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k  )&
                                     -(fdg2(i,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,k-1)&
               )&
               !-Dz[mu*gama*<q>^(z)]
              -gama_8*mu_8*half*Ver_idz_8%m(k)*(fdg2(i,j,kp)-fdg2(i,j,km))&
               !-epsi*gama*<Jz^(-1) Dz[q]>^(z)
              -epsi_8*gama_8*(   Ver_wp_8%m(k)*(fdg2(i ,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k ) &
                                +Ver_wm_8%m(k)*(fdg2(i ,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,km) &
               )&
               !+epsi*gama*mu*<q>^(z)^(z)
              +epsi_8*gama_8*mu_8*half*( Ver_wp_8%m(k)*(fdg2(i ,j,kp)+fdg2(i,j,k )) &
                                        +Ver_wm_8%m(k)*(fdg2(i ,j,k )+fdg2(i,j,km)) &
               )&
               !<Dx[q]>^(x)*Dx(ln(Jz))
              +half*GVM%mc_Ix_8(i,j,k)*geomh_invDX_8(j)*( (fdg2(i+1,j,k)-fdg2(i  ,j,k))*matbc_e_8(i,j,k) &
                                                         +(fdg2(i,j,k)  -fdg2(i-1,j,k))*matbc_w_8(i,j,k) &
               )&
               !-<Jx*<Jz^(-1)*Dz[q]>^(xz)>^(x)*Dx(ln(Jz))
              -GVM%mc_Ix_8(i,j,k)*(&
                        half*GVM%mc_Jx_8(i,j,k)*matbc_e_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i+1,j,kp)-fdg2(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )  &
                                                         +(fdg2(i , j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i+1,j,k )-fdg2(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)  &
                                                         +(fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       +half*GVM%mc_Jx_8(i-1,j,k)*matbc_w_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k )  &
                                                         +(fdg2(i-1,j,kp)-fdg2(i-1,j,k ))*GVM%mc_iJz_8(i-1,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km)  &
                                                         +(fdg2(i-1,j,k )-fdg2(i-1,j,km))*GVM%mc_iJz_8(i-1,j,km))&
                                         )&
               )&
               !+<Dy[q]>^(y)*Dy(ln(Jz))
              +half*GVM%mc_Iy_8(i,j,k)*(geomh_invDYMv_8(j  )*(fdg2(i,j+1,k)-fdg2(i,j  ,k)*matbc_n_8(i,j,k)) &
                                       +geomh_invDYMv_8(j-1)*(fdg2(i,j,k)  -fdg2(i,j-1,k)*matbc_s_8(i,j,k)))&
               !-<Jy*<Jz^(-1)*Dz[q]>^(yz)>^(y)*Dy(ln(Jz))
              -GVM%mc_Iy_8(i,j,k)*(&
                        half*GVM%mc_Jy_8(i,j,k)*matbc_n_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i ,j+1,kp)-fdg2(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )  &
                                                         +(fdg2(i ,j  ,kp)-fdg2(i,j  ,k ))*GVM%mc_iJz_8(i  ,j,k )) &
                                     +Ver_wm_8%m(k)*half*((fdg2(i ,j+1,k )-fdg2(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)  &
                                                         +(fdg2(i ,j  ,k )-fdg2(i,j  ,km))*GVM%mc_iJz_8(i  ,j,km)) &
                                         )&
                       +half*GVM%mc_Jy_8(i,j-1,k)*matbc_s_8(i,j,k) &
                                    *(Ver_wp_8%m(k)*half*((fdg2(i  ,j,kp)-fdg2(i  ,j,k ))*GVM%mc_iJz_8(i,j  ,k ) &
                                                         +(fdg2(i,j-1,kp)-fdg2(i,j-1,k ))*GVM%mc_iJz_8(i,j-1,k ))&
                                     +Ver_wm_8%m(k)*half*((fdg2(i  ,j,k )-fdg2(i,  j,km))*GVM%mc_iJz_8(i,j  ,km) &
                                                         +(fdg2(i,j-1,k )-fdg2(i,j-1,km))*GVM%mc_iJz_8(i,j-1,km))&
                                         )&
               )&
               !+gama*<Jz^(-1)*Dz[q]>^(z)*Dz(ln(Jz))
              +gama_8*GVM%mc_Iz_8(i,j,k)*(Ver_wp_8%m(k)*(fdg2(i,j,kp)-fdg2(i,j,k ))*GVM%mc_iJz_8(i,j,k  )&
                                         +Ver_wm_8%m(k)*(fdg2(i,j,k )-fdg2(i,j,km))*GVM%mc_iJz_8(i,j,k-1)&
               )&
               !-mu*gama*<q>^(z)^(z)*Dz(ln(Jz))
              -mu_8*gama_8*GVM%mc_Iz_8(i,j,k)*half*(Ver_wp_8%m(k)*(fdg2(i,j,kp)+fdg2(i,j,k ))&
                                                   +Ver_wm_8%m(k)*(fdg2(i,j,k )+fdg2(i,j,km))&
               )&
               !-gamma*q
              -gg_8*fdg2(i,j,k)
            end do
         end do
       end do

      call gtmg_stop (73)

      if (EZ_newsol) call ez_matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )

!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine matvec
