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

!** matvec_init_GY compute Sol_stencils for Matrix-vector product

      subroutine matvec_init_GY ()
      use cstv
      use geomh
      use gem_options
      use dynkernel_options
      use HORgrid_options
      use glb_ld
      use ldnh
      use opr
      use sol_mem
      use metric
      use ver
      use mem_tstp
      use lam_options
      use dyn_fisl_options
      use, intrinsic :: iso_fortran_env
      implicit none

      integer :: i,j,k, id, i0,j0, k0, km
      integer :: sol_pil_w_ext, sol_pil_e_ext, sol_pil_s_ext, sol_pil_n_ext
      real(kind=REAL64)  :: aa1, aa2, bb1, bb2, cc1
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      i0= 1+pil_w+min(pil_w,1)
      j0= 1+pil_s+min(pil_s,1)
      k0= 1+Lam_gbpil_T
      
      sol_pil_s_ext= pil_s-1
      sol_pil_n_ext= pil_n
      sol_pil_w_ext= pil_w-1
      sol_pil_e_ext= pil_e
      
      A1=0. ; A2=0. ; B1=0. ; B2=0. ; C1=0.

      k=k0
      do j=1+pil_s, l_nj-pil_n
      do i=1+pil_w, l_ni-pil_e
       !  C1(i,j,k,1)=-gama_8*(GVM%mc_iJz_8(i,j,k ) &
       !  + mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k)) - gg_8
         C1(i,j,k,5)= gama_8*(GVM%mc_iJz_8(i,j,k ) &
         - mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
      end do
      end do

      do k = k0+1,l_nk
         do j=1+pil_s, l_nj-pil_n
         do i=1+pil_w, l_ni-pil_e
         !   C1(i,j,k,1)=-gama_8*(GVM%mc_iJz_8(i,j,k ) + GVM%mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
        !    +(GVM%mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(GVM%mc_iJz_8(i,j,k-1) -mu_8*half) &
        !    -Ver_wp_8%m(k)*(GVM%mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
            C1(i,j,k,4)= gama_8*(GVM%mc_iJz_8(i,j,k-1) +  mu_8*half)*(Ver_idz_8%m(k) &
            - (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wm_8%m(k))
            C1(i,j,k,5)= gama_8*(GVM%mc_iJz_8(i,j,k ) -mu_8*half)*(Ver_idz_8%m(k) + &
            (GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k))
         end do
         end do
      end do

      do k = k0,l_nk
         km=max(k-1,1)
         do j=1+sol_pil_s_ext, l_nj-sol_pil_n
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
           !    A1(i,j,k,1)= -geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
           !                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
           !    A1(i,j,k,3)=  geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
           !                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km))
               A1(i,j,k,4)= half*GVM%mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
               A1(i,j,k,5)=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
               A1(i,j,k,8)= half*GVM%mc_Jx_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km)
               A1(i,j,k,9)=-half*GVM%mc_Jx_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k)
         end do
         end do

         do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e
           !    B1(i,j,k,1 )= -geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
           !                  (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
               B1(i,j,k,4 )= half*GVM%mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
               B1(i,j,k,5 )=-half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
               B1(i,j,k,11)= geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*  &
                             (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km))
               B1(i,j,k,14)= half*GVM%mc_Jy_8(i,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j+1,km)
               B1(i,j,k,15)= -half*GVM%mc_Jy_8(i,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j+1,k)
         end do
         end do

         do j=1+sol_pil_s_ext, l_nj-sol_pil_n
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
           ! A2(i,j,k,1)=  geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
           !               (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
           ! A2(i,j,k,2)= -geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
           !               (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km))
            A2(i,j,k,4)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            A2(i,j,k,5)= -half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            A2(i,j,k,6)= half*GVM%mc_Jx_8(i-1,j,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,km)
            A2(i,j,k,7)= -half*GVM%mc_Jx_8(i-1,j,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k)
         end do
         end do

         do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
         do i=1+sol_pil_w_ext, l_ni-sol_pil_e

           ! B2(i,j,k,1 )= geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
           !               (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km))
            B2(i,j,k,4 )= half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,km)
            B2(i,j,k,5 )=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k)
            B2(i,j,k,10)=-geomh_invDYMv_8(j-1) + half*GVM%mc_Jy_8(i,j-1,k)* &
                          (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k ) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km))
            B2(i,j,k,12)= half*GVM%mc_Jy_8(i,j-1,k)*Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j-1,km)
            B2(i,j,k,13)=-half*GVM%mc_Jy_8(i,j-1,k)*Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j-1,k)
         enddo
         enddo
      end do

      do id=4,15
         do k =k0, l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e
                  Sol_stencilh_8 (i,j,k,id) =Cstv_hco0_8* ( (A1 (i,j,k,id)-A2 (i,j,k,id))*geomh_invDXM_8(j) &
                                                 + half * ( GVM%mc_Ix_8(i,j,k)*(A1(i,j,k,id)+A2(i,j,k,id)))    &
                                                 + (B1 (i,j,k,id)*geomh_cyM_8(j)-B2 (i,j,k,id)*                &
                                                    geomh_cyM_8(j-1))*geomh_invDYM_8(j)                        &
                                                 + half * (GVM%mc_Iy_8(i,j,k)*(B1(i,j,k,id)+B2(i,j,k,id)) )    &
                                                 + C1(i,j,k,id) )
               enddo
            end do
         end do
      end do
!!$
!!$! NEW
!!$      k= 1 ; km= 1
!!$      do j=1+pil_s, l_nj-pil_n
!!$      do i=1+pil_w, l_ni-pil_e
!!$         aa1=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,1))
!!$         aa2= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,1))
!!$         bb1=-geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,1))
!!$         bb2= geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,1))
!!$         cc1= -gama_8*(GVM%mc_iJz_8(i,j,k )+ mu_8*half)*(Ver_idz_8%m(k)+(GVM%mc_Iz_8(i,j,k)-epsi_8)*Ver_wp_8%m(k)) - gg_8
!!$         Sol_stencilh_8 (i,j,k,1)=Cstv_hco0_8* ( (aa1-aa2)*geomh_invDXM_8(j) &
!!$               + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
!!$               + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1 )
!!$         aa2=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,1))
!!$         Sol_stencilh_8 (i,j,k,2)=Cstv_hco0_8* ( (-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa2)))
!!$         aa1= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
!!$              (Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,km))
!!$         Sol_stencilh_8 (i,j,k,3)=Cstv_hco0_8* ( (aa1)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1)))
!!$      end do
!!$      end do
!!$      
!!$      do k= 2, l_nk
!!$         do j=1+pil_s, l_nj-pil_n
!!$         do i=1+pil_w, l_ni-pil_e
!!$            aa1=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
!!$                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,k-1))
!!$            aa2= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
!!$                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,k-1))
!!$            bb1=-geomh_invDYMv_8(j) + half*GVM%mc_Jy_8(i,j,k)*     &
!!$                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,k-1))
!!$            bb2= geomh_invDYMv_8(j-1)  + half*GVM%mc_Jy_8(i,j-1,k)* &
!!$                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i,j,k-1))
!!$            cc1=-gama_8*(GVM%mc_iJz_8(i,j,k ) + GVM%mc_iJz_8(i,j,k-1) )*Ver_idz_8%m(k) &
!!$                +(GVM%mc_Iz_8(i,j,k)-epsi_8)*gama_8*( Ver_wm_8%m(k)*(GVM%mc_iJz_8(i,j,k-1) -mu_8*half) &
!!$                -Ver_wp_8%m(k)*(GVM%mc_iJz_8(i,j,k )  +mu_8*half) ) - gg_8
!!$            Sol_stencilh_8 (i,j,k,1)=Cstv_hco0_8* ( (aa1-aa2)*geomh_invDXM_8(j) &
!!$               + half*(GVM%mc_Ix_8(i,j,k)*(aa1+aa2)) + (bb1*geomh_cyM_8(j)-bb2*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
!!$               + half*(GVM%mc_Iy_8(i,j,k)*(bb1+bb2)) + cc1 )
!!$            aa2=-geomh_invDX_8(j) + half*GVM%mc_Jx_8(i-1,j,k)* &
!!$                (Ver_wp_8%m(k)*GVM%mc_iJz_8(i-1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i-1,j,k-1))
!!$            Sol_stencilh_8 (i,j,k,2)=Cstv_hco0_8* ( (-aa2)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa2)))
!!$            aa1= geomh_invDX_8(j) + half*GVM%mc_Jx_8(i,j,k)*   &
!!$                 (Ver_wp_8%m(k)*GVM%mc_iJz_8(i+1,j,k) - Ver_wm_8%m(k)*GVM%mc_iJz_8(i+1,j,k-1))
!!$            Sol_stencilh_8 (i,j,k,3)=Cstv_hco0_8* ( (aa1)*geomh_invDXM_8(j) + half*(GVM%mc_Ix_8(i,j,k)*(aa1)))
!!$         end do
!!$         end do
!!$      end do
!  
!     ---------------------------------------------------------------
!
      return
      end subroutine matvec_init_GY

