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

!** matvec3D_init compute Sol_stencils for Matrix-vector product subroutines (P & H coordinates)
!
      subroutine SW_matvec_init ()
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

      integer j, jj, i, ii, id, k
      real(kind=REAL64)  :: di_8
      real(kind=REAL64)  :: xxx, yyy
      real(kind=REAL64), parameter :: one=1.d0, zero=0.d0, half=0.5d0
      integer, parameter :: IDX_POINT=1, IDX_WEST=2, IDX_EAST=3, IDX_NORTH=4, IDX_SOUTH=5, IDX_TOP=6, IDX_BOTTOM=7

      integer  km, kp,k0,k0t
      integer sol_pil_w_ext, sol_pil_e_ext, sol_pil_s_ext, sol_pil_n_ext
!
!     ---------------------------------------------------------------
!
      k0=1+Lam_gbpil_T
      k0t=k0
      if (Schm_opentop_L) k0t=k0-1
      
      sol_pil_s_ext=sol_pil_s-1
      sol_pil_n_ext= sol_pil_n
      sol_pil_w_ext=sol_pil_w-1
      sol_pil_e_ext=sol_pil_e
      
      if(.not.Grd_yinyang_L) then
         if (l_west) sol_pil_w_ext=sol_pil_w
         if (l_east) sol_pil_e_ext=sol_pil_e+1
         if (l_south) sol_pil_s_ext=sol_pil_s
         if (l_north) sol_pil_n_ext= sol_pil_n+1
      endif
      
      k=k0
!!$omp do
         do j=1+sol_pil_s, l_nj-sol_pil_n
            do i=1+sol_pil_w, l_ni-sol_pil_e
               !---original---
               C1(i,j,k,1)=-gg_sw_8
               C1(i,j,k,5)= 0.
            end do
         end do
!!$omp enddo nowait

!!$omp do
         do k = k0+1,l_nk
            do j=1+sol_pil_s, l_nj-sol_pil_n
               do i=1+sol_pil_w, l_ni-sol_pil_e
                   !---original---
                   C1(i,j,k,1)=-gg_sw_8
                   C1(i,j,k,4)=0.
                   C1(i,j,k,5)=0.
               end do
            end do
         end do
!!$omp enddo nowait

!!$omp do
         do k = k0,l_nk
            km=max(k-1,1)
            kp=k+1
            do j=1+sol_pil_s_ext, l_nj-sol_pil_n
               do i=1+sol_pil_w_ext, l_ni-sol_pil_e_ext
                  A1(i,j,k,1)= -geomh_invDX_8(j) 
                  A2(i,j,k,1)=  geomh_invDX_8(j) 
                  A2(i,j,k,2)= -geomh_invDX_8(j) 
                  A1(i,j,k,3)=  geomh_invDX_8(j) 
                  A1(i,j,k,4)= 0.
                  A2(i,j,k,4)= 0.
                  A1(i,j,k,5)= 0.
                  A2(i,j,k,5)= 0.
                  A2(i,j,k,6)= 0.
                  A2(i,j,k,7)= 0.
                  A1(i,j,k,8)= 0.
                  A1(i,j,k,9)= 0.
               end do
            end do

            do j=1+sol_pil_s_ext, l_nj-sol_pil_n_ext
               do i=1+sol_pil_w_ext, l_ni-sol_pil_e
                  B1(i,j,k,1)= -geomh_invDYMv_8(j) 
                  B2(i,j,k,1)=  geomh_invDYMv_8(j-1) 
                  B1(i,j,k,4)= 0.
                  B2(i,j,k,4)= 0.
                  B1(i,j,k,5)= 0.
                  B2(i,j,k,5)= 0.
                  B2(i,j,k,10)=-geomh_invDYMv_8(j-1) 
                  B1(i,j,k,11)= geomh_invDYMv_8(j) 
                  B2(i,j,k,12) = 0.
                  B2(i,j,k,13)= 0.
                  B1(i,j,k,14)= 0.
                  B1(i,j,k,15)= 0.
               end do
            end do
         end do

!!$omp enddo

!!$omp do collapse(3)
         do id=1,15
            do k =k0, l_nk
               do j=1+sol_pil_s, l_nj-sol_pil_n
                  do i=1+sol_pil_w, l_ni-sol_pil_e
                     !---original stencil---
                     Sol_stencilh_8 (i,j,k,id) =Cstv_hco0_8* ( (A1 (i,j,k,id)-A2 (i,j,k,id))*geomh_invDXM_8(j) &
                                                 + (B1 (i,j,k,id)*geomh_cyM_8(j)-B2 (i,j,k,id)*                &
                                                    geomh_cyM_8(j-1))*geomh_invDYM_8(j)  + C1(i,j,k,id)        )
                  enddo
               end do
            end do
         end do
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_matvec_init

