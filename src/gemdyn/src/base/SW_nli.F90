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
!
!*s/r SW_nli - compute nonlinear terms for current time

      subroutine SW_nli ( F_dt_8 )
      use HORgrid_options
      use gem_options
      use dyn_fisl_options
      use dynkernel_options
      use coriolis
      use geomh
      use gmm_geof
      use tdpack
      use gmm_phy
      use glb_ld
      use gmm_contiguous
      use gmm_vt0
      use adz_mem
      use mem_tstp
      use mem_tracers
      use cstv
      use ver
      use metric
      use step_options

      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      real(kind=REAL64) :: invT_n_8
      
      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, km, kp, n
      real, dimension(:,:,:), pointer :: logT, tots
      real(kind=REAL64) :: div, barz, barzp, u_interp, v_interp, &
               t_interp, w2, w3, w4, invT_8, invT_nh_8,invT_m_8, &
               dwdz, dudz, dvdz, dudx, dvdy !this is to check each term of nlc
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      invT_n_8 = 3.0/(2.0*F_dt_8)
      HLT_j0 = 1
      HLT_jn = l_nj
      HLT_nk = l_nk
      HLT_nj = HLT_jn - HLT_j0 + 1
      call HLT_split (1, HLT_nj*HLT_nk, HLT_np, HLT_start, HLT_end)

!*********************************************************
! Compute Nu, Nv: Nonlinear terms of U, V equations      *
! Compute Nw, Nt: Nonlinear terms of w, T equations      *
!*********************************************************

      w3= 1.d0*half/(cpd_8*Cstv_Tstr_8)
      w4= 1.d0*epsi_8/grav_8

      !print *," i, j, k,     t1: |  iJz(i+1)  |  iJz(i)    |  Jx(i+1)   |  Jx(i)  |"

      !print *, "   dudx   |   dvdy   |   dudz   |   dvdz   |   dwdz   |   nl_t1"

      do n= HLT_start, HLT_end
         k= (n-1)/HLT_nj
         j= n - k*HLT_nj + HLT_j0 - 1
         k= k+1
         
         km=max(k-1,1)
         kp=min(k+1,l_nk)

         do i= 1, l_ni

            !assuming T at momentum level 1 equal to T at thermo level 3/2

            !--------nu term---------------
            v_interp = 0.25d0*(vt0(i,j,k)+vt0(i,j-1,k)+vt0(i+1,j,k)+vt0(i+1,j-1,k))

            nlu_t1(i,j,k) = - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v_interp 

            !--------nv term---------------
            u_interp = 0.25d0*(ut0(i,j,k)+ut0(i-1,j,k)+ut0(i,j+1,k)+ut0(i-1,j+1,k))

            nlv_t1(i,j,k) = ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp 

            !---nonlinear term for depth---
            div = (ut0 (i,j,k)- ut0 (i-1,j,k))*geomh_invDXM_8(j)     &
                + (vt0 (i,j,k)*geomh_cyM_8(j)-vt0 (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) 
            nlq_t1(i,j,k) = (1.d0-Cstv_swln_8)*(qt0(i,j,k)-fis0(i,j)) * div &                                                                                         
                            -Cstv_swln_8*invT_n_8*(Cstv_h0inv_8*qt0(i,j,k)-log(Cstv_h0inv_8*(qt0(i,j,k)-fis0(i,j))+1.d0))
        end do
      end do
      
!$OMP BARRIER

      call HLT_split (1, 5*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( nlu_t1(Adz_lminx,Adz_lminy,HLT_start),&
                 Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)

!     
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_nli
