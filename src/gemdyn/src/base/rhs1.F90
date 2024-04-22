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
!*s/r rhs1 - Prepare rhs* variables

      subroutine rhs1 ( F_dt_8 )
      use dyn_fisl_options
      use tdpack
      use glb_ld
      use gmm_contiguous
      use gmm_vt1
      use gmm_vt2
      use adz_mem
      use mem_tstp
      use mem_tracers
      use ver
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      
      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, km, kp, n
      real, dimension(:,:,:), pointer :: logT1, logT2
      real(kind=REAL64) :: div, barz, barzp, u_interp, v_interp,&
               t_interp, w2, w3, w4, invT_8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
!!$omp do
      do n= 1, ubound(dynt0,1)
         dynt0(n) = dynt1(n)
      end do
!!$omp end do nowait
!!$omp do
      do n= 1, ubound(trt0,1)
         trt0(n) = trt1(n)
      end do
!!$omp end do nowait

      HLT_j0 = 1
      HLT_jn = l_nj
      HLT_nk = l_nk
      HLT_nj = HLT_jn - HLT_j0 + 1
      call HLT_split (1, HLT_nj*HLT_nk, HLT_np, HLT_start, HLT_end)

!***************************************************************
! Compute rhs of terms that will be interpolated               *
!***************************************************************

      w3= half/(cpd_8*Cstv_Tstr_8)
      w4= epsi_8/grav_8

      logT1(1:l_ni,1:l_nj,1:l_nk) => WS1(1:)
      logT2(1:l_ni,1:l_nj,1:l_nk) => WS1(l_ni*l_nj*l_nk+1:)
!!$omp do collapse(2)
      do k=1, l_nk
         do j=1, l_nj
         do i= 1, l_ni
            logT1(i,j,k)= log (tt1(i,j,k)/Cstv_Tstr_8)
            logT2(i,j,k)= log (tt2(i,j,k)/Cstv_Tstr_8)
        end do
      end do
      end do
!!$omp enddo

      do n= HLT_start, HLT_end
         k= (n-1)/HLT_nj
         j= n - k*HLT_nj + HLT_j0 - 1
         k= k+1

         do i= 1, l_ni

            rhsu_bdf_t1(i,j,k) = ut1(i,j,k)
            rhsu_bdf_t2(i,j,k) = ut2(i,j,k)

            rhsv_bdf_t1(i,j,k) = vt1(i,j,k)
            rhsv_bdf_t2(i,j,k) = vt2(i,j,k)

            rhsw_bdf_t1(i,j,k) = wt1(i,j,k)
            rhsw_bdf_t2(i,j,k) = wt2(i,j,k)

            !---rhs terms of z---
            rhsf_bdf_t1(i,j,k) = GVM%ztht_8(i,j,k) - Ver_z_8%t(k)
            rhsf_bdf_t2(i,j,k) = rhsf_bdf_t1(i,j,k)

            !---rhs terms of t---
            rhst_bdf_t1(i,j,k) =  logT1(i,j,k) - w3*(qt1(i,j,k+1)+qt1(i,j,k))  
            rhst_bdf_t2(i,j,k) =  logT2(i,j,k) - w3*(qt2(i,j,k+1)+qt2(i,j,k))  

            !---rhs terms of q---
            rhsc_bdf_t1(i,j,k) = w4*qt1(i,j,k) 
            rhsc_bdf_t2(i,j,k) = w4*qt2(i,j,k) 

        end do
      end do
      
!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhst_bdf_t1(Adz_lminx,Adz_lminy,HLT_start),&
                 Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)
      call gem_xch_halo ( rhst_bdf_t2(Adz_lminx,Adz_lminy,HLT_start),&
                 Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)

!     ---------------------------------------------------------------
!
      return
      end subroutine rhs1
