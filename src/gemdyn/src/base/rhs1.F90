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
      integer :: i, j, k, km, kp, n, km2, km1, kp1, kp2, kp3
      real, dimension(:,:,:), pointer :: logT1, logT2
      real(kind=REAL64) :: w3, w4, qbzt1, qbzt2
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

      w3= one/(cpd_8*Cstv_Tstr_8)
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

         km1=max(k-1,1)
         km2=max(k-2,1)
         kp1=min(k+1,l_nk+1)
         kp2=min(k+2,l_nk+1)
         kp3=min(k+3,l_nk+1)

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
            qbzt1  = qt1(i,j,km2) * QWm2t(1,k)&
                   + qt1(i,j,km1) * QWm2t(2,k)&
                   + qt1(i,j,k  ) * QWm2t(3,k)&
                   + qt1(i,j,kp1) * QWm2t(4,k)&
                   + qt1(i,j,kp2) * QWm2t(5,k)&
                   + qt1(i,j,kp3) * QWm2t(6,k)
            qbzt2  = qt2(i,j,km2) * QWm2t(1,k)&
                   + qt2(i,j,km1) * QWm2t(2,k)&
                   + qt2(i,j,k  ) * QWm2t(3,k)&
                   + qt2(i,j,kp1) * QWm2t(4,k)&
                   + qt2(i,j,kp2) * QWm2t(5,k)&
                   + qt2(i,j,kp3) * QWm2t(6,k)
            rhst_bdf_t1(i,j,k) =  logT1(i,j,k) - w3*qbzt1
            rhst_bdf_t2(i,j,k) =  logT2(i,j,k) - w3*qbzt2

            !---rhs terms of q---
            rhsc_bdf_t1(i,j,k) = w4*qt1(i,j,k) 
            rhsc_bdf_t2(i,j,k) = w4*qt2(i,j,k) 
        end do
      end do

      do j=1, l_nj
         do i= 1, l_ni
            rhsc_bdf_t1(i,j,l_nk+1) = w4*qt1(i,j,l_nk+1)
            rhsc_bdf_t2(i,j,l_nk+1) = w4*qt2(i,j,l_nk+1)
         end do
      end do

!$OMP BARRIER

      call HLT_split (1, 12*l_nk+2, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhst_bdf_t1(Adz_lminx,Adz_lminy,HLT_start),&
                 Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)

!     ---------------------------------------------------------------
!
      return
      end subroutine rhs1
