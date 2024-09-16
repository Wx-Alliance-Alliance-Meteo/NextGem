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
!*s/r SW_rhs1 - compute the right-hand sides needed for interpolation
!               of bdf portion

      subroutine SW_rhs1 ( F_dt_8 )
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
      use gmm_vt1
      use gmm_vt2
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

!NOTE: with the new derivation of bdf, I don't need to
!compute a rhs term for T or q, only z
!however in the mean time, instead of erasing the arrays,
!I will just set them to 0            

      w3=half/(cpd_8*Cstv_Tstr_8)
      w4=epsi_8/grav_8

      !print *,"t1: |  rhst  |  rhsf  |  rhsc  |"

      print *,"---inside bdf rhs---"

      do n= HLT_start, HLT_end
         k= (n-1)/HLT_nj
         j= n - k*HLT_nj + HLT_j0 - 1
         k= k+1

         do i= 1, l_ni

            orhsu_ext(i,j,k)  = ut1(i,j,k)
            orhsv_ext(i,j,k)  = vt1(i,j,k)
            orhsw_ext(i,j,k)  = qt1(i,j,k)-fis0(i,j)
            orhst_ext(i,j,k)  = ut2(i,j,k)
            orhsc_ext(i,j,k)  = vt2(i,j,k)
            orhsf_ext(i,j,k)  = qt2(i,j,k)-fis0(i,j)

        end do
      end do
      
!!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo (  orhsu_ext(Adz_lminx,Adz_lminy,HLT_start),&
                Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)


!     
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_rhs1
