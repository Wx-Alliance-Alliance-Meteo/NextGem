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

!**s/r psadj_init - Estimate area for Yin-Yang/LAM

      subroutine psadj_init ()
      use glb_ld
      use dyn_fisl_options
      use geomh
      use gmm_geof
      use psadjust
      use ptopo
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none

      include 'mpif.h'
      integer :: err,i,j

      real(kind=REAL64) :: sum_lcl(2),l_avg_8(2),gathS(2,Ptopo_numproc*Ptopo_ncolors)
      logical :: almost_zero
!
!     ---------------------------------------------------------------
!
      sum_lcl(1)=0.d0 ; sum_lcl(2)=0.d0
!!$omp do
      do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            sum_lcl(1)= sum_lcl(1) + geomh_area_mask_8(i,j)
            if (fis0(i,j) > 1.) sum_lcl(2) = sum_lcl(2) + geomh_area_mask_8(i,j)
         end do
      end do
!!$omp end do nowait
      psadj_thread_sum(1,OMP_get_thread_num()) = sum_lcl(1)
      psadj_thread_sum(2,OMP_get_thread_num()) = sum_lcl(2)
!$OMP BARRIER

!!$omp single
      l_avg_8(1) = sum(psadj_thread_sum(1,:))
      l_avg_8(2) = sum(psadj_thread_sum(2,:))
      call MPI_Allgather(l_avg_8,2,MPI_DOUBLE_PRECISION,gathS,2,MPI_DOUBLE_PRECISION,COMM_multigrid,err)
      PSADJ_scale_8 = 1.0d0/sum(gathS(1,:))
      PSADJ_fact_8  = sum(gathS(2,:)) * PSADJ_scale_8
      if (.not.almost_zero(PSADJ_fact_8)) PSADJ_fact_8 = (1.0d0-(1.0d0-PSADJ_fact_8)*Cstv_psadj_8)/PSADJ_fact_8
!!$omp end single
!
!     ---------------------------------------------------------------
!
      return
      end subroutine psadj_init
