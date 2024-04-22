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

!**s/r psadj_p0init  - Store air mass at initial time for Yin-Yang

      subroutine psadj_p0init (F_kount,F_done_L)
      use adz_mem
      use cstv
      use dyn_fisl_options
      use init_options
      use geomh
      use gmm_pw
      use HORgrid_options
      use mem_tstp
      use psadjust
      use masshlt
      use rstr
      use ptopo
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(INOUT) :: F_done_L
      integer, intent(in   ) :: F_kount

      include 'mpif.h'
      integer :: err,i,j,dim,dimH
      real(kind=REAL64), dimension(:,:), pointer :: p0_dry_8
      real(kind=REAL64) :: sum_lcl,gathS(Ptopo_numproc*Ptopo_ncolors)
!
!     ---------------------------------------------------------------
!
      if ((Schm_psadj==0) .or. (.not.Grd_yinyang_L) .or. Rstri_rstn_L) goto 999
      
      if (Cstv_dt_8*F_kount <= Iau_period) return

      dimH= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      p0_1_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(    1:) ; dim=dimH
      p0_0_8(l_minx:l_maxx,l_miny:l_maxy) => WS1_8(dim+1:) ; dim=dim+dimH
      thread_sum(1:2,0:OMP_max_threads-1) => WS1_8(dim+1:) ; dim= dim+2*OMP_max_threads
      g_avg_8(1:2) => WS1_8(dim+1:) ; dim= dim+2


      p0_dry_8 => p0_0_8

      !Obtain Wet surface pressure/Pressure Momentum at TIME P
      !-------------------------------------------------------

      !Obtain Surface pressure minus Cstv_pref_8
      !-----------------------------------------
      if (Schm_psadj==1) then
!!$omp do
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               p0_1_8(i,j) = pw_p0_plus_8(i,j) - Cstv_pref_8
            end do
         end do
!!$omp end do
      endif

      !Compute Dry surface pressure at TIME P (Schm_psadj==2)
      !------------------------------------------------------
      if (Schm_psadj==2) then
         call dry_sfc_pressure_hlt_8(p0_dry_8,pw_pm_plus_8,pw_p0_plus_8,Adz_k0t,'P')
!!$omp do
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               p0_1_8(i,j) = p0_dry_8(i,j) - Cstv_pref_8
            end do
         end do
!!$omp end do
      endif

      !Store air mass at initial time
      !------------------------------
      sum_lcl= 0.0d0
!!$omp do
      do j=1+pil_s,l_nj-pil_n
         do i=1+pil_w,l_ni-pil_e
            sum_lcl = sum_lcl + p0_1_8(i,j) * geomh_area_mask_8(i,j)
         end do
      end do
!!$omp end do nowait
      psadj_thread_sum(1,OMP_get_thread_num()) = sum_lcl
!$OMP BARRIER

!!$omp single
      sum_lcl = sum(psadj_thread_sum(1,:))
      call MPI_Allgather(sum_lcl,1,MPI_DOUBLE_PRECISION,gathS,1,&
                         MPI_DOUBLE_PRECISION,COMM_multigrid,err)
      PSADJ_g_avg_ps_initial_8 = sum(gathS) * PSADJ_scale_8
!!$omp end single

! Schm_psadj_print_L should always be .F. for performance purposes
      if (Schm_psadj_print_L) call stat_psadj_hlt(1,"BEFORE DYNSTEP")

 999  F_done_L= .true.
!     ---------------------------------------------------------------
!
      return
      end subroutine psadj_p0init
