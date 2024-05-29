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

!**s/r tstpdyn -  Performs a dynamical timestep of the model

      subroutine tstpdyn ( F_dt_8 )
      use, intrinsic :: iso_fortran_env
      use adz_mem
      use cstv
      use ctrl
      use dyn_fisl_options
      use gmm_vt0
      use gmm_vt1
      use gmm_vt2
      use HORgrid_options
      use lam_options
      use step_options
      use sol_mem
      use mem_tstp
      use metric
      use omp_timing
      use gmm_pw
      use gmm_geof
      use gem_options
      use ptopo
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      
      logical :: print_conv, first_time_L=.true.
      integer i,j,k,i0, in, j0, jn, k0, k0t, itpc, iter, io
      integer :: HLT_np, HLT_start, HLT_end!HLT_start2, lcl2, HLT_start5, lcl5, HLT_end
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
      real(kind=REAL64) :: dt_8, invT_m_8
      !-----picard stopping criteria-----
      integer, parameter :: max_iter = 10
      real(kind=REAL64), parameter :: tol=1e-6 !tol in single precision
      real(kind=REAL64), dimension(11) :: err !error defined by bdf method with infinity-norm
      real(kind=REAL64), dimension(11) :: sucsv_err !error of successive iterates for inf-nrm
      real(kind=REAL64), dimension(11) :: err2nrm !error defined by bdf method with 2-norm
      real(kind=REAL64), dimension(11) :: sucsv_2nrm_err !error defined by bdf method
!     
!     ---------------------------------------------------------------
!
      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
      k0= 1+Lam_gbpil_T
      k0t=k0
      if (Schm_opentop_L) k0t=k0-1
      print_conv = (Ptopo_couleur == 0  ) .and. (Lun_out > 0)

      dt_8 = F_dt_8
      first_time_L= (Step_kount.le.1).and.first_time_L
      
      call gtmg_start (20, 'TSTPDYN', 10)

      call HLT_split (1, 6*l_nk+2, HLT_np, HLT_start, HLT_end)

      if(Step_kount.le.2) then
         call set_dync ( .true., dt_8 )
         call vertical_metric_omp (GVM, fis0, orols,&
                       l_minx,l_maxx,l_miny,l_maxy)
      endif

      call rhs1 (dt_8) ! first guess into t0

      if ( .not. Grd_yinyang_L ) then
         call nest_bcs (dt_8,rhsu,rhsv,l_minx,l_maxx,l_miny,l_maxy,l_nk)
      end if

      err(1:11) = 0.d0
      sucsv_err = 0.d0

      err2nrm        = 0.d0
      sucsv_2nrm_err = 0.d0

!-----Begin Picard Iterations
    
      do iter = 1, max_iter 

      if (Lun_out > 0) print *, "IT_PICARD:", iter

      call gem_xch_halo ( wt0(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      call gtmg_start (25, 'ADVECTION', 20)
      call adz_main (dt_8, iter, first_time_L)
      call gtmg_stop (25)

      call gtmg_start (27, 'PRE', 20)
      call oro_adj ()

      call elliptic_rhs (dt_8, k0, k0t)

      call rhs2 (dt_8, k0, k0t )

      call pre (dt_8, i0, j0, k0, in, jn, k0t )
      call gtmg_stop (27)

      call gtmg_start (28, 'NLI', 20)
      call nli (dt_8, i0, j0, k0, in, jn, k0t)
      call gtmg_stop (28)
      
      if (.not.ctrl_testcases_adv_L) then
      call gtmg_start (29, 'SOL', 20)
      !call statf_dm (Sol_rhs, 'RHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1,1,1,G_ni,G_nj,l_nk,8)

      call sol_fgmres (print_conv)
         
      !call statf_dm (Sol_lhs, 'LHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1,1,1,G_ni,G_nj,l_nk,8)
      call gtmg_stop (25)
      endif

!10.  Back subtitution; same back sub!
      call gtmg_start (30, 'BAC', 20)
      call bac (dt_8, i0, j0, k0, in, jn, k0t)
      call gtmg_stop (30)

      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut0, vt0,&
                                   l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call yyg_xchng_hlt (tt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
         call yyg_xchng_hlt (zdt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
         call yyg_xchng_hlt (st0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         1,    .false., 'CUBIC', .false.)
         call yyg_xchng_hlt (qt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk+1, .false., 'CUBIC', .false.)
         call yyg_xchng_hlt (wt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
      end if

      !---check convergence---
      call check_picard_stop(dt_8, err(iter), err2nrm(iter))

      if (print_conv) then
         print *,"PICARD L_inf ITERATION ERROR = ", err(iter)
         print *,"PICARD L_2   ITERATION ERROR = ", err2nrm(iter)
      endif

      if (err2nrm(iter) < tol) then 
         exit
      endif

      if (iter > 1) then
        sucsv_err(iter)      = abs(err(iter) - err(iter-1))
        sucsv_2nrm_err(iter) = abs(err2nrm(iter) - err2nrm(iter-1))

        if (print_conv) then
        print *,"  "
        print *,"L_inf: ERR(ITER) - ERR(PREV-ITER) = ", sucsv_err(iter)
        print *,"L_2  : ERR(ITER) - ERR(PREV-ITER) = ", sucsv_2nrm_err(iter)
        endif 
        if (sucsv_2nrm_err(iter) < err2nrm(iter)) exit
      endif

      enddo                     !Picard iter.
      
      first_time_L= .false.                      
      call gtmg_stop (20)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine tstpdyn
