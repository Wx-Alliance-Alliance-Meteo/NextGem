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
      use glb_ld
      use ctrl
      use dyn_fisl_options
      use gmm_vt0
      use HORgrid_options
      use step_options
      use sol_mem
      use mem_tstp
      use metric
      use omp_timing
      use gmm_geof
      use ptopo
      use glb_pil
      use ldnh
      use stat_mpi
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8

      logical, external :: picard_stop
      logical :: print_conv, first_time_L=.true.
      integer i,j,k, itpc, iter, ni,nj
      integer :: HLT_np, HLT_start, HLT_end
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
      real(kind=REAL64) :: dt_8, invT_m_8
!
!     ---------------------------------------------------------------
!
      ni=ldnh_maxx-ldnh_minx+1
      nj=ldnh_maxy-ldnh_miny+1
      
      print_conv = (Ptopo_couleur == 0  ) .and. (Lun_out > 0)

      dt_8 = F_dt_8
      first_time_L= (Step_kount.le.1).and.first_time_L
      
      call gtmg_start (20, 'TSTPDYN', 10)

      call HLT_split (1, 6*(l_nk+3), HLT_np, HLT_start, HLT_end)

      if(Step_kount.le.2) then
         call set_dync ( .true., dt_8 )
         call vertical_metric ()
      endif

      call rhs1 (dt_8) ! also first guess into t0

      if ( .not. Grd_yinyang_L ) then
         call nest_bcs (dt_8,Ruu,Rvv,l_minx,l_maxx,l_miny,l_maxy,l_nk)
      end if

!-----Begin Picard Iterations
    
      do iter = 1, Schm_itpc
         
        ! call fill_Vhalo (wt0,l_minx,l_maxx,l_miny,l_maxy,G_nk+3,6)
         call gem_xch_halo ( wt0(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

         call gtmg_start (25, 'ADVECTION', 20)
         call adz_main (dt_8, iter, first_time_L, Schm_advec_type_S)
         call gtmg_stop (25)

         call gtmg_start (27, 'PRE', 20)
         
         call oro_adj ()
         
         call elliptic_rhs (dt_8, ds_k0, ds_k0)
         
         call gtmg_stop (27)

         if (.not.ctrl_testcases_adv_L) then
            call gtmg_start (29, 'SOL', 20)
!call statf_dm (Sol_rhs, 'RHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,l_nk,8)

            call sol_fgmres (print_conv)
            
!call statf_dm (Sol_lhs, 'LHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,l_nk,8)
            call gtmg_stop (29)
         endif

         call gtmg_start (30, 'BAC', 20)
         call bac (dt_8)
         call gtmg_stop (30)

         if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut0, vt0,&
            l_minx,l_maxx,l_miny,l_maxy,G_nk)
            call yyg_xchng_hlt (tt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (zdt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (qt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
            G_nk+1, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (wt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
            G_nk, .false., 'CUBIC', .false.)
         end if

!         call blocstat (.false.)

!---  check convergence---
         if (picard_stop(dt_8,iter,print_conv)) exit

      enddo
      
      first_time_L= .false.                      
      call gtmg_stop (20)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine tstpdyn
