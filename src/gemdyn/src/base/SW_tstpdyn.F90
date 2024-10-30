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

!**s/r SW_tstpdyn -  Performs a dynamical timestep of the model

      subroutine SW_tstpdyn ( F_dt_8 )
      use adz_mem
      use cstv
      use dyn_fisl_options
      use gmm_vt0
      use gmm_vt1
      use gmm_vt2
      use HORgrid_options
      use lam_options
      use step_options
      use ldnh
      use sol_mem
      use mem_tstp
      use mem_nest
      use omp_timing
      use glb_ld
      use gmm_pw


      use gmm_geof
      use gmm_pw
      use inp_mod
      use inp_options
      use lun
      use theo_options
      use tr3d
      use mem_tracers
      use rstr
      use mem_nest
      use init_options
      use gem_options
      use wil_options
      use ptopo

      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8

      include 'mpif.h'
      logical, external :: picard_stop
      logical :: print_conv, first_time_L=.true.
      integer i0, in, j0, jn, k0, k0t, ni, nj, iln, icln
      integer i,j,k
      real(kind=REAL64) :: dt_8, invT_m_8, inf_err,glb
      integer :: HLT_np, HLT_start, HLT_end, itpc, local_np, ierr
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
      ni = ldnh_maxx-ldnh_minx+1
      nj = ldnh_maxy-ldnh_miny+1
      
      dt_8 = F_dt_8
      first_time_L= (Step_kount.le.1).and.first_time_L

      call gtmg_start (20, 'TSTPDYN', 10)

      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut0(l_minx,l_miny,1), vt0(l_minx,l_miny,1),&
                                   l_minx,l_maxx,l_miny,l_maxy,G_nk)
      endif

!     call HLT_split (1, 6*l_nk+2, HLT_np, HLT_start, HLT_end)
      call HLT_split (-2, 6*(G_nk+6)-3, HLT_np, HLT_start, HLT_end)

      call set_dync ( .true., dt_8 )

      !save the initialcondition
      !print *,"Step kount = ", Step_kount,"Step slopi = ", Step_slopi_ini 
      if (Step_kount == 1 .and. Step_slopi_ini == 1) then 
        
        do k = k0, l_nk
          do j = j0, jn
            do i = i0, in
               initq(i,j,k) = qt0(i,j,k)
               !print *, i,j,k,initq(i,j,k), qt0(i,j,k)
            end do
          end do
        end do

      end if 


!2.	Compute bdf terms that will be on rhs for current and previous time levels
      call SW_rhs1(dt_8)

      if ( .not. Grd_yinyang_L ) then
         call nest_bcs (dt_8,Ruu,Rvv,l_minx,l_maxx,l_miny,l_maxy,l_nk)
      end if

      do itpc=1, Schm_itpc

         call gem_xch_halo ( wt0(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

         !4. Perform Semi-Lagrangian advection using standard 3-level displacement
         call gtmg_start (25, 'ADVECTION', 20)
         call SW_adz_main (dt_8,itpc,first_time_L)
         call gtmg_stop (25)

         call gtmg_start (27, 'PRE', 20)

    
         !5. Form rhs using the bdf terms
         call SW_elliptic_rhs (dt_8, k0, k0t)

         print_conv = .TRUE.
         call gtmg_start (29, 'SOL', 20)
!!$!$omp single
!!$      call statf_dm (Sol_rhs, 'RHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1,1,1,G_ni,G_nj,l_nk,8)
!!$!$omp end single

         call sol_fgmres (print_conv)
         
!!$!$omp single
!!$      call statf_dm (Sol_lhs, 'LHS', 1, 'TSTP', 1,ni,1,nj,1,l_nk,1,1,1,G_ni,G_nj,l_nk,8)
!!$!$omp end single
         call gtmg_stop (25)

         !10.  Back subtitution; same back sub!
         call gtmg_start (30, 'BAC', 20)
         call SW_bac (dt_8, i0, j0, k0, in, jn, k0t)
         call gtmg_stop (30)

         if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut0(l_minx,l_miny,1), vt0(l_minx,l_miny,1),&
                                      l_minx,l_maxx,l_miny,l_maxy,G_nk)
            call yyg_xchng_hlt (tt0(l_minx,l_miny,1) , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (zdt0(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (qt0(l_minx,l_miny,1) , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk+1, .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (wt0(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
         end if

         !---  check convergence---
         if (picard_stop(dt_8,itpc,print_conv)) exit


      enddo !Picard iter

      Step_slopi_ini = Step_slopi_ini + 1                

      first_time_L= .false.
      call gtmg_stop (20)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_tstpdyn
