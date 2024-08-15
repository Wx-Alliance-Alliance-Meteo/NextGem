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
!	same as adz_main_h, expect the source for the interpolation is different.
!	Trajectories code and staggering is the same.
!
!
!------------------------------------------------------------------------------
      subroutine SW_adz_main (F_dt_8, itpc, F_euler_L)
      use ISO_C_BINDING
      use dyn_fisl_options
      use cstv
      use gem_options
      use adz_mem
      use adz_interp_mod
      use dynkernel_options
      use mem_tstp
      use gmm_vt1
      use gmm_vt2
      use step_options
      implicit none

      logical, intent(IN) :: F_euler_L
      integer, intent(IN) :: itpc
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: n,i,j,k
      integer :: HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8
      type(Adz_pntr_stack), dimension(3), target :: stack
!
!     ---------------------------------------------------------------
!
      if(F_euler_L) then
         call SW_euler_adz_traject (F_dt_8,itpc) !<--bdf trajectory
      else
         call SW_adz_traject (F_dt_8,itpc) !<--bdf trajectory
      endif

!==========RHS terms======

!---midpoint---
      stack(1)%src => orhsu_ext
      stack(1)%dst => rhsu_mid
      call adz_tricub ( stack,1,Adz_pmu,Adz_cpntr_q,Adz_num_u,&
                            Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )
      stack(1)%src => orhsv_ext
      stack(1)%dst => rhsv_mid
      call adz_tricub ( stack,1,Adz_pmv,Adz_cpntr_q,Adz_num_v,&
                            Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )
      stack(1)%src => orhsw_ext
      stack(1)%dst => rhsc_mid
      call adz_tricub ( stack,1,Adz_pm ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

!---departure---
      stack(1)%src => orhst_ext
      stack(1)%dst => rhsu_dep
      call adz_tricub ( stack,1,Adz_pdu,Adz_cpntr_q,Adz_num_u,&
                            Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )
      stack(1)%src => orhsc_ext
      stack(1)%dst => rhsv_dep
      call adz_tricub ( stack,1,Adz_pdv,Adz_cpntr_q,Adz_num_v,&
                            Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )
      stack(1)%src => orhsf_ext
      stack(1)%dst => rhsc_dep
      call adz_tricub ( stack,1,Adz_dep ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

!==========NL terms========
!---no advection for NL terms with Picard iterations---

!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_mid(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_dep(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

!     ---------------------------------------------------------------
!
      return
      end subroutine SW_adz_main
