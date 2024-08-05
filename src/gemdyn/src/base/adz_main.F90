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
      subroutine adz_main (F_dt_8, itpc, F_euler_L)
      use ISO_C_BINDING
      use dyn_fisl_options
      use cstv
      use gem_options
      use adz_mem
      use adz_interp_hlt_mod
      use SL_interp_mod
      use dynkernel_options
      use mem_tstp
      use gmm_vt1
      use gmm_vt2
      use step_options
      implicit none

      logical, intent(IN) :: F_euler_L
      integer, intent(IN) :: itpc
      real(kind=REAL64), intent(IN) :: F_dt_8
      
      integer :: n,i,j,k, HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8
      type(Adz_pntr_stack), dimension(3), target :: stack
!
!     ---------------------------------------------------------------
!
      if (F_euler_L) then
         call euler_adz_trapezoid_pic_stop (F_dt_8)
      else
         call adz_traject_pic_stop (F_dt_8,itpc) 
      endif

!---midpoint---
      stack(1)%src => rhsu_bdf_t1
      stack(1)%dst => rhsu_mid
      !call adz_tricub_hlt ( stack,1,Adz_pmu,Adz_cpntr_q,Adz_num_u,&
      !                      Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )

      call SL_interp ( stack,1, Adz_pmu ,Adz_cpntr_q, Adz_num_u,&
                       Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0,l_nk )

      stack(1)%src => rhsv_bdf_t1
      stack(1)%dst => rhsv_mid
      !call adz_tricub_hlt ( stack,1,Adz_pmv,Adz_cpntr_q,Adz_num_v,&
      !                      Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )

      call SL_interp ( stack,1, Adz_pmv ,Adz_cpntr_q, Adz_num_v,&
                       Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0,l_nk )


      stack(1)%src => rhsc_bdf_t1
      stack(1)%dst => rhsc_mid
      if (SL_sfc) then
         call SL_interp ( stack,1, Adz_pm ,Adz_cpntr_q, Adz_num_q,&
            Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk+1 )
      else
         call adz_tricub_hlt ( stack,1,Adz_pm ,Adz_cpntr_q,Adz_num_q,&
                               Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
      endif

      stack(1)%src => rhst_bdf_t1 
      stack(1)%dst => rhst_mid
      stack(2)%src => rhsf_bdf_t1
      stack(2)%dst => rhsf_mid

      stack(3)%src => rhsw_bdf_t1
      stack(3)%dst => rhsw_mid
      call adz_tricub_hlt ( stack,3,Adz_pt ,Adz_cpntr_t,Adz_num_t,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t )

      !call SL_interp ( stack,3, Adz_pt ,Adz_cpntr_t, Adz_num_t,&
      !                 Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t,l_nk )


!---departure---
      stack(1)%src => rhsu_bdf_t2
      stack(1)%dst => rhsu_dep
      !call adz_tricub_hlt ( stack,1,Adz_pdu,Adz_cpntr_q,Adz_num_u,&
      !                      Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )

      call SL_interp ( stack,1, Adz_pdu ,Adz_cpntr_q, Adz_num_u,&
                       Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0,l_nk )



      stack(1)%src => rhsv_bdf_t2
      stack(1)%dst => rhsv_dep
      !call adz_tricub_hlt ( stack,1,Adz_pdv,Adz_cpntr_q,Adz_num_v,&
      !                      Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )

      call SL_interp ( stack,1, Adz_pdv ,Adz_cpntr_q, Adz_num_v,&
                       Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0,l_nk )


      stack(1)%src => rhsc_bdf_t2
      stack(1)%dst => rhsc_dep

      if (SL_sfc) then
         call SL_interp ( stack,1, Adz_dep ,Adz_cpntr_q, Adz_num_q,&
              Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk+1 )
      else
         call adz_tricub_hlt ( stack,1,Adz_dep ,Adz_cpntr_q,Adz_num_q,&
                              Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )
      endif
                           
      stack(1)%src => rhst_bdf_t2
      stack(1)%dst => rhst_dep
      stack(2)%src => rhsf_bdf_t2
      stack(2)%dst => rhsf_dep
      stack(3)%src => rhsw_bdf_t2
      stack(3)%dst => rhsw_dep
      call adz_tricub_hlt ( stack,3,Adz_pt2 ,Adz_cpntr_t,Adz_num_t,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t )

      !call SL_interp ( stack,3, Adz_pt2 ,Adz_cpntr_t, Adz_num_t,&
      !                 Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t,l_nk )



!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_mid(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_dep(l_minx,l_miny,HLT_start),&
      l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_main
