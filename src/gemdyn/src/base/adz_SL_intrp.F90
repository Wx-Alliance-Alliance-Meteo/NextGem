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
      subroutine adz_SL_intrp (F_order)
      use ISO_C_BINDING
      use dyn_fisl_options
      use adz_mem
      use SL_interp_mod
      use mem_tstp
      implicit none

      integer, intent(IN) :: F_order
      
      type(Adz_pntr_stack), dimension(3), target :: stack
!
!     ---------------------------------------------------------------
!
!---midpoint---
      stack(1)%src => rhsu_bdf_t1
      stack(1)%dst => rhsu_mid
      call SL_interp ( stack,1, Adz_pmu ,Adz_cpntr_q, Adz_num_u,&
                       Adz_i0u,Adz_inu,Adz_j0,Adz_jn,1,l_nk,F_order,'m')

      stack(1)%src => rhsv_bdf_t1
      stack(1)%dst => rhsv_mid
      call SL_interp ( stack,1, Adz_pmv ,Adz_cpntr_q, Adz_num_v,&
                       Adz_i0,Adz_in,Adz_j0v,Adz_jnv,1,l_nk,F_order,'m')

      stack(1)%src => rhsc_bdf_t1
      stack(1)%dst => rhsc_mid
      if (SL_sfc) then
         stop 'adz_SL_intrp'
      else
         call SL_interp ( stack,1, Adz_pm ,Adz_cpntr_q, Adz_num_q,&
                          Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'m')
      endif

      stack(1)%src => rhst_bdf_t1 
      stack(1)%dst => rhst_mid
      call SL_interp ( stack,1, Adz_pt ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')
      stack(1)%src => rhsf_bdf_t1
      stack(1)%dst => rhsf_mid
      call SL_interp ( stack,1, Adz_pt ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')
      stack(1)%src => rhsw_bdf_t1
      stack(1)%dst => rhsw_mid
      call SL_interp ( stack,1, Adz_pt ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')

!---departure---
      stack(1)%src => rhsu_bdf_t2
      stack(1)%dst => rhsu_dep
      call SL_interp ( stack,1, Adz_pdu ,Adz_cpntr_q, Adz_num_u,&
                       Adz_i0u,Adz_inu,Adz_j0,Adz_jn,1,l_nk,F_order,'m')

      stack(1)%src => rhsv_bdf_t2
      stack(1)%dst => rhsv_dep
      call SL_interp ( stack,1, Adz_pdv ,Adz_cpntr_q, Adz_num_v,&
                       Adz_i0,Adz_in,Adz_j0v,Adz_jnv,1,l_nk,3,'m')

      stack(1)%src => rhsc_bdf_t2
      stack(1)%dst => rhsc_dep
      if (SL_sfc) then
         stop 'adz_SL_intrp'
      else
         call SL_interp ( stack,1, Adz_dep ,Adz_cpntr_q, Adz_num_q,&
                          Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'m')
      endif
                           
      stack(1)%src => rhst_bdf_t2
      stack(1)%dst => rhst_dep
      call SL_interp ( stack,1, Adz_pt2 ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')
      stack(1)%src => rhsf_bdf_t2
      stack(1)%dst => rhsf_dep
      call SL_interp ( stack,1, Adz_pt2 ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')
      stack(1)%src => rhsw_bdf_t2
      stack(1)%dst => rhsw_dep
      call SL_interp ( stack,1, Adz_pt2 ,Adz_cpntr_t, Adz_num_t,&
                       Adz_i0,Adz_in,Adz_j0,Adz_jn,1,l_nk,F_order,'t')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_SL_intrp
