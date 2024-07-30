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
!  New arrays added:
!	1. arrays for interpolating rhs terms to midpoint:	
!	  a. rhsq_t_mid - interpolating qt1, subsript t implies it will be used for
!				  the rhs of the Temp equation, used to compute \overline{q}^z
!	  b. rhsu_mid -
!	  c. rhsv_mid -
!	  d. rhst_mid -
!	  e. rhsf_mid -
!	  f. rhsw_mid -
!
!
!	2. arrays for interpolating rhs terms to departure point:
!	  a. rhsq_t_dep - interpolating qt2, subsript t implies it will be used for
!				  the rhs of the Temp equation, used to compute \overline{q}^z
!	  b. rhsu_mid -
!	  c. rhsv_mid -
!	  d. rhst_mid -
!	  e. rhsf_mid -
!	  f. rhsw_mid -
!
!
!	3. arrays for interpolating nonlinear terms to midoint:
!
!
!	4. arrays for interpolating nonlinear terms to departure point:
!
!     5. extra arrays to hold true rhs and nonlinear terms for T, needed in back-substitution
!        right now, the rhst and nl_t terms hold the L_T ' and Nl_t' terms for the
!        bvp S
!       a. true_rhst
!       b. true_nlt
!
!--------------------------------------------------------------------------------
module mem_tstp
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      integer :: ds_i0, ds_in, ds_j0, ds_jn, ds_k0, ds_kn
      
      real, allocatable, target, dimension (:) :: LT, var_init
 
      real, allocatable, target, dimension (:) :: rhs_mid
      real, allocatable, target, dimension (:) :: rhs_dep
      real, allocatable, target, dimension (:) :: sw_frc, sw_rhs
      
      real, dimension (:,:,:), pointer :: Ruu,Rvv,Rww,Rtt,Rzz
      real, dimension (:,:,:), pointer :: Nuu,Nvv
      
      real, dimension (:,:,:), pointer :: &
         rhsu,rhsv,rhst,rhsc,rhsw,rhsf                             ,&
         rhsu_mid, rhsv_mid, rhst_mid, rhsf_mid, rhsw_mid, rhsc_mid,&
         rhsu_dep, rhsv_dep, rhst_dep, rhsf_dep, rhsw_dep, rhsc_dep

      real, dimension (:,:,:), pointer :: sw_f1,sw_f2,sw_f3,sw_f4,&
               sw_f5,sw_f6,sw_f7,sw_f8,sw_f9,sw_f10,sw_f11,sw_f12
                                          
      real, allocatable,  dimension (:,:,:),target :: &
                    orhsu,orhsv,orhst,orhsc,orhsw,orhsf

      real, allocatable, target, dimension (:) :: &
                   orhs_extended, nl_terms, rhs_bdf

      real, dimension (:,:,:), pointer :: orhsu_ext,orhsv_ext,&
                      orhst_ext,orhsc_ext,orhsw_ext,orhsf_ext, &
                      nlu_t1, nlv_t1, nlw_t1, nlt_t1, nlq_t1,  &
                      nlu_t2, nlv_t2, nlw_t2, nlt_t2, nlq_t2,  &
                      rhsf_bdf_t1, rhst_bdf_t1, rhsc_bdf_t1,   &
                      rhsf_bdf_t2, rhst_bdf_t2, rhsc_bdf_t2, &
                      rhsu_bdf_t1, rhsv_bdf_t1, &
                      rhsu_bdf_t2, rhsv_bdf_t2, &
                      rhsw_bdf_t1, rhsw_bdf_t2
      real, allocatable,  dimension (:,:,:) :: &
                   nl_u,nl_v,nl_t,nl_c,nl_w,nl_f
      real             , allocatable, dimension (:), target :: WS1
      real(kind=REAL64), allocatable, dimension (:), target :: WS1_8

end module mem_tstp
