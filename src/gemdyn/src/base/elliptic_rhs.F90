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
!**s/r elliptic_rhs - Compute right hand side of the elliptic problem

      subroutine elliptic_rhs ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use glb_ld
      use glb_pil
      use sol_mem
      use stat_mpi
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8
!
!     ---------------------------------------------------------------
!
      if (Schm_POSO == 5) then
         call elliptic_rhs5th ( F_dt_8, k0, k0t )
         goto 999
      else if (Schm_POSO == 3) then
         call elliptic_rhs3rd ( F_dt_8, k0, k0t )
         goto 999
      else
         call elliptic_rhs2nd ( F_dt_8, k0, k0t )
      endif
 999  call statf_dm (Sol_rhs,'RHS',1,'TSTP',1,ubound(Sol_rhs,1),&
                   1,ubound(Sol_rhs,2),1,l_nk,1+Glb_pil_w,1+Glb_pil_s,&
                   1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,l_nk,8)     
!     ---------------------------------------------------------------
!
      return
      end subroutine elliptic_rhs
