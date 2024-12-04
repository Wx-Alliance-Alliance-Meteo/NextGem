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

!**   s/r set_precon - initialize preconditionner parameters 

      subroutine set_precon ()
      use dynkernel_options
      use glb_ld
      use sol_mem
      use prec
      use, intrinsic :: iso_fortran_env
      implicit none

!     ---------------------------------------------------------------

      if(Dynamics_sw_L) then
           call SW_set_oprz ()
      else
           call set_oprz ()
      endif

      call eigenabc_local (Prec_xeval_8,Prec_xevec_8,Prec_ai_8,&
                           Prec_bi_8,Prec_invbi_8,Prec_ci_8   ,&
                           sol_niloc,sol_njloc,G_nk)
!
!     ---------------------------------------------------------------
!
      return
      end
