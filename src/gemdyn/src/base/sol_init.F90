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

!**s/r sol_init - initialize solver configuration
!
      integer function sol_init ()
      use lun
      use sol_mem
      use sol_options
      implicit none

      character(len=16) :: dumc_S
!
!-------------------------------------------------------------------
!
      sol_init = -1
      call low2up  (Sol_precond3D_S ,dumc_S)
      Sol_precond3D_S = trim(dumc_S)

      isol_i = 1.0d0
      isol_d = 0.0d0

      if (Sol_krylov3D_S /= 'FGMRES') then
         if (Lun_out > 0) &
         write(Lun_out, *) 'ABORT: WRONG CHOICE OF KRYLOV METHOD FOR 3D ITERATIVE SOLVER: Sol_krylov3D_S =', Sol_krylov3D_S
         return
      endif
      sol_init = 0
!
!-------------------------------------------------------------------
!
      return
      end function sol_init

