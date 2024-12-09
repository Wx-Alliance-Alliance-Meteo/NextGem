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

!**   s/r set_params - initialize some constant parameters

      subroutine set_params (F_check_and_stop_L)
      use dynkernel_options
      use dyn_fisl_options
      use tdpack
      use cstv
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
      
      logical, intent(in) :: F_check_and_stop_L

      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0
!
!     ---------------------------------------------------------------
!
      call set_dync (Cstv_dt_8)

      Ver_igt_8    = Cstv_invT_8/grav_8
      Ver_ikt_8    = Cstv_invT_m_8/cappa_8
      Ver_igt2_8   = Ver_igt_8*(Cstv_invT_nh_8/grav_8)

!
!     ---------------------------------------------------------------
!
      return
      end subroutine set_params
