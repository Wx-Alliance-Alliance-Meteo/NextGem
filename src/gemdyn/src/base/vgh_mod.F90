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

module vgh
  use, intrinsic :: iso_fortran_env
  use gem_options
  use glb_ld
  implicit none
  private
  public :: fill_Vhalo

  interface fill_Vhalo
     module procedure fill_Vhalo_r4
     module procedure fill_Vhalo_r8
  end interface

!object
!        Fill vertical ghosts halo

contains

      subroutine fill_Vhalo_r4 (F_dst,Minx,Maxx,Miny,Maxy,F_cte)
      implicit none
      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real(kind=REAL64), intent(IN) :: F_cte
      real, dimension(Minx:Maxx,Miny:Maxy,-2:l_nk+3), intent(INOUT) :: F_dst

      include 'fill_Vhalo.inc'
!
!----------------------------------------------------------------
!
      return
      end subroutine fill_Vhalo_r4

      subroutine fill_Vhalo_r8 (F_dst,Minx,Maxx,Miny,Maxy,F_cte)
      implicit none
      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real(kind=REAL64), intent(IN) :: F_cte
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,-2:l_nk+3), intent(INOUT) :: F_dst

      include 'fill_Vhalo.inc'
!
!----------------------------------------------------------------
!
      return
      end subroutine fill_Vhalo_r8

end module vgh

