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

!**s/r var_topo - varies topography incrementally from analysis
!                 topography to target topography

      subroutine var_topo (F_topo,F_topo_LS, Minx,Maxx,Miny,Maxy)
      use gmm_geof
      use gem_options
      use tdpack
      use mem_nest
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension(Minx:Maxx,Miny:Maxy), intent(OUT)::F_topo,F_topo_LS

      integer i,j
      real(kind=REAL64), parameter :: one = 1.0d0
      real(kind=REAL64), parameter :: two = 2.0d0
      real(kind=REAL64)  pio2, f, b
!     __________________________________________________________________
!
      if ( .not. Vtopo_L ) then
         F_topo   (:,:) = topo_high(:,:,1)
         F_topo_LS(:,:) = topo_high(:,:,2)
         return
      end if

      f= max(0.d0,min(Nest_current_dayfrac_8-Vtopo_start_8,Vtopo_ndt_8))
      b= f / Vtopo_ndt_8

      pio2= pi_8 / two
      b = (cos(pio2 * (one-b) ))**2

      do j= 1-G_haloy, l_nj+G_haloy
      do i= 1-G_halox, l_ni+G_halox
         F_topo   (i,j) = (one - b)*topo_low(i,j,1) + b*topo_high(i,j,1)
         F_topo_LS(i,j) = (one - b)*topo_low(i,j,2) + b*topo_high(i,j,2)
      end do
      end do
!     __________________________________________________________________
!
      return
      end subroutine var_topo
