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

!**s/r fill_Vhalo - Fill vertical bottom halo

      subroutine fill_Vhalo (F_dst,Minx,Maxx,Miny,Maxy,NK,F_nf)
      use gem_options
      use glb_ld
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,NK,F_nf
      real, dimension(Minx:Maxx,Miny:Maxy,NK,F_nf), intent(INOUT) :: F_dst
      integer :: i,j,n
!
!     ---------------------------------------------------------------
!
      do n=1,F_nf
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            F_dst(i,j,G_nk+1,n) = F_dst(i,j,G_nk  ,n)
            F_dst(i,j,G_nk+2,n) = F_dst(i,j,G_nk-1,n)
            F_dst(i,j,G_nk+3,n) = F_dst(i,j,G_nk-2,n)
         end do
      end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine fill_Vhalo
