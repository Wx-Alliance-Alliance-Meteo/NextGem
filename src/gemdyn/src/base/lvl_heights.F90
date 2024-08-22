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

!**s/r lvl_heights - Compute level heights

      subroutine lvl_heights (zmom_8,ztht_8, F_topo, F_orols, Minx,Maxx,Miny,Maxy)
      use, intrinsic :: iso_fortran_env
      use gem_options
      use tdpack
      use glb_ld
      use ver
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_orols
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,0:G_nk+1), intent(OUT) :: zmom_8,ztht_8

      integer :: i,j,k
!
!     ---------------------------------------------------------------
!
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               zmom_8(i,j,k)=ver_z_8%m(k)+(Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_orols(i,j))/grav_8
               ztht_8(i,j,k)=ver_z_8%t(k)+(Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_orols(i,j))/grav_8
            end do
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            ztht_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,G_nk+1)= F_topo(i,j)/grav_8
            ztht_8(i,j,G_nk+1)= zmom_8(i,j,G_nk+1)
         end do
      end do
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine lvl_heights
