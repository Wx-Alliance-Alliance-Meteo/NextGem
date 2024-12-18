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

      subroutine lvl_heightsVH (zmom_8,ztht_8, F_topo, F_orols, &
                                Minx,Maxx,Miny,Maxy,F_k0,F_kn)
      use, intrinsic :: iso_fortran_env
      use gem_options
      use tdpack
      use glb_ld
      use ver
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,F_k0,F_kn
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_orols
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,F_k0:F_kn), intent(OUT) :: zmom_8
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,F_k0-1:F_kn), intent(OUT) :: ztht_8

      integer :: i,j,k,klm
!
!     ---------------------------------------------------------------
!
      do k=F_k0,F_kn
         klm=min(max(k,1),G_nk+1)
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               zmom_8(i,j,k)=ver_ext%m(k)+(Ver_b_8%m(klm)*F_topo(i,j)+Ver_c_8%m(klm)*F_orols(i,j))/grav_8
               ztht_8(i,j,k)=ver_ext%t(k)+(Ver_b_8%t(klm)*F_topo(i,j)+Ver_c_8%t(klm)*F_orols(i,j))/grav_8
            end do
         end do
      end do
      k= F_k0-1
      klm=min(max(k,1),G_nk+1)
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            ztht_8(i,j,k)=ver_ext%t(k)+(Ver_b_8%t(klm)*F_topo(i,j)+Ver_c_8%t(klm)*F_orols(i,j))/grav_8
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine lvl_heightsVH
