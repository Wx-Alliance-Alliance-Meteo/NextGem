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
!**s/r tt2wnd - Pre-compute stuff for elliptic_rhs

      subroutine SW_prerhs3rd ( F_v2u, F_u2v, Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use gmm_vt0
      use mem_tstp
      use ver
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_v2u, F_u2v

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(-1:2) :: v2q, u2q
!
!     ---------------------------------------------------------------
!
      do k=1, l_nk
         do j= 1, l_nj
           do i= 1, l_ni
             do n=-1,2
               v2q(n) = Hstag(vt0(i+n,j-2,k),vt0(i+n,j-1,k),&
                              vt0(i+n,j,k  ),vt0(i+n,j+1,k))
               u2q(n) = Hstag(ut0(i-2,j+n,k),ut0(i-1,j+n,k),&
                              ut0(i,j+n,k  ),ut0(i+1,j+n,k))
             end do
             F_v2u(i,j,k)= Hstag8( v2q(-1), v2q(0), v2q(1), v2q(2))
             F_u2v(i,j,k)= Hstag8( u2q(-1), u2q(0), u2q(1), u2q(2))
             !F_v2u(i,j,k)= 0.25d0*(vt0(i,j,k)+vt0(i,j-1,k)+vt0(i+1,j,k)+vt0(i+1,j-1,k)) 
             !F_u2v(i,j,k)= 0.25d0*(ut0(i,j,k)+ut0(i-1,j,k)+ut0(i,j+1,k)+ut0(i-1,j+1,k))
             !F_v2u(i,j,k)= 0.d0
             !F_u2v(i,j,k)= 0.d0
           end do
         end do
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( F_v2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_u2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine SW_prerhs3rd
