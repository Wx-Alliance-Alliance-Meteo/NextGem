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

      subroutine tt2wnd ( F_t2u, F_v2u, F_t2v, F_u2v, Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use gmm_vt0
      use mem_tstp
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v

      integer :: i, j, k, n
      integer :: km1,km2,km3,kp1,kp2,kp3,i00,inn,j00,jnn
      real(kind=REAL64), dimension(-2:3) :: t2q, v2q, u2q
!
!     ---------------------------------------------------------------
!
      i00= ds_i0-1 ; inn= ds_in
      j00= ds_j0-1 ; jnn= ds_jn
      if (.not.Grd_yinyang_L) then
         i00= ds_i0-1+min(pil_w,1)
         inn= ds_in  -min(pil_e,1)
         j00= ds_j0-1+min(pil_s,1)
         jnn= ds_jn  -min(pil_n,1)
      endif
      
      do k=1, l_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         do j= ds_j0, ds_jn
         do i= i00, inn
            do n=-2,3
            t2q(n) = tt0(i+n,j,km3) * QWt2m(1,k) + tt0(i+n,j,km2) * QWt2m(2,k)&
                    +tt0(i+n,j,km1) * QWt2m(3,k) + tt0(i+n,j,k  ) * QWt2m(4,k)&
                    +tt0(i+n,j,kp1) * QWt2m(5,k) + tt0(i+n,j,kp2) * QWt2m(6,k)
            v2q(n) = Hstag(vt0(i,j-3,k),vt0(i,j-2,k),vt0(i,j-1,k),&
                           vt0(i,j,k  ),vt0(i,j+1,k),vt0(i,j+2,k))
            end do
            F_t2u(i,j,k)= Hstag8(t2q(-2),t2q(-1),t2q(0),t2q(1),t2q(2),t2q(3))/Cstv_Tstr_8-1.d0
            F_v2u(i,j,k)= Hstag8(v2q(-2),v2q(-1),v2q(0),v2q(1),v2q(2),v2q(3))
         end do
         end do
         do j= j00, jnn
         do i= ds_i0, ds_in
            do n=-2,3
            t2q(n) = tt0(i,j+n,km3) * QWt2m(1,k) + tt0(i,j+n,km2) * QWt2m(2,k)&
                    +tt0(i,j+n,km1) * QWt2m(3,k) + tt0(i,j+n,k  ) * QWt2m(4,k)&
                    +tt0(i,j+n,kp1) * QWt2m(5,k) + tt0(i,j+n,kp2) * QWt2m(6,k)
            u2q(n) = Hstag(ut0(i-3,j,k),ut0(i-2,j,k),ut0(i-1,j,k),&
                           ut0(i,j,k  ),ut0(i+1,j,k),ut0(i+2,j,k))
            F_t2v(i,j,k)= Hstag8(t2q(-2),t2q(-1),t2q(0),t2q(1),t2q(2),t2q(3))/Cstv_Tstr_8-1.d0
            F_u2v(i,j,k)= Hstag8(u2q(-2),u2q(-1),u2q(0),u2q(1),u2q(2),u2q(3))
           end do
         end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      include 'H5th_ope.inc'
      end subroutine tt2wnd
