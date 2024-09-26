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
!**s/r prerhs5 - 5th order pre-compute stuff for elliptic_rhs

      subroutine prerhs5 ( F_t2u, F_v2u, F_t2v, F_u2v, F_dq2u, F_dq2v,&
                           F_dq2w, Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use gmm_vt0
      use mem_tstp
      use ver
      use vgh
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_dq2u, F_dq2v, F_dq2w

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(-2:3) :: t2qv, t2qu, v2q, u2q, dq2u, dq2v
      real(kind=REAL64) :: duu, dvv 
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) :: dqzv, dqzu
!
!     ---------------------------------------------------------------
!
      call fill_Vhalo (tt0,l_minx,l_maxx,l_miny,l_maxy,1.d0)
      call fill_Vhalo (qt0,l_minx,l_maxx,l_miny,l_maxy,1.d0)

      do k=1, l_nk
         do j= 1, l_nj
           do i= 1, l_ni
             do n=-2,3
               t2qv(n) = tt0(i+n,j,k-3) * VS5t2m(1,k) + tt0(i+n,j,k-2) * VS5t2m(2,k)&
                        +tt0(i+n,j,k-1) * VS5t2m(3,k) + tt0(i+n,j,k  ) * VS5t2m(4,k)&
                        +tt0(i+n,j,k+1) * VS5t2m(5,k) + tt0(i+n,j,k+2) * VS5t2m(6,k)

               v2q(n) = Hstag(vt0(i,j-3,k),vt0(i,j-2,k),vt0(i,j-1,k),&
                              vt0(i,j,k  ),vt0(i,j+1,k),vt0(i,j+2,k))
               t2qu(n) = tt0(i,j+n,k-3) * VS5t2m(1,k) + tt0(i,j+n,k-2) * VS5t2m(2,k)&
                        +tt0(i,j+n,k-1) * VS5t2m(3,k) + tt0(i,j+n,k  ) * VS5t2m(4,k)&
                        +tt0(i,j+n,k+1) * VS5t2m(5,k) + tt0(i,j+n,k+2) * VS5t2m(6,k)

               u2q(n) = Hstag(ut0(i-3,j,k),ut0(i-2,j,k),ut0(i-1,j,k),&
                              ut0(i,j,k  ),ut0(i+1,j,k),ut0(i+2,j,k))
               dq2u(n) = qt0(i+n,j,k-2) * VD5m2t(1,k)&
                        +qt0(i+n,j,k-1) * VD5m2t(2,k)&
                        +qt0(i+n,j,k  ) * VD5m2t(3,k)&
                        +qt0(i+n,j,k+1) * VD5m2t(4,k)&
                        +qt0(i+n,j,k+2) * VD5m2t(5,k)&
                        +qt0(i+n,j,k+3) * VD5m2t(6,k)

               dq2v(n) = qt0(i,j+n,k-2) * VD5m2t(1,k)&
                        +qt0(i,j+n,k-1) * VD5m2t(2,k)&
                        +qt0(i,j+n,k  ) * VD5m2t(3,k)&
                        +qt0(i,j+n,k+1) * VD5m2t(4,k)&
                        +qt0(i,j+n,k+2) * VD5m2t(5,k)&
                        +qt0(i,j+n,k+3) * VD5m2t(6,k)
                     end do
             F_t2u(i,j,k)= Hstag8(t2qv(-2),t2qv(-1),t2qv(0),t2qv(1),t2qv(2),t2qv(3))/Cstv_Tstr_8-1.d0
             F_v2u(i,j,k)= Hstag8( v2q(-2), v2q(-1), v2q(0), v2q(1), v2q(2), v2q(3))
             F_t2v(i,j,k)= Hstag8(t2qu(-2),t2qu(-1),t2qu(0),t2qu(1),t2qu(2),t2qu(3))/Cstv_Tstr_8-1.d0
             F_u2v(i,j,k)= Hstag8( u2q(-2), u2q(-1), u2q(0), u2q(1), u2q(2), u2q(3))
             F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0) !for Rw
             dqzu(i,j,k)= Hstag8(dq2u(-2),dq2u(-1),dq2u(0),dq2u(1),dq2u(2),dq2u (3)) !values in thermo, u-grid
             dqzv(i,j,k)= Hstag8(dq2v(-2),dq2v(-1),dq2v(0),dq2v(1),dq2v(2),dq2v (3)) !values in thermo, v-grid
         end do
         end do
      end do

      !---now interpolate back to momentum level---
      call fill_Vhalo (dqzu,l_minx,l_maxx,l_miny,l_maxy,1.d0)
      call fill_Vhalo (dqzv,l_minx,l_maxx,l_miny,l_maxy,1.d0)
      
      do k=1,l_nk
         do j= 1, l_nj
            do i= 1, l_ni
               duu = dqzu(i,j,k-3) * VS5t2m(1,k)&
                   + dqzu(i,j,k-2) * VS5t2m(2,k)&
                   + dqzu(i,j,k-1) * VS5t2m(3,k)&
                   + dqzu(i,j,k  ) * VS5t2m(4,k)&
                   + dqzu(i,j,k+1) * VS5t2m(5,k)&
                   + dqzu(i,j,k+2) * VS5t2m(6,k)

               dvv = dqzv(i,j,k-3) * VS5t2m(1,k)&
                   + dqzv(i,j,k-2) * VS5t2m(2,k)&
                   + dqzv(i,j,k-1) * VS5t2m(3,k)&
                   + dqzv(i,j,k  ) * VS5t2m(4,k)&
                   + dqzv(i,j,k+1) * VS5t2m(5,k)&
                   + dqzv(i,j,k+2) * VS5t2m(6,k)

               F_dq2u(i,j,k)= M_Jxozu(i,j,k) * duu 
               F_dq2v(i,j,k)= M_Jyozv(i,j,k) * dvv
            end do
         end do
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( F_t2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_v2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_t2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_u2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2w(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!     
!     ---------------------------------------------------------------
!
      return
      include 'H5th_ope.inc'
      end subroutine prerhs5
