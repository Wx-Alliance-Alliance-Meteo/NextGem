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
!**s/r dqdz3rd - 

      subroutine dqdz3rd ( F_q, F_dq2u , F_dq2v, F_dq2w, &
                           Minx, Maxx, Miny, Maxy, Nk,F_k0, F_kn )
      use, intrinsic :: iso_fortran_env
      use HORgrid_options
      use glb_ld
      use geomh
      use sol_mem
      use ver
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk, F_k0, F_kn
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,F_k0:F_kn), intent(IN) :: F_q
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,1:Nk), &
                               intent(OUT) :: F_dq2u, F_dq2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,-1:Nk+1),&
                               intent(OUT) :: F_dq2w

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(l_ni,l_nj) :: u,v,Jx,Jy
      real(kind=REAL64), dimension(-1:2) :: dq2u, dq2v
      real(kind=REAL64) :: dqdx, dqdy
      real(kind=REAL64), dimension(:,:,:), allocatable :: dqzv, dqzu
!
!     ---------------------------------------------------------------
!
      ! NOTE: F_q must be defined from -1:l_nk+2
      allocate (dqzu (l_minx:l_maxx,l_miny:l_maxy,-1:l_nk+1),&
                dqzv (l_minx:l_maxx,l_miny:l_maxy,-1:l_nk+1))

      do k=0, l_nk
         do j= 1, l_nj
            do i= 1, l_ni
               do n=-1,2
               dq2u(n) = (F_q(i+n,j,k-1) * VD3m2t(1,k)&
                         +F_q(i+n,j,k  ) * VD3m2t(2,k)&
                         +F_q(i+n,j,k+1) * VD3m2t(3,k)&
                         +F_q(i+n,j,k+2) * VD3m2t(4,k))*M_iJzq(i+n,j,k)
               dq2v(n) = (F_q(i,j+n,k-1) * VD3m2t(1,k)&
                         +F_q(i,j+n,k  ) * VD3m2t(2,k)&
                         +F_q(i,j+n,k+1) * VD3m2t(3,k)&
                         +F_q(i,j+n,k+2) * VD3m2t(4,k))*M_iJzq(i,j+n,k)
               end do
               F_dq2w(i,j,k) = dq2u(0)
               dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2))
               dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2))
            end do
         end do
      end do

      k= -1
         do j= 1, l_nj
            do i= 1, l_ni
               do n=-1,2
               dq2u(n) = (F_q(i+n,j,k  ) * VD3m2t(1,k)&
                         +F_q(i+n,j,k+1) * VD3m2t(2,k)&
                         +F_q(i+n,j,k+2) * VD3m2t(3,k)&
                         +F_q(i+n,j,k+3) * VD3m2t(4,k))*M_iJzq(i+n,j,k)
               dq2v(n) = (F_q(i,j+n,k  ) * VD3m2t(1,k)&
                         +F_q(i,j+n,k+1) * VD3m2t(2,k)&
                         +F_q(i,j+n,k+2) * VD3m2t(3,k)&
                         +F_q(i,j+n,k+3) * VD3m2t(4,k))*M_iJzq(i,j+n,k)
               end do
               F_dq2w(i,j,k) = dq2u(0)
               dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2))
               dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2))
            end do
         end do

      k= l_nk+1
         do j= 1, l_nj
            do i= 1, l_ni
               do n=-1,2
               dq2u(n) = (F_q(i+n,j,k-2) * VD3m2t(1,k)&
                         +F_q(i+n,j,k-1) * VD3m2t(2,k)&
                         +F_q(i+n,j,k  ) * VD3m2t(3,k)&
                         +F_q(i+n,j,k+1) * VD3m2t(4,k))*M_iJzq(i+n,j,k)
               dq2v(n) = (F_q(i,j+n,k-2) * VD3m2t(1,k)&
                         +F_q(i,j+n,k-1) * VD3m2t(2,k)&
                         +F_q(i,j+n,k  ) * VD3m2t(3,k)&
                         +F_q(i,j+n,k+1) * VD3m2t(4,k))*M_iJzq(i,j+n,k)
               end do 
               F_dq2w(i,j,k) = dq2u(0)
               dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2))
               dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2))
            end do
         end do
      !---now interpolate back to momentum level---
      do k=1, l_nk
         do j= 1, l_nj 
            do i= 1, l_ni
               u(i,j) = dqzu(i,j,k-2) * VS3t2m(1,k)&
                      + dqzu(i,j,k-1) * VS3t2m(2,k)&
                      + dqzu(i,j,k  ) * VS3t2m(3,k)&
                      + dqzu(i,j,k+1) * VS3t2m(4,k)

               v(i,j) = dqzv(i,j,k-2) * VS3t2m(1,k)&
                      + dqzv(i,j,k-1) * VS3t2m(2,k)&
                      + dqzv(i,j,k  ) * VS3t2m(3,k)&
                      + dqzv(i,j,k+1) * VS3t2m(4,k)
               Jx(i,j)= Hderiv8(VM3%zmom(i-1,j,k),VM3%zmom(i  ,j,k),&
                                VM3%zmom(i+1,j,k),VM3%zmom(i+2,j,k),geomh_invDX_8(j))
               Jy(i,j)= Hderiv8(VM3%zmom(i,j-1,k),VM3%zmom(i  ,j,k),&
                                VM3%zmom(i,j+1,k),VM3%zmom(i,j+2,k),geomh_invDY_8)
               dqdx = Hderiv8(F_q(i-1,j,k), F_q(i,j,k), &
                              F_q(i+1,j,k), F_q(i+2,j,k), geomh_invDX_8(j))
               dqdy = Hderiv8(F_q(i,j-1,k), F_q(i,j,k)  , &
                              F_q(i,j+1,k), F_q(i,j+2,k), geomh_invDY_8)
               F_dq2u(i,j,k)= dqdx - Jx(i,j) * u(i,j)
               F_dq2v(i,j,k)= dqdy - Jy(i,j) * v(i,j)
            end do
         end do
         if ( .not. Grd_yinyang_L) then
         if (l_west) then
            i=pil_w
            do j= 1, l_nj
               F_dq2u(  i,j,k)= - Jx(i,j) * u(i,j)
               F_dq2u(1+i,j,k)= (F_q(i+2,j,k)-F_q(i+1,j,k))*geomh_invDXM_8(j) - Jx(i+1,j) * u(i+1,j)
            end do
         endif
         if (l_south) then
            j=pil_s
            do i= 1, l_ni
               F_dq2v(i,  j,k)= - Jy(i,j) * v(i,j)
               F_dq2v(i,1+j,k)=  geomh_cyv_8(j+1)*(F_q(i,j+2,k)-F_q(i,j+1,k))*geomh_invDYM_8(j+1) - Jy(i,j+1) * v(i,j+1)
            end do
         endif
         if (l_east) then
            i=l_ni-pil_e
            do j= 1, l_nj
               F_dq2u(  i,j,k)= - Jx(i,j) * u(i,j)
               F_dq2u(i-1,j,k)= (F_q(i,j,k)-F_q(i-1,j,k))*geomh_invDXM_8(j) - Jx(i-1,j) * u(i-1,j)
            end do
         endif
         if (l_north) then
            j=l_nj-pil_n
            do i= 1, l_ni
               F_dq2v(i,  j,k)= - Jy(i,j) * v(i,j)
               F_dq2v(i,j-1,k)=  geomh_cyv_8(j-1)*(F_q(i,j,k)-F_q(i,j-1,k))*geomh_invDYM_8(j-1) - Jy(i,j-1) * v(i,j-1)
            end do
         endif
         endif
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( F_dq2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( F_dq2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call HLT_split (-1, l_nk+1, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( F_dq2w(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      deallocate (dqzv, dqzu)
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine dqdz3rd
