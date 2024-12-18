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
!**s/r prerhs3rd - Pre-compute stuff for elliptic_rhs

      subroutine prerhs3rdVH ( F_t2u, F_v2u, F_t2v, F_u2v, &
               F_dq2u, F_dq2v, F_dq2w, Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use glb_ld
      use gmm_geof
      use gmm_vt0
      use mem_tstp
      use tdpack
      use ver
      use vgh
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_dq2u, F_dq2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,-1:Nk+1), intent(OUT) :: F_dq2w

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(-1:2) :: t2qv, t2qu, v2q, u2q, dq2u, dq2v
      real(kind=REAL64) :: duu, dvv ,iJzq
      real(kind=REAL64), dimension(:,:,:), allocatable :: dqzv, dqzu, ext_t, ext_q
!
!     ---------------------------------------------------------------
!
      allocate (dqzu (l_minx:l_maxx,l_miny:l_maxy,-1:l_nk+1),&
                dqzv (l_minx:l_maxx,l_miny:l_maxy,-1:l_nk+1),&
                ext_t(l_minx:l_maxx,l_miny:l_maxy,-1:l_nk+1),&
                ext_q(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3))

      ext_t(:,:,1:l_nk) = tt0(:,:,1:l_nk)
      ext_q(:,:,1:l_nk) = qt0(:,:,1:l_nk)
      call fill_Vhalo (ext_t,l_minx,l_maxx,l_miny,l_maxy,lbound(ext_t,3),ubound(ext_t,3),1.d0)
      call fill_Vhalo (ext_q,l_minx,l_maxx,l_miny,l_maxy,lbound(ext_q,3),ubound(ext_q,3),1.d0)

      do k=1, l_nk
         do j= 1, l_nj
           do i= 1, l_ni
             do n=-1,2
                t2qv(n) = ext_t(i+n,j,k-2)*VS3t2m(1,k) + ext_t(i+n,j,k-1)*VS3t2m(2,k) &
                        + ext_t(i+n,j,k  )*VS3t2m(3,k) + ext_t(i+n,j,k+1)*VS3t2m(4,k)
               v2q(n) = Hstag(vt0(i+n,j-2,k),vt0(i+n,j-1,k),&
                              vt0(i+n,j,k  ),vt0(i+n,j+1,k))
               t2qu(n) = ext_t(i,j+n,k-2)*VS3t2m(1,k) + ext_t(i,j+n,k-1)*VS3t2m(2,k) &
                       + ext_t(i,j+n,k  )*VS3t2m(3,k) + ext_t(i,j+n,k+1)*VS3t2m(4,k)
               u2q(n) = Hstag(ut0(i-2,j+n,k),ut0(i-1,j+n,k),&
                              ut0(i,j+n,k  ),ut0(i+1,j+n,k))
               dq2u(n) = ext_q(i+n,j,k-1) * VD3m2t(1,k)&
                        +ext_q(i+n,j,k  ) * VD3m2t(2,k)&
                        +ext_q(i+n,j,k+1) * VD3m2t(3,k)&
                        +ext_q(i+n,j,k+2) * VD3m2t(4,k)
               dq2v(n) = ext_q(i,j+n,k-1) * VD3m2t(1,k)&
                        +ext_q(i,j+n,k  ) * VD3m2t(2,k)&
                        +ext_q(i,j+n,k+1) * VD3m2t(3,k)&
                        +ext_q(i,j+n,k+2) * VD3m2t(4,k)
             end do
             F_t2u(i,j,k)= Hstag8(t2qv(-1),t2qv(0),t2qv(1),t2qv(2))/Cstv_Tstr_8-1.d0
             F_v2u(i,j,k)= Hstag8( v2q(-1), v2q(0), v2q(1), v2q(2))
             F_t2v(i,j,k)= Hstag8(t2qu(-1),t2qu(0),t2qu(1),t2qu(2))/Cstv_Tstr_8-1.d0
             F_u2v(i,j,k)= Hstag8( u2q(-1), u2q(0), u2q(1), u2q(2))
             iJzq = (Ver_ext%m(k-1)+fis0(i,j)/grav_8) * VD3m2t(1,k)&
                   +(Ver_ext%m(k  )+fis0(i,j)/grav_8) * VD3m2t(2,k)&
                   +(Ver_ext%m(k+1)+fis0(i,j)/grav_8) * VD3m2t(3,k)&
                   +(Ver_ext%m(k+2)+fis0(i,j)/grav_8) * VD3m2t(4,k)
             F_dq2w(i,j,k) = dq2u(0) / iJzq
           !  F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
             dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
             dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
           !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k),dq2u(0)
         end do
         end do
      end do
      k=-1
      do j= 1, l_nj
         do i= 1, l_ni
            do n=-1,2
               dq2u(n) = ext_q(i+n,j,-1) * VD3m2t(1,k)&
                        +ext_q(i+n,j, 0) * VD3m2t(2,k)&
                        +ext_q(i+n,j, 1) * VD3m2t(3,k)&
                        +ext_q(i+n,j, 2) * VD3m2t(4,k)
               dq2v(n) = ext_q(i,j+n,-1) * VD3m2t(1,k)&
                        +ext_q(i,j+n, 0) * VD3m2t(2,k)&
                        +ext_q(i,j+n, 1) * VD3m2t(3,k)&
                        +ext_q(i,j+n, 2) * VD3m2t(4,k)
            end do
            iJzq = (Ver_ext%m(-1)+fis0(i,j)/grav_8) * VD3m2t(1,k)&
                  +(Ver_ext%m( 0)+fis0(i,j)/grav_8) * VD3m2t(2,k)&
                  +(Ver_ext%m( 1)+fis0(i,j)/grav_8) * VD3m2t(3,k)&
                  +(Ver_ext%m( 2)+fis0(i,j)/grav_8) * VD3m2t(4,k)
            F_dq2w(i,j,k) = dq2u(0) / iJzq
          !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k),dq2u(0)
          !         F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
           ! if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k)
            dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
            dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
          end do
      end do
      k=0
      do j= 1, l_nj
         do i= 1, l_ni
            do n=-1,2
               dq2u(n) = ext_q(i+n,j,k-1) * VD3m2t(1,k)&
                        +ext_q(i+n,j,k  ) * VD3m2t(2,k)&
                        +ext_q(i+n,j,k+1) * VD3m2t(3,k)&
                        +ext_q(i+n,j,k+2) * VD3m2t(4,k)
               dq2v(n) = ext_q(i,j+n,k-1) * VD3m2t(1,k)&
                        +ext_q(i,j+n,k  ) * VD3m2t(2,k)&
                        +ext_q(i,j+n,k+1) * VD3m2t(3,k)&
                        +ext_q(i,j+n,k+2) * VD3m2t(4,k)
            end do
            iJzq = (Ver_ext%m(-1)+fis0(i,j)/grav_8) * VD3m2t(1,k)&
                  +(Ver_ext%m( 0)+fis0(i,j)/grav_8) * VD3m2t(2,k)&
                  +(Ver_ext%m( 1)+fis0(i,j)/grav_8) * VD3m2t(3,k)&
                  +(Ver_ext%m( 2)+fis0(i,j)/grav_8) * VD3m2t(4,k)
            F_dq2w(i,j,k) = dq2u(0) / iJzq
          !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k),dq2u(0)
          !         F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
          !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k)
            !F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
            !if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo12: ',k,F_dq2w(i,j,k),M_iJzq(i,j,k)
            dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
            dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
          end do
      end do
      k=l_nk+1
      do j= 1, l_nj
         do i= 1, l_ni
            do n=-1,2
               dq2u(n) = ext_q(i+n,j,l_nk-1) * VD3m2t(1,k)&
                        +ext_q(i+n,j,l_nk  ) * VD3m2t(2,k)&
                        +ext_q(i+n,j,l_nk+1) * VD3m2t(3,k)&
                        +ext_q(i+n,j,l_nk+2) * VD3m2t(4,k)
               dq2v(n) = ext_q(i,j+n,l_nk-1) * VD3m2t(1,k)&
                        +ext_q(i,j+n,l_nk  ) * VD3m2t(2,k)&
                        +ext_q(i,j+n,l_nk+1) * VD3m2t(3,k)&
                        +ext_q(i,j+n,l_nk+2) * VD3m2t(4,k)
            end do
            iJzq = (Ver_ext%m(l_nk-1)+fis0(i,j)/grav_8) * VD3m2t(1,k)&
                  +(Ver_ext%m(l_nk  )+fis0(i,j)/grav_8) * VD3m2t(2,k)&
                  +(Ver_ext%m(l_nk+1)+fis0(i,j)/grav_8) * VD3m2t(3,k)&
                  +(Ver_ext%m(l_nk+2)+fis0(i,j)/grav_8) * VD3m2t(4,k)
            F_dq2w(i,j,k) = dq2u(0) / iJzq
          !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k),dq2u(0)
          !         F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
          !  if ((i==l_ni/2).and.(j==L_nj/2+1)) print*, 'allo266: ',k,F_dq2w(i,j,k)
          !  F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0)
            dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
            dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
          end do
      end do
      i=l_ni/2
      j=l_nj/2+1
      do k= -1, l_nk+1
!         print*, k,F_dq2w(i,j,k)
      enddo
      !---now interpolate back to momentum level---
      !Note: dqdzu and dqdzv are already staggered appropriatly from previous loop
      do k=1,l_nk
         do j= 1, l_nj 
            do i= 1, l_ni
               duu = dqzu(i,j,k-2) * VS3t2m(1,k)&
                   + dqzu(i,j,k-1) * VS3t2m(2,k)&
                   + dqzu(i,j,k  ) * VS3t2m(3,k)&
                   + dqzu(i,j,k+1) * VS3t2m(4,k)

               dvv = dqzv(i,j,k-2) * VS3t2m(1,k)&
                   + dqzv(i,j,k-1) * VS3t2m(2,k)&
                   + dqzv(i,j,k  ) * VS3t2m(3,k)&
                   + dqzv(i,j,k+1) * VS3t2m(4,k)

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
      deallocate (dqzv, dqzu, ext_t, ext_q)
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine prerhs3rdVH
