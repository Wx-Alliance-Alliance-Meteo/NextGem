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

      subroutine prerhs3rd ( F_t2u, F_v2u, F_t2v, F_u2v, &
               F_dq2u, F_dq2v, F_dq2w, Minx, Maxx, Miny, Maxy, Nk )
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
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), intent(OUT) :: F_dq2u, F_dq2v, F_dq2w

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      integer :: km1,km2,km3,kp1,kp2,kp3
      real(kind=REAL64), dimension(-1:2) :: t2qv, t2qu, v2q, u2q, dq2u, dq2v
      real(kind=REAL64) :: duu, dvv 
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) :: dqzv, dqzu
!
!     ---------------------------------------------------------------
!
      do k=1, l_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         kp3=min(k+3,G_nk)
         do j= 1, l_nj
           do i= 1, l_ni
             do n=-1,2
                t2qv(n) = tt0(i+n,j,km2)*CWt2m(1,k) + tt0(i+n,j,km1)*CWt2m(2,k) &
                        + tt0(i+n,j,k  )*CWt2m(3,k) + tt0(i+n,j,kp1)*CWt2m(4,k)

               v2q(n) = Hstag(vt0(i,j-2,k),vt0(i,j-1,k),&
                              vt0(i,j,k  ),vt0(i,j+1,k))
               t2qu(n) = tt0(i,j+n,km2)*CWt2m(1,k) + tt0(i,j+n,km1)*CWt2m(2,k) &
                       + tt0(i,j+n,k  )*CWt2m(3,k) + tt0(i,j+n,kp1)*CWt2m(4,k)

               u2q(n) = Hstag(ut0(i-2,j,k),ut0(i-1,j,k),&
                              ut0(i,j,k  ),ut0(i+1,j,k))
               dq2u(n) = qt0(i+n,j,km1) * CDm2t(1,k)&
                        +qt0(i+n,j,k  ) * CDm2t(2,k)&
                        +qt0(i+n,j,k+1) * CDm2t(3,k)&
                        +qt0(i+n,j,kp2) * CDm2t(4,k)

               dq2v(n) = qt0(i,j+n,km1) * CDm2t(1,k)&
                        +qt0(i,j+n,k  ) * CDm2t(2,k)&
                        +qt0(i,j+n,k+1) * CDm2t(3,k)&
                        +qt0(i,j+n,kp2) * CDm2t(4,k)
             end do
             F_t2u(i,j,k)= Hstag8(t2qv(-1),t2qv(0),t2qv(1),t2qv(2))/Cstv_Tstr_8-1.d0
             F_v2u(i,j,k)= Hstag8( v2q(-1), v2q(0), v2q(1), v2q(2))
             F_t2v(i,j,k)= Hstag8(t2qu(-1),t2qu(0),t2qu(1),t2qu(2))/Cstv_Tstr_8-1.d0
             F_u2v(i,j,k)= Hstag8( u2q(-1), u2q(0), u2q(1), u2q(2))
             F_dq2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0) !for Rw
             dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
             dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
         end do
         end do
      end do

      !---now interpolate back to momentum level---
      !Note: dqdzu and dqdzv are already staggered appropriatly from previous loop
      do k=1,l_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         do j= 1, l_nj 
            do i= 1, l_ni
               duu = dqzu(i,j,km2) * CWt2m(1,k)&
                   + dqzu(i,j,km1) * CWt2m(2,k)&
                   + dqzu(i,j,k  ) * CWt2m(3,k)&
                   + dqzu(i,j,kp1) * CWt2m(4,k)

               dvv = dqzv(i,j,km2) * CWt2m(1,k)&
                   + dqzv(i,j,km1) * CWt2m(2,k)&
                   + dqzv(i,j,k  ) * CWt2m(3,k)&
                   + dqzv(i,j,kp1) * CWt2m(4,k)

               F_dq2u(i,j,k)= M_Jxozu(i,j,k)! * duu 
               F_dq2v(i,j,k)= M_Jyozv(i,j,k)! * dvv
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
      include 'H3rd_ope.inc'
      end subroutine prerhs3rd
