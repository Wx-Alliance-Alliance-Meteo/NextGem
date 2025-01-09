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

      subroutine prerhs3rdVH ( F_t2u  , F_v2u , F_t2v , F_u2v,&
                               F_dq2u , F_dq2v, F_dq2w       ,&
                               Minx, Maxx, Miny, Maxy, Nk )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use glb_ld
      use gmm_vt0
      use sol_mem
      use ver
      use metric
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), &
                    intent(OUT) :: F_t2u, F_v2u, F_t2v, F_u2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,Nk), &
                    intent(OUT) :: F_dq2u, F_dq2v
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,-1:Nk+1), intent(OUT) :: F_dq2w

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      real(kind=REAL64), dimension(-1:2) :: t2qv, t2qu, v2q, u2q
!
!     ---------------------------------------------------------------
!
      do k=1, l_nk
         do j= 1, l_nj
            do i= 1, l_ni
               do n=-1,2
                t2qv(n) = ext_t(i+n,j,k-2)*VS3t2m(1,k) &
                        + ext_t(i+n,j,k-1)*VS3t2m(2,k) &
                        + ext_t(i+n,j,k  )*VS3t2m(3,k) &
                        + ext_t(i+n,j,k+1)*VS3t2m(4,k)
                v2q(n)  = Hstag(vt0(i+n,j-2,k),vt0(i+n,j-1,k),&
                                vt0(i+n,j,k  ),vt0(i+n,j+1,k))
                t2qu(n) = ext_t(i,j+n,k-2)*VS3t2m(1,k) &
                        + ext_t(i,j+n,k-1)*VS3t2m(2,k) &
                        + ext_t(i,j+n,k  )*VS3t2m(3,k) &
                        + ext_t(i,j+n,k+1)*VS3t2m(4,k)
                u2q(n)  = Hstag(ut0(i-2,j+n,k),ut0(i-1,j+n,k),&
                                ut0(i,j+n,k  ),ut0(i+1,j+n,k))
               end do
               F_t2u(i,j,k)= Hstag8(t2qv(-1),t2qv(0),t2qv(1),t2qv(2))/Cstv_Tstr_8-1.d0
               F_v2u(i,j,k)= Hstag8( v2q(-1), v2q(0), v2q(1), v2q(2))
               F_t2v(i,j,k)= Hstag8(t2qu(-1),t2qu(0),t2qu(1),t2qu(2))/Cstv_Tstr_8-1.d0
               F_u2v(i,j,k)= Hstag8( u2q(-1), u2q(0), u2q(1), u2q(2))
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

      call dqdz3rd ( Sol_lhs, F_dq2u , F_dq2v, F_dq2w, Minx, Maxx, &
                     Miny, Maxy, Nk, -3, l_nk+4 )
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine prerhs3rdVH
