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

!** matvec5 - 5th order 3D Matrix-vector product

      subroutine matvec5 ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                           F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use ldnh
      use metric
      use omp_timing
      use sol_mem
      use mem_tstp
      use tdpack
      use ver
      use vgh
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(IN ) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn        ,F_nk), intent(OUT) :: F_prod

      integer :: i, j, k
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )

      do k= 1, l_nk
         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            vgh_q(i,j,k)= F_vector(i,j,k)
         end do
         end do
      end do

      call delQ5 (vgh_q, l_minx,l_maxx,l_miny,l_maxy, Qu,Qv,Qw,Qq,-2,l_nk+3)
      call fill_Vhalo (Qw,l_minx,l_maxx,l_miny,l_maxy,1.d0)

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
      do k= 1, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in

               dxQu = Hderiv8(Qu(i-3,j,k), Qu(i-2,j,k), Qu(i-1,j,k), &
                              Qu(i  ,j,k), Qu(i+1,j,k), Qu(i+2,j,k), geomh_invDXM_8(j))

               dyQv = Hderiv8(Qv(i,j-3,k)*geomh_cyM_8(j-3), &
                              Qv(i,j-2,k)*geomh_cyM_8(j-2), &
                              Qv(i,j-1,k)*geomh_cyM_8(j-1), &
                              Qv(i,j  ,k)*geomh_cyM_8(j  ), &
                              Qv(i,j+1,k)*geomh_cyM_8(j+1), &
                              Qv(i,j+2,k)*geomh_cyM_8(j+2), &
                              geomh_invDYM_8(j) )

               barxQu = Hstag8(Qu(i-3,j,k), Qu(i-2,j,k), Qu(i-1,j,k), &
                               Qu(i  ,j,k), Qu(i+1,j,k), Qu(i+2,j,k))

               baryQv = Hstag8(Qv(i,j-3,k), Qv(i,j-2,k), Qv(i,j-1,k), &
                               Qv(i,j  ,k), Qv(i,j+1,k), Qv(i,j+2,k))

               barzQw =  Qw(i,j,k-3) * VS5t2m(1,k) & 
                      +  Qw(i,j,k-2) * VS5t2m(2,k) &
                      +  Qw(i,j,k-1) * VS5t2m(3,k) &
                      +  Qw(i,j,k  ) * VS5t2m(4,k) &
                      +  Qw(i,j,k+1) * VS5t2m(5,k) &
                      +  Qw(i,j,k+2) * VS5t2m(6,k) 

               dzQw =  Qw(i,j,k-3) * VD5t2m(1,k) & 
                    +  Qw(i,j,k-2) * VD5t2m(2,k) &
                    +  Qw(i,j,k-1) * VD5t2m(3,k) &
                    +  Qw(i,j,k  ) * VD5t2m(4,k) &
                    +  Qw(i,j,k+1) * VD5t2m(5,k) &
                    +  Qw(i,j,k+2) * VD5t2m(6,k) 

               F_prod(i,j,k)= -gg_8*vgh_q(i,j,k) + dxQu + dyQv + gama_8*dzQw &
                              +gama_8*barzQw*(M_logJzq(i,j,k)-epsi_8) &
                              +barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k)
            end do
         end do
      end do
      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H5th_ope.inc'
      end subroutine matvec5
