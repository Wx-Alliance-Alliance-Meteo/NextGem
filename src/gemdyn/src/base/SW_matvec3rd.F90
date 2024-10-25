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

!** matvec - 3D Matrix-vector product

      subroutine SW_matvec3rd ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
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
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(IN ) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn        ,F_nk), intent(OUT) :: F_prod

      integer :: i, j, k, i0, j0, in, jn, k0, km, kp
      integer :: km1,km2,km3,kp1,kp2,kp3
!     real(kind=REAL64) :: r1(l_ni), r3(l_ni)
      real(kind=REAL64) :: r1(l_ni), r2(l_ni), r3(l_ni), r4(l_ni)
!     real(kind=REAL64) :: s1(l_ni),s2(l_ni),s3(l_ni)
      real(kind=REAL64) :: s1(l_ni),s2(l_ni),s3(l_ni),s4(l_ni),&
                           s5(l_ni),s6(l_ni),s7(l_ni),s8(l_ni)
!     real(kind=REAL64) :: w1(l_ni),w2(l_ni),w3(l_ni),w4(l_ni),&
!                          w5(l_ni),w6(l_ni),w9(l_ni)
      real(kind=REAL64) :: w1(l_ni),w2(l_ni),w3(l_ni),w4(l_ni),&
          w5(l_ni),w6(l_ni),w7(l_ni),w8(l_ni),w9(l_ni),w10(l_ni)
      real(kind=REAL64) :: dxQu, dyQv
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )

      do k= 1,l_nk
       do j= 1,l_nj
        do i= 1,l_ni
           ext_q(i,j,k)= F_vector(i,j,k)
        end do
       end do
      end do

      call SW_delQ3rd (ext_q, l_minx,l_maxx,l_miny,l_maxy,Qu,Qv,0,l_nk+1)

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
      do k= 1, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in

               dxQu = Hderiv8(Qu(i-2,j,k), Qu(i-1,j,k), &
                              Qu(i  ,j,k), Qu(i+1,j,k), geomh_invDXM_8(j))

               dyQv = Hderiv8(Qv(i,j-2,k)*geomh_cyM_8(j-2), &
                              Qv(i,j-1,k)*geomh_cyM_8(j-1), &
                              Qv(i,j  ,k)*geomh_cyM_8(j  ), &
                              Qv(i,j+1,k)*geomh_cyM_8(j+1), &
                              geomh_invDYM_8(j) )

               F_prod(i,j,k)= -gg_sw_8*ext_q(i,j,k) + dxQu + dyQv 

              !F_prod(i,j,k) = &
              !       -2.0d0*geomh_invDX_8(j)*geomh_invDXM_8(j)*ext_q(i,j,k)                 &
              !       +geomh_invDX_8(j)*geomh_invDXM_8(j)*ext_q(i+1,j,k)                     &
              !       +geomh_invDX_8(j)*geomh_invDXM_8(j)*ext_q(i-1,j,k)                     &
              !       -geomh_invDYMv_8(j)*geomh_cyM_8(j)*geomh_invDYM_8(j)*ext_q(i,j,k)      &
              !       -geomh_invDYMv_8(j-1)*geomh_cyM_8(j-1)*geomh_invDYM_8(j)*ext_q(i,j,k)  &
              !       +geomh_invDYMv_8(j)*geomh_cyM_8(j)*geomh_invDYM_8(j)*ext_q(i,j+1,k)    &
              !       +geomh_invDYMv_8(j-1)*geomh_cyM_8(j-1)*geomh_invDYM_8(j)*ext_q(i,j-1,k)&
              !       - gg_sw_8*ext_q(i,j,k)
            end do
         end do
      end do
      call gtmg_stop (92)

!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine SW_matvec3rd
