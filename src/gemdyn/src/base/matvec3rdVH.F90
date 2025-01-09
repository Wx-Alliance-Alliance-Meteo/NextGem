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

      subroutine matvec3rdVH ( F_prod,F_i0,F_in,F_j0,F_jn, F_nk, F_flag )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_pil
      use glb_ld
      use ldnh
      use metric
      use omp_timing
      use gmm_vt0
      use sol_mem
      use mem_tstp
      use tdpack
      use ver
      use vgh
      use stat_mpi
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_i0,F_in,F_j0,F_jn,F_nk,F_flag
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(OUT) :: F_prod

      integer :: i, j, k, kk, ub, dim
      integer :: HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw,b1,q
      real(kind=REAL64), parameter :: half=0.5d0, one=1.0d0
      real(kind=REAL64), dimension(:,:,:), pointer :: dqz2u, dqz2v, dqz2w
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )
      matvq => Sol_lhs
      if ( F_flag > 0 ) then
      matvq => kryq
      do k= 0,-3,-1
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            q= (VM3%zmom(i,j,k+1)-VM3%zmom(i,j,k))/(VM3%zmom(i,j,k+2)-VM3%zmom(i,j,k+1))
            matvq(i,j,k)= (1+q)*matvq(i,j,k+1) - q*matvq(i,j,k+2)
      end do
      end do
      end do
      do k= l_nk+1, l_nk+4
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            q= (VM3%zmom(i,j,k)-VM3%zmom(i,j,k-1))/(VM3%zmom(i,j,k-1)-VM3%zmom(i,j,k-2))
            matvq(i,j,k)= (1+q)*matvq(i,j,k-1) - q*matvq(i,j,k-2)
      end do
      end do
      end do
      if ( .not. Grd_yinyang_L) then
         if (l_west) then
            do i=1,pil_w
               matvq(i, 1+pil_s:l_nj-pil_s , :) = matvq(1+pil_w, 1+pil_s:l_nj-pil_s , :)
            end do
         endif
         if (l_east) then
            do i=l_ni-pil_e+1,l_ni
               matvq(i, 1+pil_s:l_nj-pil_s  , :) = matvq(l_ni-pil_e, 1+pil_s:l_nj-pil_s , :)
            end do
         endif
         if (l_south) then
            do j=1,pil_s
               matvq(1:l_ni , j , :) = matvq(1:l_ni, 1+pil_s , :) 
            end do
         endif
         if (l_north) then
            do j=l_nj-pil_n+1,l_nj
               matvq(1:l_ni, j , :) = matvq(1:l_ni, l_nj-pil_n , :) 
            end do
         endif
      endif
      endif

      !---vertical derivative of q to points u,v,w in rhs---
      ub=0
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      dqz2u (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2v (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2w (l_minx:l_maxx, l_miny:l_maxy, -1:l_nk+1) => WS1_8(ub+1:); ub=ub+dim*(l_nk+3)

      call dqdz3rd ( matvq, dqz2u , dqz2v, dqz2w, l_minx,l_maxx,&
                     l_miny,l_maxy, G_nk, -3, l_nk+4 )
      
      do k= 0, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               b1=  matvq(i,j,k-1) * VS3m2t(1,k) & 
                   +matvq(i,j,k  ) * VS3m2t(2,k) & 
                   +matvq(i,j,k+1) * VS3m2t(3,k) & 
                   +matvq(i,j,k+2) * VS3m2t(4,k)
               dqz2w(i,j,k) = dqz2w(i,j,k) - mu_8*b1
            end do
         end do
      end do
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            b1= matvq(i,j,-1) * VS3m2t(1,-1) & 
               +matvq(i,j, 0) * VS3m2t(2,-1) & 
               +matvq(i,j, 1) * VS3m2t(3,-1) & 
               +matvq(i,j, 2) * VS3m2t(4,-1)                      
            dqz2w(i,j,-1) = dqz2w(i,j,-1) - mu_8*b1
            b1= matvq(i,j,l_nk-1) * VS3m2t(1,l_nk+1) & 
               +matvq(i,j,l_nk  ) * VS3m2t(2,l_nk+1) & 
               +matvq(i,j,l_nk+1) * VS3m2t(3,l_nk+1) & 
               +matvq(i,j,l_nk+2) * VS3m2t(4,l_nk+1)
            dqz2w(i,j,l_nk+1) = dqz2w(i,j,l_nk+1) - mu_8*b1
         end do
      end do

      call gtmg_stop (91)
      call gtmg_start (92, 'MATVEC2', 29 )
       
      do k= 1, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in

               dxQu = Hderiv8(dqz2u(i-2,j,k), dqz2u(i-1,j,k), &
                              dqz2u(i  ,j,k), dqz2u(i+1,j,k), geomh_invDXM_8(j))

               dyQv = Hderiv8(dqz2v(i,j-2,k)*geomh_cyM_8(j-2), &
                              dqz2v(i,j-1,k)*geomh_cyM_8(j-1), &
                              dqz2v(i,j  ,k)*geomh_cyM_8(j  ), &
                              dqz2v(i,j+1,k)*geomh_cyM_8(j+1), &
                              geomh_invDYM_8(j) )

               barxQu = Hstag8(dqz2u(i-2,j,k), dqz2u(i-1,j,k),&
                               dqz2u(i  ,j,k), dqz2u(i+1,j,k) )

               baryQv = Hstag8(dqz2v(i,j-2,k), dqz2v(i,j-1,k),&
                               dqz2v(i,j  ,k), dqz2v(i,j+1,k) )

               dzQw  =  VD3t2m(1,k)*dqz2w(i,j,k-2)+VD3t2m(2,k)*dqz2w(i,j,k-1)+VD3t2m(3,k)*dqz2w(i,j,k)+VD3t2m(4,k)*dqz2w(i,j,k+1)
               barzQw=  VS3t2m(1,k)*dqz2w(i,j,k-2)+VS3t2m(2,k)*dqz2w(i,j,k-1)+VS3t2m(3,k)*dqz2w(i,j,k)+VS3t2m(4,k)*dqz2w(i,j,k+1)
                    
               F_prod(i,j,k)= -gg_8*matvq(i,j,k) + dxQu + dyQv + gama_8*dzQw &
                              -gama_8*barzQw*epsi_8 &
                              +barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k)&
                              +gama_8*barzQw*M_logJzq(i,j,k)
            end do
         end do
      end do

      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine matvec3rdVH
