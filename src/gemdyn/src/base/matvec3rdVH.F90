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

      subroutine matvec3rdVH ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                             F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use geomh
      use dyn_fisl_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use glb_pil
      use ldnh
      use metric
      use omp_timing
      use sol_mem
      use mem_tstp
      use tdpack
      use ver
      use vgh
      use stat_mpi
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(IN ) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn        ,F_nk), intent(OUT) :: F_prod

      integer :: i, j, k, km1,km2,km3,kp1,kp2,kp3
      integer :: HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: r1(l_ni), r3(l_ni), s1(l_ni),s2(l_ni),s3(l_ni)
      real(kind=REAL64) :: w1(l_ni),w2(l_ni),w3(l_ni),w4(l_ni),&
                           w5(l_ni),w6(l_ni),w9(l_ni),c1,c2
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw, zero_8
      real(kind=REAL64), parameter :: half=0.5d0
!
!     ---------------------------------------------------------------
!
      call gtmg_start (91, 'MATVEC1', 29 )

      ext_q=0.
      ext_q(ds_i0:ds_in,ds_j0:ds_jn,1:l_nk)= F_vector(ds_i0:ds_in,ds_j0:ds_jn,1:l_nk)
      call fill_Vhalo (ext_q,l_minx,l_maxx,l_miny,l_maxy,lbound(ext_q,3),ubound(ext_q,3),1.d0)
      if ( .not. Grd_yinyang_L) then
         if (l_west) then
            do i=1,pil_w
               ext_q(i, 1+pil_s:l_nj-pil_s , :) = ext_q(1+pil_w, 1+pil_s:l_nj-pil_s , :)
            end do
         endif
         if (l_east) then
            do i=l_ni-pil_e+1,l_ni
               ext_q(i, 1+pil_s:l_nj-pil_s  , :) = ext_q(l_ni-pil_e, 1+pil_s:l_nj-pil_s , :)
            end do
         endif
         if (l_south) then
            do j=1,pil_s
               ext_q(1:l_ni , j , :) = ext_q(1:l_ni, 1+pil_s , :) 
            end do
         endif
         if (l_north) then
            do j=l_nj-pil_n+1,l_nj
               ext_q(1:l_ni, j , :) = ext_q(1:l_ni, l_nj-pil_n , :) 
            end do
         endif
      endif
      
      call delQ (ext_q,l_minx,l_maxx,l_miny,l_maxy, Qu,Qv,Qw,Qq,lbound(ext_q,3),ubound(ext_q,3))
      ext_q=0.
      ext_q(ds_i0:ds_in,ds_j0:ds_jn,1:l_nk)= Qw(ds_i0:ds_in,ds_j0:ds_jn,1:l_nk)
      call fill_Vhalo (ext_q,l_minx,l_maxx,l_miny,l_maxy,lbound(ext_q,3),ubound(ext_q,3),1.d0)
      
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

               barxQu = Hstag8(Qu(i-2,j,k), Qu(i-1,j,k),&
                               Qu(i  ,j,k), Qu(i+1,j,k) )

               baryQv = Hstag8(Qv(i,j-2,k), Qv(i,j-1,k),&
                               Qv(i,j  ,k), Qv(i,j+1,k) )

               barzQw =  ext_q(i,j,k-2) * VS3t2m(1,k) &
                      +  ext_q(i,j,k-1) * VS3t2m(2,k) &
                      +  ext_q(i,j,k  ) * VS3t2m(3,k) &
                      +  ext_q(i,j,k+1) * VS3t2m(4,k)

               !--- remains second order for now ---
               dzQw= Ver_idz_8%m(k)*(ext_q(i,j,k)-ext_q(i,j,k-1))
!!$               dzQw =  Qw(i,j,k-2) * VD3t2m(1,k) &
!!$                    +  Qw(i,j,k-1) * VD3t2m(2,k) &
!!$                    +  Qw(i,j,k  ) * VD3t2m(3,k) &
!!$                    +  Qw(i,j,k+1) * VD3t2m(4,k)
                    
               F_prod(i,j,k)= -gg_8*F_vector(i,j,k) + dxQu + dyQv + gama_8*dzQw &
                              +gama_8*barzQw*(M_logJzq(i,j,k)-epsi_8) &
                              +barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k)
            end do
         end do
      end do
!!$      do k= 1, l_nk
!!$         call statf_dm (F_prod(F_i0:,F_j0:,k:k), 'PROD', k, 'MATV', F_i0,F_in,F_j0,F_jn,1,1,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,1,8)
!!$      end do
!!$      call statf_dm (F_prod, 'PROD', k, 'MATV', F_i0,F_in,F_j0,F_jn,1,l_nk,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,l_nk,8)
      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine matvec3rdVH
