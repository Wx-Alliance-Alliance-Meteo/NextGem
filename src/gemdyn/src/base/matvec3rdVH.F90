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
      use glb_pil
      use glb_ld
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

      integer :: i, j, k, kk
      integer :: HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: dxQu, dyQv, barxQu, baryQv, barzQw, dzQw, zero_8, d1,d2,d3,d4,b1
      real(kind=REAL64), parameter :: half=0.5d0
      real(kind=REAL64), dimension(:,:,:), allocatable :: barz,delz
!
!     ---------------------------------------------------------------
!
!      print*, 'F_prod matvec3rdVH'
      call gtmg_start (91, 'MATVEC1', 29 )
      allocate (barz(l_ni,l_nj,-1:l_nk+1),delz(l_ni,l_nj,-1:l_nk+1))

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

      do k= -1, l_nk+2
     !    print*, 'QT0: ',k,ext_q(l_ni/2,l_nj/2+1,k  )
      end do
      
      do k= 0, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
                        b1=  ext_q(i,j,k-1) * VS3m2t(1,k) & 
                            +ext_q(i,j,k  ) * VS3m2t(2,k) & 
                            +ext_q(i,j,k+1) * VS3m2t(3,k) & 
                            +ext_q(i,j,k+2) * VS3m2t(4,k)                      
               delz(i,j,k)=  ext_q(i,j,k-1) * VD3m2t(1,k) & 
                            +ext_q(i,j,k  ) * VD3m2t(2,k) & 
                            +ext_q(i,j,k+1) * VD3m2t(3,k) & 
                            +ext_q(i,j,k+2) * VD3m2t(4,k) - mu_8*b1
            end do
         end do
      end do
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
                      b1= ext_q(i,j,-1) * VS3m2t(1,-1) & 
                         +ext_q(i,j, 0) * VS3m2t(2,-1) & 
                         +ext_q(i,j, 1) * VS3m2t(3,-1) & 
                         +ext_q(i,j, 2) * VS3m2t(4,-1)                      
            delz(i,j,-1)= ext_q(i,j,-1) * VD3m2t(1,-1) & 
                         +ext_q(i,j, 0) * VD3m2t(2,-1) & 
                         +ext_q(i,j, 1) * VD3m2t(3,-1) & 
                         +ext_q(i,j, 2) * VD3m2t(4,-1) - mu_8*b1
                      b1= ext_q(i,j,l_nk-1) * VS3m2t(1,l_nk+1) & 
                         +ext_q(i,j,l_nk  ) * VS3m2t(2,l_nk+1) & 
                         +ext_q(i,j,l_nk+1) * VS3m2t(3,l_nk+1) & 
                         +ext_q(i,j,l_nk+2) * VS3m2t(4,l_nk+1)                      
            delz(i,j,l_nk+1)= ext_q(i,j,l_nk-1) * VD3m2t(1,l_nk+1) & 
                             +ext_q(i,j,l_nk  ) * VD3m2t(2,l_nk+1) & 
                             +ext_q(i,j,l_nk+1) * VD3m2t(3,l_nk+1) & 
                             +ext_q(i,j,l_nk+2) * VD3m2t(4,l_nk+1) - mu_8*b1
         end do
      end do
      i=l_ni/2
      j=l_nj/2+1
      do k= -1, l_nk+1
     !    print*, k,delz(i,j,k)
      enddo

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

               dzQw  =  VD3t2m(1,k)*delz(i,j,k-2)+VD3t2m(2,k)*delz(i,j,k-1)+VD3t2m(3,k)*delz(i,j,k)+VD3t2m(4,k)*delz(i,j,k+1)
               barzQw=  VS3t2m(1,k)*delz(i,j,k-2)+VS3t2m(2,k)*delz(i,j,k-1)+VS3t2m(3,k)*delz(i,j,k)+VS3t2m(4,k)*delz(i,j,k+1)
                    
               F_prod(i,j,k)= -gg_8*ext_q(i,j,k) + dxQu + dyQv + gama_8*dzQw &
                              +gama_8*barzQw*(M_logJzq(i,j,k)-epsi_8) &
                              +barxQu*M_logJzu(i,j,k) + baryQv*M_logJzv(i,j,k)
!if ((i==l_ni/2).and.(j==l_nj/2+1)) write(6,'(a,i3,5(1pe22.12))') 'SOL_LHS: ',k, F_prod(i,j,k)!,-gg_8*ext_q(i,j,k) ,gama_8*dzQw &
                              !,gama_8*barzQw*(M_logJzq(i,j,k)-epsi_8)
                           end do
         end do
      end do
      deallocate (barz,delz)

      call gtmg_stop (92)
!     
!     ---------------------------------------------------------------
!     
      return
      include 'H3rd_ope.inc'
      end subroutine matvec3rdVH
