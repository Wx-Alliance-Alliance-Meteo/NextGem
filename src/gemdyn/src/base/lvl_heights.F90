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

!**s/r lvl_heights - Compute level heights

      subroutine lvl_heights (zmom_8,ztht_8, F_topo, F_orols, Minx,Maxx,Miny,Maxy)
      use, intrinsic :: iso_fortran_env
      use gem_options
      use tdpack
      use glb_ld
      use ver
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_orols
      real(kind=REAL64), dimension(Minx:Maxx,Miny:Maxy,0:G_nk+1), intent(OUT) :: zmom_8,ztht_8

      integer :: i,j,k
!
!     ---------------------------------------------------------------
!
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               zmom_8(i,j,k)=ver_z_8%m(k)+(Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_orols(i,j))/grav_8
               ztht_8(i,j,k)=ver_z_8%t(k)+(Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_orols(i,j))/grav_8
            end do
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            ztht_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,0)=ver_z_8%m(0)
            zmom_8(i,j,G_nk+1)= F_topo(i,j)/grav_8
            ztht_8(i,j,G_nk+1)= zmom_8(i,j,G_nk+1)
         end do
      end do
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine lvl_heights
      
      subroutine lvl_heights2 (F_metric, F_m2t, F_t2m, F_topo, F_orols, Minx,Maxx,Miny,Maxy)
      use, intrinsic :: iso_fortran_env
      use gem_options
      use metric
      use tdpack
      use glb_ld
      use ver
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy
      real, dimension (Minx:Maxx,Miny:Maxy), intent(IN) :: F_topo, F_orols
      type(Vmetric), intent(INOUT) :: F_metric
      type(Vops)   , intent(INOUT) :: F_m2t, F_t2m

      integer :: i,j,k,n
      real(kind=REAL64) :: p(6),z(6),s(6),w1
!     
!     ---------------------------------------------------------------
!
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               F_metric%zmom_8(i,j,k)=ver_z_8%m(k)+(Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_orols(i,j))/grav_8
               F_metric%ztht_8(i,j,k)=ver_z_8%t(k)+(Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_orols(i,j))/grav_8
            end do
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            F_metric%ztht_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,0)=ver_z_8%m(0)
            F_metric%zmom_8(i,j,G_nk+1)= F_topo(i,j)/grav_8
            F_metric%ztht_8(i,j,G_nk+1)= F_metric%zmom_8(i,j,G_nk+1)
         end do
      end do
!!$omp enddo
      do k=3,G_nk-3
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               do n=1,6
                  p(n)= F_metric%ztht_8(i,j,k)-F_metric%zmom_8(i,j,k+n-3)
               end do
               s(1) = p(2)*p(3)*p(4)*p(5)*p(6)
               s(2) = p(1)*p(3)*p(4)*p(5)*p(6)
               s(3) = p(1)*p(2)*p(4)*p(5)*p(6)
               s(4) = p(1)*p(2)*p(3)*p(5)*p(6)
               s(5) = p(1)*p(2)*p(3)*p(4)*p(6)
               s(6) = p(1)*p(2)*p(3)*p(4)*p(5)
               w1= F_metric%zmom_8(i,j,k-2)
               z(1)= (w1-F_metric%zmom_8(i,j,k-1))&
                    *(w1-F_metric%zmom_8(i,j,k  ))&
                    *(w1-F_metric%zmom_8(i,j,k+1))&
                    *(w1-F_metric%zmom_8(i,j,k+2))&
                    *(w1-F_metric%zmom_8(i,j,k+3))
               w1= F_metric%zmom_8(i,j,k-1)
               z(2)= (w1-F_metric%zmom_8(i,j,k-2))&
                    *(w1-F_metric%zmom_8(i,j,k  ))&
                    *(w1-F_metric%zmom_8(i,j,k+1))&
                    *(w1-F_metric%zmom_8(i,j,k+2))&
                    *(w1-F_metric%zmom_8(i,j,k+3))
               w1= F_metric%zmom_8(i,j,k)
               z(3)= (w1-F_metric%zmom_8(i,j,k-2))&
                    *(w1-F_metric%zmom_8(i,j,k-1))&
                    *(w1-F_metric%zmom_8(i,j,k+1))&
                    *(w1-F_metric%zmom_8(i,j,k+2))&
                    *(w1-F_metric%zmom_8(i,j,k+3))
               w1= F_metric%zmom_8(i,j,k+1)
               z(4)= (w1-F_metric%zmom_8(i,j,k-2))&
                    *(w1-F_metric%zmom_8(i,j,k-1))&
                    *(w1-F_metric%zmom_8(i,j,k ))&
                    *(w1-F_metric%zmom_8(i,j,k+2))&
                    *(w1-F_metric%zmom_8(i,j,k+3))
               w1= F_metric%zmom_8(i,j,k+2)
               z(5)= (w1-F_metric%zmom_8(i,j,k-2))&
                    *(w1-F_metric%zmom_8(i,j,k-1))&
                    *(w1-F_metric%zmom_8(i,j,k  ))&
                    *(w1-F_metric%zmom_8(i,j,k+1))&
                    *(w1-F_metric%zmom_8(i,j,k+3))
               w1= F_metric%zmom_8(i,j,k+3)
               z(6)= (w1-F_metric%zmom_8(i,j,k-2))&
                    *(w1-F_metric%zmom_8(i,j,k-1))&
                    *(w1-F_metric%zmom_8(i,j,k  ))&
                    *(w1-F_metric%zmom_8(i,j,k+1))&
                    *(w1-F_metric%zmom_8(i,j,k+2))
               F_m2t%L1(i,j,k) = s(1) / z(1)
               F_m2t%L2(i,j,k) = s(2) / z(2)
               F_m2t%L3(i,j,k) = s(3) / z(3)
               F_m2t%L4(i,j,k) = s(4) / z(4)
               F_m2t%L5(i,j,k) = s(5) / z(5)
               F_m2t%L6(i,j,k) = s(6) / z(6)
            end do
         end do
      end do

      F_m2t%L1(:,:,G_nk)= 0.d0;F_m2t%L2(:,:,G_nk)= 0.d0;F_m2t%L3(:,:,G_nk)= 0.d0
      F_m2t%L4(:,:,G_nk)= 0.d0;F_m2t%L5(:,:,G_nk)= 0.d0;F_m2t%L6(:,:,G_nk)= 0.d0
      F_m2t%L1(:,:,       1:2)= 0.d0 ; F_m2t%L6(:,:,       1:2)= 0.d0
      F_m2t%L1(:,:,  G_nk-2: )= 0.d0 ; F_m2t%L6(:,:,  G_nk-2: )= 0.d0
      F_m2t%L2(:,:,   1)= 0.d0 ; F_m2t%L5(:,:,   1)= 0.d0
      F_m2t%L2(:,:,G_nk-1)= 0.d0 ; F_m2t%L5(:,:,G_nk-1)= 0.d0
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            do n=2,5
               p(n)= F_metric%ztht_8(i,j,2)-F_metric%zmom_8(i,j,n-1)
            end do
            s(2) = p(3)*p(4)*p(5)
            s(3) = p(2)*p(4)*p(5)
            s(4) = p(2)*p(3)*p(5)
            s(5) = p(2)*p(3)*p(4)
            w1= F_metric%zmom_8(i,j,1)
            z(2)= (w1-F_metric%zmom_8(i,j,2))&
                 *(w1-F_metric%zmom_8(i,j,3))&
                 *(w1-F_metric%zmom_8(i,j,4))
            w1= F_metric%zmom_8(i,j,2)
            z(3)= (w1-F_metric%zmom_8(i,j,1))&
                 *(w1-F_metric%zmom_8(i,j,3))&
                 *(w1-F_metric%zmom_8(i,j,4))
            w1= F_metric%zmom_8(i,j,3)
            z(4)= (w1-F_metric%zmom_8(i,j,1))&
                 *(w1-F_metric%zmom_8(i,j,2))&
                 *(w1-F_metric%zmom_8(i,j,4))
            w1= F_metric%zmom_8(i,j,4)
            z(5)= (w1-F_metric%zmom_8(i,j,1))&
                 *(w1-F_metric%zmom_8(i,j,2))&
                 *(w1-F_metric%zmom_8(i,j,3))
            F_m2t%L2(i,j,2) = s(2) / z(2)
            F_m2t%L3(i,j,2) = s(3) / z(3)
            F_m2t%L4(i,j,2) = s(4) / z(4)
            F_m2t%L5(i,j,2) = s(5) / z(5)
            F_m2t%L3(i,j,1) = (F_metric%zmom_8(i,j,2)-F_metric%ztht_8(i,j,1))&
                             /(F_metric%zmom_8(i,j,2)-F_metric%zmom_8(i,j,1))
            F_m2t%L4(i,j,1) = 1.d0 - F_m2t%L3(i,j,1)
            do n=2,5
               p(n)= F_metric%ztht_8(i,j,G_nk-2)-F_metric%zmom_8(i,j,G_nk+n-5)
            end do
            s(2) = p(3)*p(4)*p(5)
            s(3) = p(2)*p(4)*p(5)
            s(4) = p(2)*p(3)*p(5)
            s(5) = p(2)*p(3)*p(4)
            w1= F_metric%zmom_8(i,j,G_nk-3)
            z(2)= (w1-F_metric%zmom_8(i,j,G_nk-2))&
                 *(w1-F_metric%zmom_8(i,j,G_nk-1))&
                 *(w1-F_metric%zmom_8(i,j,G_nk  ))
            w1= F_metric%zmom_8(i,j,G_nk-2)
            z(3)= (w1-F_metric%zmom_8(i,j,G_nk-3))&
                 *(w1-F_metric%zmom_8(i,j,G_nk-1))&
                 *(w1-F_metric%zmom_8(i,j,G_nk  ))
            w1= F_metric%zmom_8(i,j,G_nk-1)
            z(4)= (w1-F_metric%zmom_8(i,j,G_nk-3))&
                 *(w1-F_metric%zmom_8(i,j,G_nk-2))&
                 *(w1-F_metric%zmom_8(i,j,G_nk  ))
            w1= F_metric%zmom_8(i,j,G_nk  )
            z(5)= (w1-F_metric%zmom_8(i,j,G_nk-3))&
                 *(w1-F_metric%zmom_8(i,j,G_nk-2))&
                 *(w1-F_metric%zmom_8(i,j,G_nk-1))
            F_m2t%L2(i,j,G_nk-2) = s(2) / z(2)
            F_m2t%L3(i,j,G_nk-2) = s(3) / z(3)
            F_m2t%L4(i,j,G_nk-2) = s(4) / z(4)
            F_m2t%L5(i,j,G_nk-2) = s(5) / z(5)
            F_m2t%L3(i,j,G_nk-1) = (F_metric%zmom_8(i,j,G_nk)-F_metric%ztht_8(i,j,G_nk-1))&
                                  /(F_metric%zmom_8(i,j,G_nk)-F_metric%zmom_8(i,j,G_nk-1))
            F_m2t%L4(i,j,G_nk-1) = 1.d0 - F_m2t%L3(i,j,G_nk-1)
         end do
      end do

      do k=4,G_nk-2
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               do n=1,6
                  p(n)= F_metric%zmom_8(i,j,k)-F_metric%ztht_8(i,j,k+n-4)
               end do
               s(1) = p(2)*p(3)*p(4)*p(5)*p(6)
               s(2) = p(1)*p(3)*p(4)*p(5)*p(6)
               s(3) = p(1)*p(2)*p(4)*p(5)*p(6)
               s(4) = p(1)*p(2)*p(3)*p(5)*p(6)
               s(5) = p(1)*p(2)*p(3)*p(4)*p(6)
               s(6) = p(1)*p(2)*p(3)*p(4)*p(5)
               w1= F_metric%ztht_8(i,j,k-3)
               z(1)= (w1-F_metric%ztht_8(i,j,k-2))&
                    *(w1-F_metric%ztht_8(i,j,k-1))&
                    *(w1-F_metric%ztht_8(i,j,k ))&
                    *(w1-F_metric%ztht_8(i,j,k+1))&
                    *(w1-F_metric%ztht_8(i,j,k+2))
               w1= F_metric%ztht_8(i,j,k-2)
               z(2)= (w1-F_metric%ztht_8(i,j,k-3))&
                    *(w1-F_metric%ztht_8(i,j,k-1))&
                    *(w1-F_metric%ztht_8(i,j,k  ))&
                    *(w1-F_metric%ztht_8(i,j,k+1))&
                    *(w1-F_metric%ztht_8(i,j,k+2))
               w1= F_metric%ztht_8(i,j,k-1)
               z(3)= (w1-F_metric%ztht_8(i,j,k-3))&
                    *(w1-F_metric%ztht_8(i,j,k-2))&
                    *(w1-F_metric%ztht_8(i,j,k  ))&
                    *(w1-F_metric%ztht_8(i,j,k+1))&
                    *(w1-F_metric%ztht_8(i,j,k+2))
               w1= F_metric%ztht_8(i,j,k  )
               z(4)= (w1-F_metric%ztht_8(i,j,k-3))&
                    *(w1-F_metric%ztht_8(i,j,k-2))&
                    *(w1-F_metric%ztht_8(i,j,k-1))&
                    *(w1-F_metric%ztht_8(i,j,k+1))&
                    *(w1-F_metric%ztht_8(i,j,k+2))
               w1= F_metric%ztht_8(i,j,k+1)
               z(5)= (w1-F_metric%ztht_8(i,j,k-3))&
                    *(w1-F_metric%ztht_8(i,j,k-2))&
                    *(w1-F_metric%ztht_8(i,j,k-1))&
                    *(w1-F_metric%ztht_8(i,j,k ))&
                    *(w1-F_metric%ztht_8(i,j,k+2))
               w1= F_metric%ztht_8(i,j,k+2)
               z(6)= (w1-F_metric%ztht_8(i,j,k-3))&
                    *(w1-F_metric%ztht_8(i,j,k-2))&
                    *(w1-F_metric%ztht_8(i,j,k-1))&
                    *(w1-F_metric%ztht_8(i,j,k  ))&
                    *(w1-F_metric%ztht_8(i,j,k+1))
               F_t2m%L1(i,j,k) = s(1) / z(1)
               F_t2m%L2(i,j,k) = s(2) / z(2)
               F_t2m%L3(i,j,k) = s(3) / z(3)
               F_t2m%L4(i,j,k) = s(4) / z(4)
               F_t2m%L5(i,j,k) = s(5) / z(5)
               F_t2m%L6(i,j,k) = s(6) / z(6)
            end do
         end do
      end do
      
      F_t2m%L1(:,:,1)= 0.d0;F_t2m%L2(:,:,1)= 0.d0;F_t2m%L3(:,:,1)= 0.d0
      F_t2m%L4(:,:,1)= 0.d0;F_t2m%L5(:,:,1)= 0.d0;F_t2m%L6(:,:,1)= 0.d0
      F_t2m%L1(:,:,       2:3)= 0.d0 ; F_t2m%L6(:,:,       2:3)= 0.d0
      F_t2m%L1(:,:,  G_nk-1: )= 0.d0 ; F_t2m%L6(:,:,  G_nk-1: )= 0.d0
      F_t2m%L2(:,:,   2)= 0.d0 ; F_t2m%L5(:,:,   2)= 0.d0
      F_t2m%L2(:,:,G_nk)= 0.d0 ; F_t2m%L5(:,:,G_nk)= 0.d0
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            do n=2,5
               p(n)= F_metric%zmom_8(i,j,3)-F_metric%ztht_8(i,j,n-1)
            end do
            s(2) = p(3)*p(4)*p(5)
            s(3) = p(2)*p(4)*p(5)
            s(4) = p(2)*p(3)*p(5)
            s(5) = p(2)*p(3)*p(4)
            w1= F_metric%ztht_8(i,j,1)
            z(2)= (w1-F_metric%ztht_8(i,j,2))&
                 *(w1-F_metric%ztht_8(i,j,3))&
                 *(w1-F_metric%ztht_8(i,j,4))
            w1= F_metric%ztht_8(i,j,2)
            z(3)= (w1-F_metric%ztht_8(i,j,1))&
                 *(w1-F_metric%ztht_8(i,j,3))&
                 *(w1-F_metric%ztht_8(i,j,4))
            w1= F_metric%ztht_8(i,j,3)
            z(4)= (w1-F_metric%ztht_8(i,j,1))&
                 *(w1-F_metric%ztht_8(i,j,2))&
                 *(w1-F_metric%ztht_8(i,j,4))
            w1= F_metric%ztht_8(i,j,4)
            z(5)= (w1-F_metric%ztht_8(i,j,1))&
                 *(w1-F_metric%ztht_8(i,j,2))&
                 *(w1-F_metric%ztht_8(i,j,3))
            F_t2m%L2(i,j,3) = s(2) / z(2)
            F_t2m%L3(i,j,3) = s(3) / z(3)
            F_t2m%L4(i,j,3) = s(4) / z(4)
            F_t2m%L5(i,j,3) = s(5) / z(5)
            F_t2m%L3(i,j,2) = (F_metric%ztht_8(i,j,2)-F_metric%zmom_8(i,j,2))&
                             /(F_metric%ztht_8(i,j,2)-F_metric%ztht_8(i,j,1))
            F_t2m%L4(i,j,2) = 1.d0 - F_t2m%L3(i,j,2)
            do n=2,5
               p(n)= F_metric%zmom_8(i,j,G_nk-1)-F_metric%ztht_8(i,j,G_nk+n-5)
            end do
            s(2) = p(3)*p(4)*p(5)
            s(3) = p(2)*p(4)*p(5)
            s(4) = p(2)*p(3)*p(5)
            s(5) = p(2)*p(3)*p(4)
            w1= F_metric%ztht_8(i,j,G_nk-3)
            z(2)= (w1-F_metric%ztht_8(i,j,G_nk-2))&
                 *(w1-F_metric%ztht_8(i,j,G_nk-1))&
                 *(w1-F_metric%ztht_8(i,j,G_nk  ))
            w1= F_metric%ztht_8(i,j,G_nk-2)
            z(3)= (w1-F_metric%ztht_8(i,j,G_nk-3))&
                 *(w1-F_metric%ztht_8(i,j,G_nk-1))&
                 *(w1-F_metric%ztht_8(i,j,G_nk  ))
            w1= F_metric%ztht_8(i,j,G_nk-1)
            z(4)= (w1-F_metric%ztht_8(i,j,G_nk-3))&
                 *(w1-F_metric%ztht_8(i,j,G_nk-2))&
                 *(w1-F_metric%ztht_8(i,j,G_nk  ))
            w1= F_metric%ztht_8(i,j,G_nk  )
            z(5)= (w1-F_metric%ztht_8(i,j,G_nk-3))&
                 *(w1-F_metric%ztht_8(i,j,G_nk-2))&
                 *(w1-F_metric%ztht_8(i,j,G_nk-1))
            F_t2m%L2(i,j,G_nk-1) = s(2) / z(2)
            F_t2m%L3(i,j,G_nk-1) = s(3) / z(3)
            F_t2m%L4(i,j,G_nk-1) = s(4) / z(4)
            F_t2m%L5(i,j,G_nk-1) = s(5) / z(5)
            F_t2m%L3(i,j,G_nk) = (F_metric%ztht_8(i,j,G_nk)-F_metric%zmom_8(i,j,G_nk  ))&
                                /(F_metric%ztht_8(i,j,G_nk)-F_metric%ztht_8(i,j,G_nk-1))
            F_t2m%L4(i,j,G_nk) = 1.d0 - F_t2m%L3(i,j,G_nk)            
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine lvl_heights2
