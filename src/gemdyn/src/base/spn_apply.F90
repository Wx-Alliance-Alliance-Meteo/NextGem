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
!---------------------------------- LICENCE END --------------------------------

!*s/r spn_fld - 2D forward FFT + filter + backward FFT 
!               and apply nudging tendency

      subroutine spn_apply ( F_ft1, F_nest, Minx, Maxx, Miny, Maxy, Nk )
      use spn_options
      use glb_ld
      use glb_pil
      use HORgrid_options
      use ldnh
      use gem_fft
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      !! Locally-computed field, to be modified in-place with the nudging
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(INOUT) :: F_ft1
      !! Field information from the outer grid
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(IN   ) :: F_nest 

      integer i,j,k
      integer bk, bj, bi ! Loop bound variables
      ! real(kind=REAL64) pri
!
!----------------------------------------------------------------------
!      
      !Spn_wrk(1:l_ni,1:l_nj,1:l_nk) = F_nest(1:l_ni,1:l_nj,1:l_nk) - F_ft1(1:l_ni,1:l_nj,1:l_nk)

      ! Subtract the current field from the nesting field to determine the unmodified
      ! correction.  Both fields are defined on the tracer grid, and this grid potentially
      ! includes a piloting region that we want to exclude.  The required offset has been
      ! precomputed in the Spn_pil_* variables, and we need Spn_pil_w for the i (x) index
      ! and Spn_pil_s for the j (y) index.
      bk = Spn_zgrid_lnz
      bj = Spn_zgrid_lny
      bi = Spn_zgrid_lnx
      !$omp do collapse(2)
      do k=1,Spn_zgrid_lnz
         do j=1,Spn_zgrid_lny
            do i=1,Spn_zgrid_lnx
               Spn_wrk(i,j,k) = F_nest(i+Spn_pil_w,j+Spn_pil_s,k) - F_ft1(i+Spn_pil_w,j+Spn_pil_s,k)
            end do
         end do
      end do

      ! Perform the z->x transpose.
      ! zx_tranpose%src_array is aliased to Spn_work, and that is the source array
      ! of the transpose.  zx_transpose%dst_array is aliased to Spn_xgrid, and this
      ! is the destination array of the transpose.
      call transpose_forward(zx_transpose)

      ! Perform the x-FFT in place, by executing the precomputed plan
      !$omp single
      call execute_r2r_dft_plan(fft_x_forward,Spn_xgrid,Spn_xgrid)
      !$omp end single

      ! Perform the x->y transpose
      ! The source array is Spn_xgrid (alised to xy_transpose%src_array), and
      ! the destination array is Spn_ygrid (aliased to %dst_array).
      call transpose_forward(xy_transpose)

      ! Perform the y-FFT in place.
      !$omp single
      call execute_r2r_dft_plan(fft_y_forward,Spn_ygrid,Spn_ygrid)
      !$omp end single

      ! Apply the filter, multiplying by Spn_flt.  Note that Spn_ygrid
      ! has (z,x,y) order, so the conventional loop indices are (k,i,j)

      !$omp do collapse(2)
      do j=1,Spn_ygrid_lny
         do i=1,Spn_ygrid_lnx
            do k=1,Spn_ygrid_lnz
               Spn_ygrid(k,i,j) = Spn_ygrid(k,i,j)*Spn_flt(i,j)
            end do
         end do
      end do

      ! Invert the y-FFT
      !$omp single
      call execute_r2r_dft_plan(fft_y_reverse,Spn_ygrid,Spn_ygrid)
      !$omp end single

      ! Invert the x->y transpose, writing the x-transformed values
      ! back to Spn_xgrid
      call transpose_reverse(xy_transpose)

      ! Invert the x-FFT
      !$omp single
      call execute_r2r_dft_plan(fft_x_reverse,Spn_xgrid,Spn_xgrid)
      !$omp end single

      ! Invert the z->x transpose, putting the increment values
      ! in Spn_wrk
      call transpose_reverse(zx_transpose)

      ! Apply the adjustment to F_ft1, modifying it in-place.
      ! As with the initial caluclation of the increment, we want to work with
      ! the grid inside the global piloting region; the piloting offset must
      ! be added when addressing F_ft1.
      !$omp do collapse(2)
      do k=2,Spn_zgrid_lnz ! Skip the surface level
         do j=1,Spn_zgrid_lny
            do i=1,Spn_xgrid_lnx
               F_ft1(i+Spn_pil_w,j+Spn_pil_s,k) = &
                  F_ft1(i+Spn_pil_w,j+Spn_pil_s,k) + prof(k)*Spn_wrk(i,j,k)*Spn_weight
            end do
         end do
      end do

!----------------------------------------------------------------------
!
      return
      end subroutine spn_apply

