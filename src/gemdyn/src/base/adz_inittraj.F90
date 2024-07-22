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

      subroutine adz_inittraj () 
      use glb_ld
      use adz_options
      use adz_mem

      integer i,j,k
!
!---------------------------------------------------------------------
!
!!$omp single
      Adz_niter = max( 6, Adz_itraj )
!!$omp end single
!!$omp do
     do k = Adz_k0m, l_nk
         do j = 1, l_nj
            do i = 1, l_ni
               Adz_pxyzm(1,i,j,k) = i + l_i0 - 1
               Adz_wpxyz(i,j,k,1) = Adz_pxyzm(1,i,j,k)
               Adz_pxyzm(2,i,j,k) = j + l_j0 - 1
               Adz_wpxyz(i,j,k,2) = Adz_pxyzm(2,i,j,k)
               Adz_pxyzm(3,i,j,k) = k
               Adz_wpxyz(i,j,k,3) = Adz_pxyzm(3,i,j,k)
               Adz_wpz(i,j,k)     = k

               !initialize displacement to 0.
               Adz_disp(i,j,k,:) = 0.d0

               !initialize previous midpoint to 0
               Adz_prev_mid(i,j,k,1) = 0.d0
               Adz_prev_mid(i,j,k,2) = 0.d0
               Adz_prev_mid(i,j,k,3) = 0.d0
            end do
         end do
      end do
!!$omp end do
      do j = 1, l_nj
         do i = 1, l_ni
            Adz_wpxyz(i,j,l_nk+1,1) = i + l_i0 - 1
            Adz_wpxyz(i,j,l_nk+1,2) = j + l_j0 - 1
            Adz_wpxyz(i,j,l_nk+1,3) = k
         end do
      end do
      
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_inittraj
