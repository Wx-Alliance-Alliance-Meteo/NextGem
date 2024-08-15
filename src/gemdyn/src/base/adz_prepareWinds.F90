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

      subroutine adz_prepareWinds (F_itpc)
      use adz_mem
      use dcst
      use glb_ld
      use gmm_vt0
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: F_itpc

      integer :: i,j,k
      real(kind=REAL64), parameter :: alpha1= -1.d0/16.d0, &
                                      alpha2=  9.d0/16.d0
!
!     ---------------------------------------------------------------
!
      do k= 1, l_nk
         do j= 1, l_nj
            do i= 1, l_ni
              Adz_uu_arr(i,j,k)=((ut0(i-2,j,k) + ut0(i+1,j,k))*alpha1 &
                                +(ut0(i  ,j,k) + ut0(i-1,j,k))*alpha2)&
                                * Dcst_inv_rayt_8 * Adz_cy_8(j)
              Adz_vv_arr(i,j,k)=((vt0(i,j-2,k) + vt0(i,j+1,k))*alpha1 &
                                +(vt0(i  ,j,k) + vt0(i,j-1,k))*alpha2)&
                                * Dcst_inv_rayt_8
            end do
         end do
      end do
      do j= 1, l_nj
         do i= 1, l_ni
            Adz_ww_arr(i,j,1)= Adz_vw5 * zdt0(i,j,1)
            Adz_ww_arr(i,j,2)= adz_vwt2m(2,2) * zdt0(i,j,1) + &
                               adz_vwt2m(3,2) * zdt0(i,j,2) + &
                               adz_vwt2m(4,2) * zdt0(i,j,3)
            Adz_ww_arr(i,j,l_nk)= adz_vwt2m(1,l_nk) * zdt0(i,j,l_nk-2) + &
                                  adz_vwt2m(2,l_nk) * zdt0(i,j,l_nk-1) + &
                                  adz_vwt2m(3,l_nk) * zdt0(i,j,l_nk  )
         end do
      end do

      do k= 3, l_nk-1
         do j = 1,l_nj
            do i = 1,l_ni
               Adz_ww_arr(i,j,k)= adz_vwt2m(1,k) * zdt0(i,j,k-2) + &
                                  adz_vwt2m(2,k) * zdt0(i,j,k-1) + &
                                  adz_vwt2m(3,k) * zdt0(i,j,k  ) + &
                                  adz_vwt2m(4,k) * zdt0(i,j,k+1)
            end do
         end do
      end do

      if (F_itpc == 1) then
         do k= 1, l_nk
            do j= 1-Adz_haloy, l_nj+Adz_haloy
            do i= 1-Adz_halox, l_ni+Adz_halox
               Adz_uvw_d(1,i,j,k)= Dcst_inv_rayt_8 * Adz_uu_ext(i,j,k) &
                                                   * Adz_cy_8(j)
               Adz_uvw_d(2,i,j,k)= Dcst_inv_rayt_8 * Adz_vv_ext(i,j,k)
            end do
            end do
         end do
         do j= 1-Adz_haloy, l_nj+Adz_haloy
         do i= 1-Adz_halox, l_ni+Adz_halox
            Adz_uvw_d(3,i,j,1)= Adz_vw5 * Adz_ww_ext(i,j,1)
            Adz_uvw_d(3,i,j,2)= adz_vwt2m(2,2) * Adz_ww_ext(i,j,1) + &
                                adz_vwt2m(3,2) * Adz_ww_ext(i,j,2) + &
                                adz_vwt2m(4,2) * Adz_ww_ext(i,j,3)
            Adz_uvw_d(3,i,j,l_nk)= adz_vwt2m(1,l_nk) * Adz_ww_ext(i,j,l_nk-2) + &
                                   adz_vwt2m(2,l_nk) * Adz_ww_ext(i,j,l_nk-1) + &
                                   adz_vwt2m(3,l_nk) * Adz_ww_ext(i,j,l_nk  )
         end do
         end do      
         do k= 3, l_nk-1
            do j =1-Adz_haloy, l_nj+Adz_haloy
            do i =1-Adz_halox, l_ni+Adz_halox
               Adz_uvw_d(3,i,j,k)= adz_vwt2m(1,k) * Adz_ww_ext(i,j,k-2) + &
                                   adz_vwt2m(2,k) * Adz_ww_ext(i,j,k-1) + &
                                   adz_vwt2m(3,k) * Adz_ww_ext(i,j,k  ) + &
                                   adz_vwt2m(4,k) * Adz_ww_ext(i,j,k+1)
            end do
            end do
         end do
      endif
!
!     ---------------------------------------------------------------
!
      return

      end subroutine adz_prepareWinds
