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
!----------------------------------LICENCE END ---------------------------------
      
      subroutine fillVH !(0F_s, Minx, Maxx, Miny, Maxy)
      use glb_ld
      use gmm_vt1
      use dyn_fisl_options
      use gem_options
      use metric
      use sol_mem
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none
      
!      integer, intent(in) :: Minx, Maxx, Miny, Maxy
!      real, dimension(Minx:Maxx,Miny:Maxy), intent(inout ):: F_s

      integer :: i, j, k
      real(kind=REAL64) :: w5,ttt,q
!
!     ________________________________________________________________
!
      ext_t(:,:,1:l_nk) = tt1(:,:,1:l_nk)
      do k=-4,0
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               w5= (VM3%ztht(i,j,k)-VM3%ztht(i,j,2))/(VM3%ztht(i,j,2)-VM3%ztht(i,j,1))
               ext_t(i,j,k)= tt1(i,j,2)*(1+w5) -w5*tt1(i,j,1)
            end do
         end do
      end do
      do k= l_nk+1,l_nk+4
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               w5= (VM3%ztht(i,j,k)-VM3%ztht(i,j,l_nk))/(VM3%ztht(i,j,l_nk)-VM3%ztht(i,j,l_nk-1))
               ext_t(i,j,k)= tt1(i,j,l_nk)*(1+w5) -w5*tt1(i,j,l_nk-1)
            !   w5= VM3%ztht(i,j,l_nk)-VM3%ztht(i,j,k)
            !   ext_t(i,j,k) = tt1(i,j,l_nk) + stlo_8 * w5 !??? Schuman-Newel lapse rate
            end do
         end do
      end do
      
      do k= 0,-3,-1
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
            ttt= ext_t(i,j,k-1)
            if (k==0) ttt=  0.5d0*(ext_t(i,j,k) + ext_t(i,j,k-1))
!            q= qt1(i,j,k+1)-grav_8*VM3%zmom(i,j,k+1)
!            qt1(i,j,k)= q + grav_8*Cstv_tstr_8*(VM3%zmom(i,j,k+1) - VM3%zmom(i,j,k))/ttt + grav_8*VM3%zmom(i,j,k)
            q= (VM3%zmom(i,j,k+1)-VM3%zmom(i,j,k))/(VM3%zmom(i,j,k+2)-VM3%zmom(i,j,k+1))
            qt1(i,j,k)= (1+q)*qt1(i,j,k+1) -  q*qt1(i,j,k+2)
         end do
      end do
      end do
!!$      k= l_nk
!!$      ps=(rgasd_8*Cstv_tstr_8)*log(F_ps(i,j)/1.e5)
!!$      q= (VM3%zmom(i,j,k)-VM3%zmom(i,j,k-1))/(VM3%zmom(i,j,k-1)-VM3%zmom(i,j,k-2))
!!$      qt1(i,j,k)= (1+q)*qt1(i,j,k-1) - q*qt1(i,j,k-2)
      do k= l_nk+1, l_nk+4
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
            ttt= ext_t(i,j,k)
            if (k==l_nk+1) ttt=  0.5d0*(ext_t(i,j,k) + ext_t(i,j,k-1))
!            q= qt1(i,j,k-1)-grav_8*VM3%zmom(i,j,k-1)
!            qt1(i,j,k)= q + grav_8*Cstv_tstr_8*(VM3%zmom(i,j,k-1) - VM3%zmom(i,j,k))/ttt + grav_8*VM3%zmom(i,j,k)
            q= (VM3%zmom(i,j,k)-VM3%zmom(i,j,k-1))/(VM3%zmom(i,j,k-1)-VM3%zmom(i,j,k-2))
            qt1(i,j,k)= (1+q)*qt1(i,j,k-1) - q*qt1(i,j,k-2)
      end do
      end do
      end do
!
!     ________________________________________________________________
!
      return
      end subroutine fillVH
