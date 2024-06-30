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

      subroutine check_picard_stop (F_dt_8, err, err2nrm)
      use glb_ld
      use glb_pil
      use geomh
      use ver
      use dyn_fisl_options
      use adz_mem
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN)  :: F_dt_8
      real(kind=REAL64), intent(OUT) :: err, err2nrm

      include 'mpif.h'
      integer :: i,j,k,k00,ierr
      real(kind=REAL64) :: inv_dt_8, N
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      !----for stopping criteria------
      real(kind=REAL64) :: max_x, max_y, max_z, tmp, glb !<-- hold max value in eachdirection
      real(kind=REAL64) :: x_sum, y_sum, z_sum, total_sum !<-- hold the squared sum in each direction
      real(kind=REAL64), dimension(l_ni) :: xerr, yerr, zerr !<-- error in each direction
!
!     ---------------------------------------------------------------
!
      inv_dt_8 = 1.0/F_dt_8

      call adz_prepareWinds ()

      k00=Adz_k0m
      if (Adz_k0>1) k00=1

      !4. estimate L_inf error
      ! || 1/tau( 3/2 x^{+} - 2x^{0} + 1/2 x^{-} ) - V^{+} ||_{inf}
      ! where x^{0} and x^{+} are the current estimates of the
      ! mid and departure points
      xerr = 0.d0 ; yerr = 0.d0 ; zerr = 0.d0 ; err = -huge(1.)
      x_sum = 0.d0; y_sum = 0.d0; z_sum = 0.d0; total_sum = 0.d0
      N = ((G_ni-Glb_pil_e)-(1+Glb_pil_w)+1)*((G_nj-Glb_pil_n)-(1+Glb_pil_s)+1)*G_nk

      do k= k00, l_nk
         do j= 1+pil_s, l_nj-pil_n
            do i= 1+pil_w, l_ni-pil_e
             !---for infinity norm---
             xerr(i) = abs( inv_dt_8*((3.0/2.0)*dble(i+l_i0-1) - 2.0*Adz_wpxyz(i,j,k,1) + (1.0/2.0)*Adz_dpxyz(i,j,k,1)) - Adz_uu_arr(i,j,k)*geomh_inv_hx_8 )
             yerr(i) = abs( inv_dt_8*((3.0/2.0)*dble(j+l_j0-1) - 2.0*Adz_wpxyz(i,j,k,2) + (1.0/2.0)*Adz_dpxyz(i,j,k,2)) - Adz_vv_arr(i,j,k)*geomh_inv_hy_8 )
             zerr(i) = abs( inv_dt_8*((3.0/2.0)*Ver_z_8%m(k)   - 2.0*Adz_wpz(i,j,k) + (1.0/2.0)*Adz_dpz(i,j,k)) - Adz_ww_arr(i,j,k) ) 

             !--for 2-norm---
             x_sum = ( inv_dt_8*((3.0/2.0)*dble(i+l_i0-1) - 2.0*Adz_wpxyz(i,j,k,1) + (1.0/2.0)*Adz_dpxyz(i,j,k,1)) - Adz_uu_arr(i,j,k)*geomh_inv_hx_8 ) **2
             y_sum = ( inv_dt_8*((3.0/2.0)*dble(j+l_j0-1) - 2.0*Adz_wpxyz(i,j,k,2) + (1.0/2.0)*Adz_dpxyz(i,j,k,2)) - Adz_vv_arr(i,j,k)*geomh_inv_hy_8 ) **2
             z_sum = ( inv_dt_8*((3.0/2.0)*Ver_z_8%m(k)   - 2.0*Adz_wpz(i,j,k) + (1.0/2.0)*Adz_dpz(i,j,k)) - Adz_ww_arr(i,j,k) )**2  
             
             total_sum =  total_sum + x_sum + y_sum + z_sum 
           end do

            !max value of x,y,z array
           max_x = maxval(xerr)
           max_y = maxval(yerr)
           max_z = maxval(zerr)

            !which direction has the max error
           tmp = max(max_x, max(max_y, max_z))
           err = max(err, tmp)
        end do
      enddo

!     take square root of total error
      call MPI_allreduce ( err       , glb, 1, MPI_DOUBLE_PRECISION,MPI_MAX,COMM_MULTIGRID,ierr )
      err= glb
      call MPI_allreduce ( total_sum , glb, 1, MPI_DOUBLE_PRECISION,MPI_SUM,COMM_MULTIGRID,ierr )
      err2nrm = sqrt(glb/N)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine check_picard_stop
