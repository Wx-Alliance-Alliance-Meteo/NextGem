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

      logical function picard_stop (F_dt_8, F_iter, print_conv)
      use dcst
      use glb_ld
      use glb_pil
      use gmm_vt0
      use geomh
      use lun
      use ver
      use dyn_fisl_options
      use adz_mem
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(IN) :: print_conv
      integer, intent(IN) :: F_iter
      real(kind=REAL64), intent(IN) :: F_dt_8

      include 'mpif.h'
      integer :: i,j,k,kk1,ierr
      real(kind=REAL64) :: inv_dt_8, N
      real :: rate
      real, dimension(l_ni) :: uu_arr,vv_arr
      real, dimension(l_ni,l_nj,l_nk) :: ww_arr
      real(kind=REAL64) :: rx,ra,rb,rc,rd,w1,w2,w3,w4,t1,t2,t3
      real(kind=REAL64) :: err, err2nrm, glb, total_sum
      real(kind=REAL64), dimension(l_ni) :: r1,r2,r3,s1,s2,s3,lenght
      real(kind=REAL64), parameter :: one=1.d0, two=2.d0, half=0.5d0, to2=3.d0/2.d0, epsi=1.d-15
      real(kind=REAL64), parameter :: alpha1= -1.d0/16.d0, alpha2= 9.d0/16.d0
      real(kind=REAL64), dimension(100), save :: sucsv_err
!     ---------------------------------------------------------------
!
      picard_stop= .false.
      
      ! Estimate L_inf error
      ! || 1/tau( 3/2 x^{+} - 2x^{0} + 1/2 x^{-} ) - V^{+} ||_{inf}
      ! where x^{0} and x^{+} are the current estimates of the
      ! mid and departure points

      ! prepare the winds at arrival point for current iterate
      rx = Ver_z_8%m(2)
      ra = Ver_z_8%t(0)
      rb = Ver_z_8%t(1)
      rc = Ver_z_8%t(2)
      rd = Ver_z_8%t(3)
      t1 = lag3(rx, rb, ra, rc, rd)
      t2 = lag3(rx, rc, ra, rb, rd)
      t3 = lag3(rx, rd, ra, rb, rc)
      rx = Ver_z_8%m(l_nk)
      ra = Ver_z_8%t(l_nk-2)
      rb = Ver_z_8%t(l_nk-1)
      rc = Ver_z_8%t(l_nk)
      rd = Ver_z_8%t(l_nk+1)
      w1 = lag3(rx, ra, rb, rc, rd)
      w2 = lag3(rx, rb, ra, rc, rd)
      w3 = lag3(rx, rc, ra, rb, rd)
      do j= 1+pil_s, l_nj-pil_n
         do i= 1+pil_w, l_ni-pil_e
            ww_arr(i,j,   1)= Adz_vw5 * zdt0(i,j,1)
            ww_arr(i,j,   2)= t1*zdt0(i,j,1) + t2*zdt0(i,j,2) + t3*zdt0(i,j,3)
            ww_arr(i,j,l_nk)= w1*zdt0(i,j,l_nk-2) + w2*zdt0(i,j,l_nk-1) + w3*zdt0(i,j,l_nk  )
         end do
      end do
      do k= 3, l_nk-1
         rx = Ver_z_8%m(k)
         ra = Ver_z_8%t(k-2)
         rb = Ver_z_8%t(k-1)
         rc = Ver_z_8%t(k)
         rd = Ver_z_8%t(k+1)
         w1 = lag3(rx, ra, rb, rc, rd)
         w2 = lag3(rx, rb, ra, rc, rd)
         w3 = lag3(rx, rc, ra, rb, rd)
         w4 = lag3(rx, rd, ra, rb, rc)         
         do j= 1+pil_s, l_nj-pil_n
            do i= 1+pil_w, l_ni-pil_e
               ww_arr(i,j,k) = w1 * zdt0(i,j,k-2) + w2 * zdt0(i,j,k-1) + &
                               w3 * zdt0(i,j,k  ) + w4 * zdt0(i,j,k+1)
            end do
         end do
      end do
      
      inv_dt_8 = 1.0/F_dt_8
      err = -huge(1.) ; total_sum = 0.d0 ; lenght=0.d0
      N = ((G_ni-Glb_pil_e)-(1+Glb_pil_w)+1)*((G_nj-Glb_pil_n)-(1+Glb_pil_s)+1)*G_nk

      do k= 1, l_nk
         do j= 1+pil_s, l_nj-pil_n
            do i= 1+pil_w, l_ni-pil_e
               uu_arr(i)=((ut0(i-2,j,k) + ut0(i+1,j,k))*alpha1 &
                         +(ut0(i  ,j,k) + ut0(i-1,j,k))*alpha2)&
                                * Dcst_inv_rayt_8 * Adz_cy_8(j)
               vv_arr(i)=((vt0(i,j-2,k) + vt0(i,j+1,k))*alpha1 &
                         +(vt0(i  ,j,k) + vt0(i,j-1,k))*alpha2)&
                                             * Dcst_inv_rayt_8
               r1(i) = to2*dble(i+l_i0-1) - two*Adz_wpxyz(i,j,k,1) + half*Adz_dpxyz(i,j,k,1)
               r2(i) = to2*dble(j+l_j0-1) - two*Adz_wpxyz(i,j,k,2) + half*Adz_dpxyz(i,j,k,2)
               r3(i) = k + to2*k - two*Adz_wpxyz(i,j,k,3) + half*Adz_dpxyz(i,j,k,3)
               s1(i) = F_dt_8*uu_arr(i)*geomh_inv_hx_8
               s2(i) = F_dt_8*vv_arr(i)*geomh_inv_hy_8
               s3(i) = Ver_z_8%m(k) + F_dt_8*ww_arr(i,j,k)
               kk1 = (s3(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
               kk1 = min(max(1,kk1),ubound(adz_search_m,1))
               kk1 = adz_search_m(kk1)
               s3(i) = kk1+(Ver_z_8%m(kk1)-s3(i))/(Ver_z_8%m(kk1)-Ver_z_8%m(kk1+1))
               lenght(i)= sqrt((r1(i) - s1(i))**2. + (r2(i) - s2(i))**2. + (r3(i) - s3(i))**2.)
           end do
           total_sum = total_sum + sum(lenght(1+pil_w:l_ni-pil_e))
           err = max(err,maxval(lenght))
        end do
      enddo

      call MPI_allreduce ( err       , glb, 1, MPI_DOUBLE_PRECISION,MPI_MAX,COMM_MULTIGRID,ierr )
      err= glb
      call MPI_allreduce ( total_sum , glb, 1, MPI_DOUBLE_PRECISION,MPI_SUM,COMM_MULTIGRID,ierr )
      err2nrm = glb/N

      sucsv_err(F_iter)= err2nrm
      rate=1.d0
      if (F_iter > 1) rate= (sucsv_err(F_iter-1)-sucsv_err(F_iter))/sucsv_err(F_iter-1)
      if ((err2nrm < Schm_tolpic).or.(rate<Schm_ratepic)) picard_stop= .true.
      if (print_conv) write(Lun_out,'(" PICARD ",i1,": DIFF_MAX,L_2,rate: ",2(1pe11.4),1pe9.2)')&
                      F_iter,err,err2nrm,rate
!
!     ---------------------------------------------------------------
!
      return
contains
      real(kind=REAL64) function lag3(zz, z1, z2, z3, z4)
      implicit none
      real(kind=REAL64), intent(in)  :: zz, z1, z2, z3, z4
      
      lag3 = ( ((zz - z2) * (zz - z3) * (zz - z4) ) / &
               ((z1 - z2) * (z1 - z3) * (z1 - z4) ) )
      end function lag3
      
      end function picard_stop
