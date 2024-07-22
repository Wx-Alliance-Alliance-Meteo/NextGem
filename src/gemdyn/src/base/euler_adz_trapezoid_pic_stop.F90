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

      subroutine euler_adz_trapezoid_pic_stop (F_dt_8)
      use ISO_C_BINDING
      use glb_ld
      use cstv
      use geomh
      use gmm_vt0
      use gmm_vt1
      use ver
      use dyn_fisl_options
      use HORgrid_options
      use adz_options
      use adz_mem
      use omp_lib
      use omp_timing
      use ptopo
      use stat_mpi, only:statf_dm
      use, intrinsic :: iso_fortran_env
      implicit none
      
      include 'mpif.h'

      real(kind=REAL64), intent(IN) :: F_dt_8

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nb,k00, ierr1, ierr2,dim
      integer :: HLT_np, HLT_start, HLT_end
      integer,dimension(l_ni) :: kk
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8,quarter_dt_8,inv_dt_8,pos
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64), dimension(:,:), pointer :: xchg
      type(C_PTR) :: Cpntr
      !----for stopping criteria----
      real(kind=REAL64), parameter :: max_iter = 10, tol=1e-6 
      real(kind=REAL64) :: max_x, max_y, max_z, local_max, einf
      real(kind=REAL64) :: x_sum, y_sum, z_sum, local_sum, e2
      real(kind=REAL64), dimension(11) :: err, err2nrm, sucsv_err, sucsv_err2nrm
      real(kind=REAL64), dimension(l_ni) :: xerr, yerr, zerr
      real(kind=REAL64) :: N, glb_dot_pic, glb_max_pic, l_avg_8_pic, l_max_8_pic
      real(kind=REAL64), dimension(:,:), allocatable :: thread_sum_pic, thread_max_pic
!
!     ---------------------------------------------------------------
!
      allocate (thread_sum_pic(1,0:OMP_get_max_threads()-1), &
                thread_max_pic(1,0:OMP_get_max_threads()-1))
      thread_sum_pic = 0.d0; thread_max_pic = 0.d0

      k00=Adz_k0m
      if (Adz_k0>1) k00=1

      quarter_dt_8= 0.25d0 *Cstv_dt_8

      N = l_nk * G_nj * G_ni

      call adz_prepareWinds()

      do iter = 1, max_iter !iterate at most up to 10

!!$omp do
        !0a.Interpolate velocity V(t) at current time to estimated midpoints x^{n}: V^{n}
        do k= k00, l_nk
           call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                   Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
        end do
!!$omp end do

!!$omp do
        !0b.Estimate midpoints as x^{n-1}=x^{n+1}- 0.25*dt*(V^{n} + V^{+})  using velocity
        !interpolated at the midpoints (0a)
           do k= k00, l_nk
              do j= 1, l_nj
                 do i= 1, l_ni
                    xm(i) = dble(i+l_i0-1)-quarter_dt_8*(Adz_uvw_dep(1,i,j,k) + Adz_uu_arr(i,j,k))*geomh_inv_hx_8
                    ym(i) = dble(j+l_j0-1)-quarter_dt_8*(Adz_uvw_dep(2,i,j,k) + Adz_vv_arr(i,j,k))*geomh_inv_hy_8
                    pos   = Ver_z_8%m(k)  -quarter_dt_8*(Adz_uvw_dep(3,i,j,k) + Adz_ww_arr(i,j,k))
                    zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                    kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                    kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                    kk1 = adz_search_m(kk1)
                    if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                    if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                    kk(i) = kk1
                 enddo
                 do i= 1, l_ni
                    Adz_wpxyz(i,j,k,1) = xm(i)
                    Adz_pxyzm(1,i,j,k) = min(max(xm(i),Adz_iminposx),&
                                                       Adz_imaxposx)
                    Adz_wpxyz(i,j,k,2) = ym(i)
                    Adz_pxyzm(2,i,j,k) = min(max(ym(i),Adz_iminposy),&
                                                     Adz_imaxposy)
                    kk1 = min(l_nk+1,max(0,kk(i)))
                    nb  = max(min(kk1,G_nk-1),1)
                    Adz_wpz(i,j,k) = zm(i)
                    Adz_dpz(i,j,k) = zm(i) !dep = midpoint for euler
                    Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                        *Adz_odelz_m(nb) + dble(nb)
                    Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
                 end do
              enddo
            enddo
!!$omp end do

         !error = succssive iterates ||x^0_{k+1} - x^0_{k}|| 
         !compute infinity and 2-nrm error
         x_sum = 0.d0; y_sum = 0.d0; z_sum = 0.d0; local_sum = 0.d0; local_max = 0.d0

!!$omp do
         do k= k00, l_nk
           do j= 1, l_nj
             do i= 1, l_ni
                !note: previous midpoint estimate = 0, set in adz_inittraj

                !---for inf-norm---
                xerr(i) = abs( Adz_wpxyz(i,j,k,1) - Adz_prev_mid(i,j,k,1) )
                yerr(i) = abs( Adz_wpxyz(i,j,k,2) - Adz_prev_mid(i,j,k,2) )
                zerr(i) = abs( Adz_wpz(i,j,k)     - Adz_prev_mid(i,j,k,3) ) 

                !---for 2-norm---
                x_sum = ( Adz_wpxyz(i,j,k,1) - Adz_prev_mid(i,j,k,1) )**2
                y_sum = ( Adz_wpxyz(i,j,k,2) - Adz_prev_mid(i,j,k,2) )**2
                z_sum = ( Adz_wpz(i,j,k)     - Adz_prev_mid(i,j,k,3) ) **2
                local_sum =  local_sum + x_sum + y_sum + z_sum 

              end do

               !max value of x,y,z array
               max_x = maxval(xerr)
               max_y = maxval(yerr)
               max_z = maxval(zerr)

               !which direction has the max error
               local_max = max(local_max, max(max_x, max(max_y,max_z)))
           end do
         enddo
!!$omp end do

      !the parallel portion is modeled after sol_fgmres for
      !computing the inner products
      thread_sum_pic(1, OMP_get_thread_num()) = local_sum
      thread_max_pic(1, OMP_get_thread_num()) = local_max

!$OMP BARRIER

!!$omp single
      l_avg_8_pic = sum(thread_sum_pic(1,:))
      l_max_8_pic = maxval(thread_max_pic(1,:))

      call MPI_allreduce(l_avg_8_pic, glb_dot_pic, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM_MULTIGRID, ierr1) 
      call MPI_allreduce(l_max_8_pic, glb_max_pic, 1, MPI_DOUBLE_PRECISION, MPI_max, COMM_MULTIGRID, ierr2) 

      e2   = sqrt(glb_dot_pic/N)
      einf = glb_max_pic
!!$omp end single copyprivate(e2, einf) 

         !now save inf-norm error and 2-norm error
         err(iter)     = einf
         err2nrm(iter) = e2

         !print *, iter, err(iter), err2nrm(iter)

         !break out of loop if we reached convergence in inf-norm
         if(err(iter) < tol) then
         !if(err2nrm(iter) < tol) then
           exit
         endif

         !or break if we've stagnated
         if (iter > 1) then 
           !1. compute successive errors
           sucsv_err(iter)     = abs(err(iter) - err(iter-1))
           sucsv_err2nrm(iter) = abs(err2nrm(iter) - err2nrm(iter-1))
           !print *, iter, err(iter),err2nrm(iter), sucsv_err(iter), sucsv_err2nrm(iter)
           if (sucsv_err(iter) < err(iter) ) then 
           !if (sucsv_err2nrm(iter) < err2nrm(iter) ) then 
             exit
           endif
         endif

         !save the current estimate for next iteration
!!$omp do 
         do k= k00, l_nk
           do j= 1, l_nj
             do i= 1, l_ni
                !save current midpoint
                Adz_prev_mid(i,j,k,1) = Adz_wpxyz(i,j,k,1)
                Adz_prev_mid(i,j,k,2) = Adz_wpxyz(i,j,k,2)
                Adz_prev_mid(i,j,k,3) = Adz_wpz(i,j,k)
              end do
           end do
         enddo
!!$omp enddo

      end do !Adz_niter

!!$omp do
      !set departure points == midpoints
         do k= k00, l_nk
            do j= 1, l_nj
               do i= 1, l_ni
                  Adz_dpxyz(i,j,k,1) = Adz_wpxyz(i,j,k,1)
                  Adz_pxyzd(1,i,j,k) = Adz_pxyzm(1,i,j,k)
                  Adz_dpxyz(i,j,k,2) = Adz_wpxyz(i,j,k,2)
                  Adz_pxyzd(2,i,j,k) = Adz_pxyzm(2,i,j,k)
                  Adz_dpxyz(i,j,k,3) = Adz_wpxyz(i,j,k,3)
                  Adz_pxyzd(3,i,j,k) = Adz_pxyzm(3,i,j,k) 
               end do
            enddo
          enddo
!!$omp end do
      
!!$omp do
      do k= Adz_k0, l_nk
         Adz_pm   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzm(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!!$omp enddo nowait

!---for departure point---      
!!$omp do
      do k= Adz_k0, l_nk
         Adz_dep  (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzd(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!!$omp enddo nowait

      dim=3*ubound(Adz_wpxyz,3)
      call HLT_split (1, dim, HLT_np, HLT_start, HLT_end)
      Cpntr= c_loc( Adz_wpxyz(-1,-1,1,1) )
      call C_F_POINTER ( Cpntr, xchg, [(l_ni+4)*(l_nj+4),dim] )
      call gem_xch_halo_8 ( xchg(1,HLT_start),&
                    -1,l_ni+2,-1,l_nj+2, HLT_np,-1)      
      dim=3*ubound(Adz_dpxyz,3)
      call HLT_split (1, dim, HLT_np, HLT_start, HLT_end)
      Cpntr= c_loc( Adz_dpxyz(-1,-1,1,1) )
      call C_F_POINTER ( Cpntr, xchg, [(l_ni+4)*(l_nj+4),dim] )
      call gem_xch_halo_8 ( xchg(1,HLT_start),&
                    -1,l_ni+2,-1,l_nj+2, HLT_np,-1)

!like adz_interp_traj, but added code to account for the departure point
!at previous time level                 
      call BDF_interp_traj (dtzD_8, dtzA_8, F_dt_8)

!     ---------------------------------------------------------------
!
      return
      end subroutine euler_adz_trapezoid_pic_stop
