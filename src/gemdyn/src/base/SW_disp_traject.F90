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
!	New arrays added:
!		
!	 1. Adz_disp: holds the displacement alpha from arrival to midpoint
!		follows the declaration and allocation as Adz_wpxyz
!	 2. Adz_dep:   follows dec and alloc of Adz_pm
!	 3. Adz_pxyzd: follows dec and alloc of Adz_pxyzm
!	 4. Adz_dpxyz: follows dec and alloc of Adz_wpxyz
!	 5. Adz_uvw_lastl: follows declaration and allocation of Adz_uvw_dep
!------------------------------------------------------------------------------
      subroutine SW_disp_traject (F_dt_8)
      use glb_ld
      use cstv
      use geomh
      use ver
      use dyn_fisl_options
      use adz_options
      use adz_mem
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nb,k00
      integer,dimension(l_ni) :: kk
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8
      real(kind=REAL64) pos_d, pos_m, x_est, y_est, z_est, &
                        xd_est, yd_est, zd_est
!
!     ---------------------------------------------------------------
!
      if (Schm_advec == 1) then ! traditional advection
         dtA_8  = F_dt_8 * 0.5d0
         dtzA_8 = F_dt_8 * 0.5d0
      end if
      if (Schm_advec == 2) then ! consistent advection
         dtA_8  = F_dt_8
         dtzA_8 = F_dt_8
      end if
      if (Schm_advec == 3) then ! reversed advection
         dtA_8  = 0.
         dtzA_8 = 0.
      end if

      dtD_8  = F_dt_8 - dtA_8
      dtzD_8 = F_dt_8 - dtzA_8

      if (Schm_advec == 0) then ! no advection
         dtA_8  = 0.d0
         dtD_8  = 0.d0
         dtzA_8 = 0.d0
         dtzD_8 = 0.d0
      end if

      call adz_prepareWinds ()

      k00=Adz_k0m
      if (Adz_k0>1) k00=1

!Using the method from Staniforth and Cote. Standard displacement
!3-time-level scheme

!First, iterate to find the displacement alpha
      do  iter = 1, Adz_niter
!!$omp do
         do k= k00, l_nk

             !1. interpolate velocity at current time to estimated midpoints
             call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                     Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)

             !2. loop through i,j indices to find displacement for midpoint,
             !i.e find alpha
             ! alpha ^(k+1) = dt*V(x - alpha^k, t)
             do j= 1, l_nj
                do i= 1, l_ni

                   !compute displacement 
                   Adz_disp(i,j,k,1) = Cstv_dt_8*Adz_uvw_dep(1,i,j,k)
                   Adz_disp(i,j,k,2) = Cstv_dt_8*Adz_uvw_dep(2,i,j,k)
                   Adz_disp(i,j,k,3) = Cstv_dt_8*Adz_uvw_dep(3,i,j,k)
                end do !end i 

                !3. estimate midpoint
                !mid = arrival - disp
                do i= 1, l_ni

                   x_est = dble(i+l_i0-1) - Adz_disp(i,j,k,1)*geomh_inv_hx_8
                   y_est = dble(j+l_j0-1) - Adz_disp(i,j,k,2)*geomh_inv_hy_8
                   pos_m = Ver_z_8%m(k)   - Adz_disp(i,j,k,3)

                   !--extra stuff for z--
                   z_est = min(max(pos_m,Ver_zmin_8),Ver_zmax_8)

                   kk1 = (z_est - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                   kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                   kk1 = adz_search_m(kk1)
                   if ( sig * z_est > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                   if ( sig * z_est < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                   kk(i) = kk1
                   !----------------------
 
                   !----set the midpoint----
                   Adz_wpxyz(i,j,k,1) = x_est
                   Adz_pxyzm(1,i,j,k) = min(max(x_est,Adz_iminposx),&
                                                     Adz_imaxposx)
                   Adz_wpxyz(i,j,k,2) = y_est
                   Adz_pxyzm(2,i,j,k) = min(max(y_est,Adz_iminposy),&
                                                     Adz_imaxposy)
                   kk1 = min(l_nk+1,max(0,kk(i)))
                   nb  = max(min(kk1,G_nk-1),1)
                   Adz_wpxyz(i,j,k,3) = (z_est-ver_z_8%m(nb))&
                                      *Adz_odelz_m(nb) + dble(nb)
                   Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
                   !----done----

                end do ! end i
            end do     ! end j
         enddo         ! end k        
!!$omp enddo
      end do          !end adz_iter
      

!once the displacement is found, compute the departure point
! departure = arrival_point - 2*displacement
!!$omp do
         do k= k00, l_nk
           do j= 1, l_nj
             do i= 1, l_ni
                !1. compute departure point
                 xd_est = dble(i+l_i0-1) - 2.d0*Adz_disp(i,j,k,1)*geomh_inv_hx_8 
                 yd_est = dble(j+l_j0-1) - 2.d0*Adz_disp(i,j,k,2)*geomh_inv_hy_8
                 pos_d  = Ver_z_8%m(k)   - 2.d0*Adz_disp(i,j,k,3)

                 !---extra stuff for z---
                 zd_est = min(max(pos_d,Ver_zmin_8),Ver_zmax_8)

                 kk1 = (zd_est - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                 kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                 kk1 = adz_search_m(kk1)
                 if ( sig * zd_est > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                 if ( sig * zd_est < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                 kk(i) = kk1
                 !----------------------

                 !---set the departure points---
                 Adz_dpxyz(i,j,k,1) = xd_est
                 Adz_pxyzd(1,i,j,k) = min(max(xd_est,Adz_iminposx),&
                                                     Adz_imaxposx)
                 Adz_dpxyz(i,j,k,2) = yd_est
                 Adz_pxyzd(2,i,j,k) = min(max(yd_est,Adz_iminposy),&
                                                    Adz_imaxposy)
                 kk1 = min(l_nk+1,max(0,kk(i)))
                 nb = max(min(kk1,G_nk-1),1)
                 Adz_dpxyz(i,j,k,3) = (zd_est-ver_z_8%m(nb))&
                                      *Adz_odelz_m(nb) + dble(nb)
                 Adz_pxyzd(3,i,j,k) = Adz_dpxyz(i,j,k,3)
                 !----done----

              end do !end i
           end do    !end j
        end do       !end k
!!$omp end do

!now compute the velocity at the departure point
!the data used for interpolation is no longer Adz_uvw_d, since that
!correspoinds to pw_plus winds (at current time), we want the velocity
!at the previous time, i.e pw_moins winds
!Note: the data Adz_uvw_d_dep follows Adz_uvw_d (original); also 
!initialized in adz_prepareWinds

!!$omp do
       do k= k00, l_nk
           call tricublin_zyx3_n ( Adz_uvw_lastl(1,1,1,k),Adz_uvw_d_dep(1,1,1,1), &
                                   Adz_pxyzd(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
       end do
!!$omp enddo


!!$OMP BARRIER

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

!---for midpoint---
!!$omp single
      Adz_niter = Adz_itraj
      call rpn_comm_xch_halo_8 (Adz_wpxyz, -1,l_ni+2, -1,l_nj+2,&
                 l_ni,l_nj, 3*l_nk, 2,2, .false.,.false., l_ni,0)
!!$omp end single

!---for departure---
!!$omp single
      call rpn_comm_xch_halo_8 (Adz_dpxyz, -1,l_ni+2, -1,l_nj+2,&
                 l_ni,l_nj, 3*l_nk, 2,2, .false.,.false., l_ni,0)
!!$omp end single

!like adz_interp_traj, but added code to account for the departure point
!at previous time level                 
      call BDF_interp_traj (dtzD_8, dtzA_8, F_dt_8)

!     ---------------------------------------------------------------
!
      return
      end subroutine SW_disp_traject
