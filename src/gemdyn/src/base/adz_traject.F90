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

      subroutine adz_traject (F_dt_8)
      use ISO_C_BINDING
      use glb_ld
      use geomh
      use ver
      use dyn_fisl_options
      use adz_options
      use adz_mem
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nb,k00,dim
      integer :: HLT_np, HLT_start, HLT_end
      integer,dimension(l_ni) :: kk
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64), dimension(:,:), pointer :: xchg
      real(kind=REAL64) pos
      type(C_PTR) :: Cpntr
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

      do  iter = 1, Adz_niter
!!$omp do
         do k= k00, l_nk
             call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                     Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - &
                    ( dtD_8 * Adz_uvw_dep(1,i,j,k)  &
                    + dtA_8 * Adz_uu_arr(i,j,k) )&
                    * geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - &
                    ( dtD_8 * Adz_uvw_dep(2,i,j,k)  &
                    + dtA_8 * Adz_vv_arr(i,j,k) )&
                    * geomh_inv_hy_8
                  pos = Ver_z_8%m(k)- dtzD_8* Adz_uvw_dep(3,i,j,k) &
                                    - dtzA_8* Adz_ww_arr(i,j,k)
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                  kk(i) = kk1
               end do
               do i= 1, l_ni
                  Adz_wpxyz(i,j,k,1) = xm(i)
                  Adz_pxyzm(1,i,j,k) = min(max(xm(i),Adz_iminposx),&
                                                     Adz_imaxposx)
                  Adz_wpxyz(i,j,k,2) = ym(i)
                  Adz_pxyzm(2,i,j,k) = min(max(ym(i),Adz_iminposy),&
                                                     Adz_imaxposy)
                  kk1 = min(l_nk+1,max(0,kk(i)))
                  nb  = max(min(kk1,G_nk-1),1)
                  Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                     *Adz_odelz_m(nb) + dble(nb)
                  Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
               end do
            end do
         enddo
!!$omp enddo
      end do
      
!!$omp do
      do k= Adz_k0, l_nk
         Adz_pm   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzm(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!!$omp enddo nowait

      dim=3*ubound(Adz_wpxyz,3)
      call HLT_split (1, dim, HLT_np, HLT_start, HLT_end)
      Cpntr= c_loc( Adz_wpxyz(-1,-1,1,1) )
      call C_F_POINTER ( Cpntr, xchg, [(l_ni+4)*(l_nj+4),dim] )
      call gem_xch_halo_8 ( xchg(1,HLT_start),&
                    -1,l_ni+2,-1,l_nj+2, HLT_np,-1)      
!!$omp single
      Adz_niter = Adz_itraj
!!$omp end single
                 
      call adz_interp_traj (dtzD_8, dtzA_8, F_dt_8)

!     ---------------------------------------------------------------
!
      return
      end subroutine adz_traject
