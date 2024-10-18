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

!**s/r  bac_H- backsubstitution: obtain new values of the variables:u,v,w,t,q,zd
!                from new q , the right-hand sides (Ru,Rv,Rw,Rt,Rf)
!                             and non-linear terms (Nu,Nv,Nw,Nt)
!             - Height-type vertical coordinate

      subroutine SW_bac ( F_dt_8, i0, j0, k0, in, jn ,k0t )
      use gem_options
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use sol_mem
      use HORgrid_options
      use tdpack
      use gmm_vt0
      use mem_tstp
      use glb_pil
      use glb_ld
      use cstv
      use ver
      use metric
      use ctrl
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, km
      integer :: HLT_np, HLT_start, HLT_end, itpc, local_np
      real :: w5,w6
      real(kind=REAL64) :: tau_m_8, tau_nh_8, invT_m_8, invT_nh_8, Buoy
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0

!
!     ---------------------------------------------------------------
!
      if (Ctrl_testcases_adv_L) then
!!$omp single
         call canonical_cases ("BAC")
!!$omp end single
         return
      end if

      tau_m_8  = (2.d0*F_dt_8) / 3.d0 
      tau_nh_8 = (2.d0*F_dt_8) / 3.d0 
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

      !print *, "---------------------------------------------"

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
            do i= i0, in
               qt0(i,j,k) = sngl(Sol_lhs(i,j,l_nk))

               !if (k == l_nk) then 
               !  print *, qt0(i,j,k)
               !end if

            end do
         end do
      end do
!!$omp enddo

!!$omp do collapse(2)
         do j= j0, jn
            do i= i0, in
               qt0(i,j,l_nk+1) = sngl(Sol_lhs(i,j,l_nk))
            end do
         end do
!!$omp enddo


!!$omp single
!     call rpn_comm_xch_halo(qt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
!                            G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
!!$omp end single
      call HLT_split (1, G_nk+1, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( qt0(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
            km=max(k-1,1)
!DIR$ SIMD
            do i= i0, l_niu-pil_e

   !           Compute U
   !           ~~~~~~~~~
               ut0(i,j,k) = tau_m_8*(Ruu(i,j,k) - grav_8*(qt0(i+1,j,k)-qt0(i,j,k))*geomh_invDX_8(j)) 
            end do
         end do
      end do
!!$omp enddo nowait
      
!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, l_njv-pil_n
            km=max(k-1,1)
!DIR$ SIMD
            do i= i0, in

   !           Compute V
   !           ~~~~~~~~~
               vt0(i,j,k) = tau_m_8*(Rvv(i,j,k) - grav_8*(qt0(i,j+1,k)-qt0(i,j,k))*geomh_invDYMv_8(j) )
            end do
         end do
      end do
!!$omp enddo nowait


!     ---------------------------------------------------------------
!
      return
      end subroutine SW_bac
