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

      subroutine SW_bac3rd ( F_dt_8, i0, j0, k0, in, jn ,k0t )
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
      use yyg_param
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, km
      integer :: HLT_np, HLT_start, HLT_end, itpc, local_np
      real :: w5,w6
      real(kind=REAL64) :: tau_m_8, tau_nh_8, invT_m_8, invT_nh_8, Buoy, dqdx, dqdy
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) :: F_q

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
               F_q(i,j,k) = sngl(Sol_lhs(i,j,l_nk))
            end do
         end do
      end do
!!$omp enddo

      if ( Grd_yinyang_L) then
         call yyg_xchng_8 (F_q, YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, l_nk, .false., 'CUBIC', .true.)
      else
        !call HLT_split (F_k0, F_kn, HLT_np, HLT_start, HLT_end)
         call HLT_split (1, G_nk, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo_8 ( F_q(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
            km=max(k-1,1)
!DIR$ SIMD
            do i= i0, l_niu-pil_e

               dqdx = Hderiv8(F_q(i-1,j,k), F_q(i,j,k), &
                             F_q(i+1,j,k), F_q(i+2,j,k), geomh_invDX_8(j))
   !           Compute U
   !           ~~~~~~~~~
               ut0(i,j,k) = tau_m_8*(Ruu(i,j,k) - grav_8*dqdx)
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

               dqdy = Hderiv8(F_q(i,j-1,k),F_q(i  ,j,k),&
                             F_q(i,j+1,k),F_q(i,j+2,k),geomh_invDY_8)
   !           Compute V
   !           ~~~~~~~~~
               vt0(i,j,k) = tau_m_8*(Rvv(i,j,k) - grav_8*dqdy)
            end do
         end do
      end do
!!$omp enddo nowait


!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine SW_bac3rd
