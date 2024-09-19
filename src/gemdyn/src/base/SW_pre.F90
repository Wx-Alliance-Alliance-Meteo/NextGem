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
!
!**s/r pre_H -  Combine some rhs obtaining Rt' and Rc", the linear
!               contributions to the rhs of Helmholtz equation
!            -  Height-type vertical coordinate
!
!	Now we have a new BVP because the metric terms in the continuity equation
!   are no longer simplified to ln(Jz). See notes for more detail on
!   exactly which part changed
!
      subroutine SW_pre ( F_dt_8, i0, j0, k0, in, jn, k0t )
      use HORgrid_options
      use gem_options
      use geomh
      use gmm_geof
      use mem_tstp
      use mem_nest
      use tdpack
      use glb_ld
      use lun
      use cstv
      use ver
      use sol_mem
      use dyn_fisl_options
      use dynkernel_options
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
      HLT_np, HLT_start, HLT_end

      integer :: i, j, k, km
      real    :: w0
      real(kind=REAL64) :: div, w1, w2, w3, tau_8, tau_m_8, tau_nh_8
      real(kind=REAL64) :: avg_ludz, avg_lvdz !for the new blocks
      real(kind=REAL64) :: invT_8, invT_m_8, invT_nh_8
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      tau_8    = (2.d0 * F_dt_8 ) / 3.d0
      tau_m_8  = tau_8
      tau_nh_8 = tau_8
      invT_8   = one/tau_8
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

      w3 = grav_8 * tau_8/ Cstv_Tstr_8

!*************************************
! Combination of governing equations *
!*************************************

!     Compute Rc' combining Ru, Rv and Rc
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in
               div  = (rhsu(i,j,k)-rhsu(i-1,j,k))*geomh_invDXM_8(j) &
                    + (rhsv(i,j,k)*geomh_cyM_8(j)-rhsv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) 
               rhsc(i,j,k) = div/grav_8 - (invT_8*Cstv_h0inv_8/grav_8)*rhsc(i,j,k)
            end do
         end do
      end do
!!$omp enddo

1000  format(3X,'UPDATE  THE RIGHT-HAND-SIDES: (S/R PRE_H)')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_pre
