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
!
!	Now we have a new BVP because the metric terms in the continuity equation
!   are no longer simplified to ln(Jz). See notes for more detail on
!   exactly which part changed
!
      subroutine pre ( F_dt_8, i0, j0, k0, in, jn, k0t )
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
      use step_options
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8
      
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

      w3 = grav_8 * tau_8

      if (Schm_opentop_L) then
!!$omp do
         do j= j0, jn
            do i= i0, in
               rhsb(i,j) = rhst(i,j,k0t)-mu_8*tau_nh_8* rhsw(i,j,k0t)
            end do
         end do
!!$omp enddo
      end if

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
!if(step_kount==2) print*,'rhsu:==>',i,j,k,rhsu(i,j,k)
!if(step_kount==2) print*,'rhsv:==>',i,j,k,rhsv(i,j,k)
!if(step_kount==2) print*,'rhsw:==>',i,j,k,rhsw(i,j,k)
!if(step_kount==2) print*,'rhst:==>',i,j,k,rhst(i,j,k)
!if(step_kount==2) print*,'rhsc:==>',i,j,k,rhsc(i,j,k)
!if(step_kount==2) print*,'rhsf:==>',i,j,k,rhsf(i,j,k)
               div  = (rhsu(i,j,k)-rhsu(i-1,j,k))*geomh_invDXM_8(j) &
                    + (rhsv(i,j,k)*geomh_cyM_8(j)-rhsv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + half*(GVM%mc_Ix_8(i,j,k)*(rhsu(i,j,k)+rhsu(i-1,j,k)) &
                           +GVM%mc_Iy_8(i,j,k)*(rhsv(i,j,k)+rhsv(i,j-1,k)) )
               rhsc(i,j,k) = div - invT_m_8*rhsc(i,j,k) 

            end do
         end do
      end do
!!$omp enddo


!!$omp do collapse(2)
      do k=k0t, l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in
!              Compute Rf'
!              ~~~~~~~~~~~
               w0 = invT_nh_8*(GVM%ztht_8(i,j,k)-Ver_z_8%t(k))
               rhsf(i,j,k) = rhsf(i,j,k) - w0


!              Compute Rt'
!              ~~~~~~~~~~~
               rhst(i,j,k) = gama_bdf_8 * ( w3 * rhst(i,j,k) + rhsw(i,j,k) )
            end do
         end do
      end do
!!$omp enddo

!     Compute Rc"', combining Rc" and Rt"
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$omp do collapse(2)
      do k= max(k0,2), l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in
               w1= (Ver_idz_8%m(k) + (GVM%mc_Iz_8(i,j,k)- epsi_8)*Ver_wp_8%m(k))
               w2= (Ver_idz_8%m(k) - (GVM%mc_Iz_8(i,j,k)- epsi_8)*Ver_wm_8%m(k))*Ver_onezero(k)

               rhsc(i,j,k) = rhsc(i,j,k) + invT_m_8* &
                            ((rhsf(i,j,k)-Ver_onezero(k)*rhsf(i,j,k-1))*Ver_idz_8%m(k) +&
                             GVM%mc_Iz_8(i,j,k)* &
                             (Ver_wp_8%m(k)*rhsf(i,j,k)+Ver_wm_8%m(k)*rhsf(i,j,k-1)*Ver_onezero(k))) &
                           +(w1 * rhst(i,j,k) - w2 * rhst(i,j,k-1))
            end do
         end do
      end do
!!$omp enddo
      if (k0 == 1) then
!!$omp do
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in
               w1= (Ver_idz_8%m(1) + (GVM%mc_Iz_8(i,j,1) - epsi_8)*Ver_wp_8%m(1))
               w2= (Ver_idz_8%m(1) - (GVM%mc_Iz_8(i,j,1) - epsi_8)*Ver_wm_8%m(1))*Ver_onezero(1)
               rhsc(i,j,1) = rhsc(i,j,1) + invT_m_8* &
                            ((rhsf(i,j,1)-Ver_onezero(1)*rhsf(i,j,1))*Ver_idz_8%m(1) +&
                             GVM%mc_Iz_8(i,j,1)* &
                             (Ver_wp_8%m(1)*rhsf(i,j,1)+Ver_wm_8%m(1)*rhsf(i,j,1)*Ver_onezero(1))) &
                           +(w1 * rhst(i,j,1) - w2 * rhst(i,j,1))
            end do
         end do
!!$omp enddo
       endif

       if (Schm_opentop_L) then
!     Apply opentop boundary conditions
!!$omp do
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in
               rhsb(i,j) = rhsb(i,j) - (one/tau_8+mu_8*tau_nh_8*grav_8)*(nest_t(i,j,k0t)&
                                  -  Cstv_Tstr_8)/nest_t(i,j,k0t)
            end do
         end do
!!$omp enddo
      endif


!     Apply lower boundary conditions
!!$omp do
      do j= j0, jn
!DIR$ SIMD
         do i= i0, in
            rhsc(i,j,l_nk) = rhsc(i,j,l_nk) - isol_d*GVM%mc_cssp_H_8(i,j) * &
                    (rhst(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhst(i,j,l_nk-1) &
                  + invT_m_8*(rhsf(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhsf(i,j,l_nk-1)))
         end do
      end do
!!$omp enddo

1000  format(3X,'UPDATE  THE RIGHT-HAND-SIDES: (S/R PRE_H)')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine pre 
