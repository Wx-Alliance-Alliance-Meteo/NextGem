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
!**s/r nli_H  - compute non-linear terms:  Nu, Nv, Nt, Nc, Nw,
!             - compute full right-hand side of Helmholtz eqn: Rp=Rc-Nc
!             - Height-type vertical coordinate
!
!	combine nonlinear terms
!**********************************************************************
!
      subroutine SW_comb_nli ( F_dt_8, i0, j0, k0, in, jn, k0t, itpc)

      use HORgrid_options
      use gem_options
      use step_options
      use dyn_fisl_options
      use dynkernel_options
      use coriolis
      use geomh
      use tdpack
      use mem_tstp
      use gmm_geof
      use gmm_vt1
      use glb_ld
      use cstv
      use dcst
      use ver
      use metric
      use sol_mem
      use mem_tstp
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8
      integer, intent(IN)  :: itpc 

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
            HLT_np, HLT_start, HLT_end

      integer :: i, j, k, km
      real(kind=REAL64) :: c0,div,w1,w2,w3,tots, i_hydro
      real(kind=REAL64) :: avg_ludz, avg_lvdz 
      real(kind=REAL64) :: tau_8, tau_m_8, tau_nh_8, invT_8, invT_m_8, invT_nh_8, tmp_nlc
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!     __________________________________________________________________
!
      i_hydro=0.d0

      tau_8    = (2.0*F_dt_8) / 3.0
      tau_m_8  = tau_8
      tau_nh_8 = tau_8
      invT_8   = one/tau_8
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

      !c0 = Dcst_rayt_8**2
      c0 = 1.d0
      w3 = grav_8 * tau_8 / Cstv_Tstr_8
      !print *,"c0 = ", c0

!***********************************************************
! The nonlinear deviation of horizontal momentum equations *
!***********************************************************

!     Compute Nu, Nv and Nc
!     ~~~~~~~~~~~~~~~~~~~~~

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0-1, jn
            km=max(k-1,1)
            do i= i0-1, in

               !---nonlinear terms for u
               nl_u(i,j,k) = nlu_t1(i,j,k)

               !---nonlinear terms for v---
               nl_v(i,j,k) = nlv_t1(i,j,k) 

               !---nonlinear terms for c 
               nl_c(i,j,k) = nlq_t1(i,j,k) 

            end do
         end do
      end do
!!$omp enddo


!     Compute Nc
!     ~~~~~~~~~~~
!NOTE: how Nc is computed is different from original
!only keep terms with the epsilon coefficient
!DO NOT  need \delta_Z N_T

   !print *,"i,j,k        rhs sol    |    rhsc     |    nl_c   "
   !print *,"i,j,k    nlc (b4)  |    div   |   w1*nlt(k)   |   w2*nlt(km)   |   nlc (a4)  |  "

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( nl_u(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( nl_v(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( nl_c(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)


!!$omp do collapse(2)
      do k=k0,l_nk
         do j = j0, jn
            do i = i0, in
               div = (nl_u(i,j,k)-nl_u(i-1,j,k)) * geomh_invDXM_8(j) &
                   + (nl_v(i,j,k)*geomh_cyM_8(j)-nl_v(i,j-1,k)* &
                      geomh_cyM_8(j-1))*geomh_invDYM_8(j) 

               nl_c(i,j,k) = div/grav_8 - invT_8*((1.d0-Cstv_swln_8)*Cstv_h0inv_8+Cstv_swln_8)/grav_8*nl_c(i,j,k)

               Sol_rhs(i,j,k) = c0*(rhsc(i,j,k)-nl_c(i,j,k))

            end do
         end do
      end do
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_comb_nli
