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
!**********************************************************************
!
      subroutine nli ( F_dt_8, i0, j0, k0, in, jn, k0t )
      use HORgrid_options
      use gem_options
      use dyn_fisl_options
      use dynkernel_options
      use coriolis
      use geomh
      use tdpack
      use mem_tstp
      use gmm_geof
      use gmm_vt0
      use gmm_vt1
      use glb_ld
      use cstv
      use dcst
      use ver
      use metric
      use adz_mem
      use sol_mem
      use step_options
      use, intrinsic :: iso_fortran_env
      use stat_mpi, only:statf_dm
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, km 
      integer :: HLT_np, HLT_start, HLT_end
      real, dimension(:,:,:), pointer :: tots, logT, logQ, tots2, logT2, logQ2
      real(kind=REAL64) :: c0,div,w1,w2,w3,w4,barz,barzp,t_interp,u_interp,v_interp
      real(kind=REAL64) :: tau_8, tau_m_8, tau_nh_8, invT_8, invT_m_8, invT_nh_8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!     __________________________________________________________________
!
      tau_8    = (2.0*F_dt_8) / 3.0
      tau_m_8  = tau_8
      tau_nh_8 = tau_8
      invT_8   = one/tau_8
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

      c0 = Dcst_rayt_8**2
      w3 = grav_8 * tau_8

      tots (1:l_ni,1:l_nj,1:l_nk) => WS1(1:)
      logT (1:l_ni,1:l_nj,1:l_nk) => WS1(  l_ni*l_nj*l_nk+1:)

!!$omp do collapse(2)
      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               tots(i,j,k)= tt0(i,j,k)/Cstv_Tstr_8
               logT(i,j,k)= log(tots(i,j,k))
            end do
         end do
      end do
!!$omp enddo

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

!----NL_u
               barz  = Ver_wp_8%m(k)*tt0(i  ,j,k)+Ver_wm_8%m(k)*tt0(i  ,j,km)
               barzp = Ver_wp_8%m(k)*tt0(i+1,j,k)+Ver_wm_8%m(k)*tt0(i+1,j,km)
               t_interp = (barz+barzp)*half/Cstv_Tstr_8

               v_interp = 0.25d0*(vt0(i,j,k)+vt0(i,j-1,k)+vt0(i+1,j,k)+vt0(i+1,j-1,k))
               nl_u(i,j,k) = (t_interp-one) * ( qt0(i+1,j,k) - qt0(i,j,k) ) * geomh_invDX_8(j) &
                           - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v_interp &

   !           Adding vertical coordinate metric terms
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                           - (t_interp - one) * GVM%mc_Jx_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i+1,j,k+1)-qt0(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )   &
                                                 +(qt0(i  ,j,k+1)-qt0(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) ) &
                            +Ver_wm_8%m(k)*half*( (qt0(i+1,j,k  )-qt0(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)   &
                                                 +(qt0(i  ,j,k  )-qt0(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km) ) )

!----NL_v
               barz  = Ver_wp_8%m(k)*tt0(i,j  ,k)+Ver_wm_8%m(k)*tt0(i,j  ,km)
               barzp = Ver_wp_8%m(k)*tt0(i,j+1,k)+Ver_wm_8%m(k)*tt0(i,j+1,km)
               t_interp = (barz+barzp)*half/Cstv_Tstr_8

               u_interp = 0.25d0*(ut0(i,j,k)+ut0(i-1,j,k)+ut0(i,j+1,k)+ut0(i-1,j+1,k))

               nl_v(i,j,k) = (t_interp-one) * ( qt0(i,j+1,k) - qt0(i,j,k) ) * geomh_invDY_8 &
                           + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp &

   !           Adding vertical coordinate metric terms
   !           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                           - (t_interp - one) * GVM%mc_Jy_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i,j+1,k+1)-qt0(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )   &
                                                 +(qt0(i,j  ,k+1)-qt0(i,j  ,k ))*GVM%mc_iJz_8(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (qt0(i,j+1,k  )-qt0(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)   &
                                                 +(qt0(i,j  ,k  )-qt0(i,j  ,km))*GVM%mc_iJz_8(i,j  ,km) ) )
!----NL_q

               nl_c(i,j,k) = 0.d0
            end do
         end do
      end do
!!$omp enddo

! Set lateral boundaries for LAM configurations
      if (.not.Grd_yinyang_L) then
!!$omp do
         do k=1, l_nk
            if (l_west ) nl_u(     pil_w,1:l_nj,k) = 0.
            if (l_east ) nl_u(l_ni-pil_e,1:l_nj,k) = 0.
            if (l_south) nl_u(1:l_ni,     pil_s,k) = 0.
            if (l_south) nl_v(1:l_ni,     pil_s,k) = 0.
            if (l_north) nl_v(1:l_ni,l_nj-pil_n,k) = 0.
            if (l_west ) nl_v(     pil_w,1:l_nj,k) = 0.
         end do
!!$omp enddo nowait
      end if

!!$omp do collapse(2)
      do k=k0t,l_nk
         do j= j0, jn
            do i= i0, in
   !           Compute Nw
   !           ~~~~~~~~~~
               nl_w(i,j,k) = (tots(i,j,k)-one)*GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) &
                           -  (tots(i,j,k)-one)*grav_8*(one-one/tots(i,j,k)) 
   !           Compute Nt
   !           ~~~~~~~~~~
               w1 = invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
               nl_t(i,j,k) = gama_bdf_8 * ( w3 * w1 + nl_w(i,j,k) )
               true_nlt(i,j,k) = w1
            end do
         end do
      end do
!!$omp enddo

!     Compute Nc
!     ~~~~~~~~~~~
!!$omp do collapse(2)
      do k=k0,l_nk
         do j = j0, jn
            km=max(k-1,1)
            do i = i0, in
               w1= (Ver_idz_8%m(k) + (GVM%mc_Iz_8(i,j,k) - epsi_8)*Ver_wp_8%m(k))
               w2= (Ver_idz_8%m(k) - (GVM%mc_Iz_8(i,j,k) - epsi_8)*Ver_wm_8%m(k))*Ver_onezero(k)

               div = (nl_u(i,j,k)-nl_u(i-1,j,k)) * geomh_invDXM_8(j) &
                   + (nl_v(i,j,k)*geomh_cyM_8(j)-nl_v(i,j-1,k)* &
                      geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                    + (half * ( GVM%mc_Ix_8(i,j,k)*(nl_u(i,j,k)+nl_u(i-1,j,k)) &
                                     + GVM%mc_Iy_8(i,j,k)*(nl_v(i,j,k)+nl_v(i,j-1,k)) ))

               nl_c(i,j,k) = div - nl_c(i,j,k)*invT_m_8 !+ (w1 * nl_t(i,j,k) - w2 * nl_t(i,j,km))

               nl_c(i,j,k) = nl_c(i,j,k) + (w1 * nl_t(i,j,k) - w2 * nl_t(i,j,km))

               Sol_rhs(i,j,k) =  c0 * ( rhsc(i,j,k) - nl_c(i,j,k) )

            end do
         end do
      end do
!!$omp enddo

!     Apply top and bottom boundary conditions
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$omp do
      do j= j0, jn
         do i= i0, in
           ! Sol_rhs(i,j,l_nk) = Sol_rhs(i,j,l_nk) + isol_d * c0 * GVM%mc_cssp_H_8(i,j) * &
           !              (nl_t(i,j,l_nk ) - Ver_wmstar_8(G_nk)*nl_t(i,j,l_nk-1))
         end do
      end do
!!$omp enddo

      if (Schm_opentop_L) then
!!$omp do
         do j= j0, jn
            do i= i0, in
               w1= invT_8*( logT(i,j,k0t) - (one-one/tots(i,j,k0t)) )
               nl_b(i,j) = w1-mu_8*tau_nh_8* nl_w(i,j,k0t)
            end do
         end do
!!$omp enddo

      endif

      call vert_boundary ( i0,j0,in,jn )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine nli
