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

!**   s/r set_dync - initialize the dynamics model configuration

      subroutine set_dync (F_dt_8)
      use dynkernel_options
      use cstv
      use dcst
      use dyn_fisl_options
      use glb_ld
      use lam_options
      use tdpack
      use ver
      use metric
      use step_options
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN)  :: F_dt_8

      integer :: k0
      real(kind=REAL64)  :: w1, w2, Nstr2_8, cstr2_8, bdf_tau
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------
!
      Nstr2_8=grav_8*grav_8/(cpd_8*Cstv_Tstr_8)
      cstr2_8=rgasd_8*Cstv_Tstr_8/(one-cappa_8)
      mu_8=Nstr2_8/grav_8
      epsi_8=grav_8/cstr2_8

      !--bdf variables--
      bdf_tau    = (2.d0 * F_dt_8)/ 3.d0
      gama_bdf_8 = bdf_tau/(bdf_tau + Nstr2_8*bdf_tau**3)
      gg_bdf_8   = epsi_8/(grav_8*bdf_tau**2)
      gg_sw_bdf_8=1.d0/(bdf_tau**2)*Cstv_h0inv_8/grav_8
      !-----------------

      gama_8         = gama_bdf_8
      gg_8           = gg_bdf_8
      gg_sw_8        = gg_sw_bdf_8

      Cstv_tau_8     = bdf_tau
      Cstv_invT_8    = one/Cstv_tau_8
      Cstv_tau_m_8   = Cstv_tau_8
      Cstv_invT_m_8  = Cstv_invT_8
      Cstv_tau_nh_8  = Cstv_tau_8
      Cstv_invT_nh_8 = Cstv_invT_8

      Cstv_swln_8 = 0.d0
      if(Dynamics_swln_L) Cstv_swln_8 = 1.d0

      Ver_css_8   = one/gama_8 / (Ver_idz_8%t(G_nk)-mu_8*half)
      Ver_alfas_8 = Ver_css_8*gama_8* &
                    (Ver_idz_8%t(G_nk  )+mu_8*half &
                    +Ver_wmstar_8(G_nk)*(Ver_idz_8%t(G_nk-1)-mu_8*half) )
      Ver_betas_8 = Ver_css_8*gama_8* &
      Ver_wmstar_8(G_nk)*(Ver_idz_8%t(G_nk-1)+mu_8*half)

      Ver_alfat_8 = one
      Ver_cst_8   = zero
      Ver_cstp_8  = zero

      k0=1+Lam_gbpil_T
      if (Schm_opentop_L) then
         w1= Ver_idz_8%t(k0-1)*(Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8)
         w2= mu_8*half*(Ver_idz_8%m(k0)+Ver_wm_8%m(k0)*epsi_8)
         Ver_cst_8   =  one / (-(mu_8* Cstv_tau_nh_8)*Ver_idz_8%t(k0-1) + half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))
         Ver_alfat_8 =(-(mu_8* Cstv_tau_nh_8)*Ver_idz_8%t(k0-1) - half* one/(Cstv_tau_8*cpd_8*Cstv_Tstr_8))*Ver_cst_8
         Ver_cstp_8  = gama_8*(w1 + w2)*Ver_cst_8
      end if

!
!     ---------------------------------------------------------------
!
      return
      end
