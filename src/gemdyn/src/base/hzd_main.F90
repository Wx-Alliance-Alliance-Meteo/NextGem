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

!**s/r hzd_main - Main controler for horizontal diffusion

      subroutine hzd_main
      use glb_ld
      use gmm_vt1
      use ens_gmm_var
      use ens_options
      use HORgrid_options
      use hvdif_options
      use ens_options
      use lun
      use tr3d
      use mem_tstp
      use mem_tracers
      use omp_timing
      use, intrinsic :: iso_fortran_env
      implicit none

      logical switch_on_UVW, switch_on_TR, switch_on_vrtspng, &
              switch_on_eqspng, switch_on_THETA
      logical xch_UV,xch_TT,xch_TR,xch_WZD
      integer i,n
!
!-------------------------------------------------------------------
!
      if (Lun_debug_L) write (Lun_out,1000)

      call gtmg_start (60, 'HZD_main', 1 )

      xch_UV = .false.
      xch_TT = .false.
      xch_TR = .false.
      xch_WZD= .false.
      if (ens_conf) then
          xch_UV = .true.
          xch_TT = .true.
         ! Save ut1,vt1 for SKEB (ens_markov_main)
         difut1 = ut1
         difvt1 = vt1
      end if
      switch_on_UVW     = Hzd_lnr       > 0.
      switch_on_TR      =(Hzd_lnr_tr    > 0.) .and. any(Tr3d_hzd)
      switch_on_THETA   = Hzd_lnr_theta > 0.
      switch_on_vrtspng =(Vspng_nk      >=1 ) .and. (Vspng_niter>0)
      switch_on_eqspng  = Eq_nlev       > 1

!**********************************
!  Horizontal diffusion on theta  *
!**********************************

      if ( switch_on_THETA ) then
         call gtmg_start (61, 'HZD_theta', 60)
         xch_TT = .true.
         call hzd_theta ()
         call gtmg_stop (61)
      end if

!**********************************
!  Horizontal diffusion on tracers*hzd_theta.F90
!**********************************

      if ( switch_on_TR ) then
         call gtmg_start (62, 'HZD_tracers', 60)
         xch_TR = .true.
         do i=1, Tr3d_ntr
            if (Tr3d_hzd(i)) then
               call hzd_exp_deln (tracers_P(i)%pntr, Hzd_pwr_tr,&
                     Hzd_lnR_tr, WS1, l_minx,l_maxx,l_miny,l_maxy,G_nk)
            end if
         end do
         call gtmg_stop (62)
      end if

!************************
!  Horizontal diffusion *
!************************

      if ( switch_on_UVW ) then
         call gtmg_start (64, 'HZD_bkgrnd', 60)
         xch_UV = .true.
         xch_TT = .true.
         xch_WZD= .true.
         call hzd_exp_deln ( wt1, Hzd_pwr, Hzd_lnR, WS1,&
                             l_minx,l_maxx,l_miny,l_maxy,4*G_nk)
         call gtmg_stop (64)
      end if

!********************
!  Vertical sponge  *
!********************

      xch_UV = xch_UV  .or. switch_on_vrtspng .or. switch_on_eqspng
      xch_TT = xch_TT  .or. switch_on_vrtspng
      xch_WZD= xch_WZD .or. switch_on_vrtspng
      if ( switch_on_vrtspng ) then
         call gtmg_start (65, 'V_SPNG', 60)
         call hzd_lid_sponge ()
         call gtmg_stop (65)
      end if

!**********************
!  Equatorial sponge  *
!**********************

      if ( switch_on_eqspng ) then
         call gtmg_start (67, 'EQUA_SPNG', 60)
         call eqspng (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call gtmg_stop (67)
      end if

      if (Ens_conf) call itf_ens_hzd () 

      call hzd_smago_main()

      call gtmg_stop (60)

 1000 format(3X,'MAIN HORIZONTAL DIFFUSION : (S/R HZD_MAIN)')
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_main
