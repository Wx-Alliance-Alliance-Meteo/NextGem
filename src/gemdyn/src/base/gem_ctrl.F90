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

!**s/r gem_ctrl - initiate the forward integration of the model

      subroutine gem_ctrl ()
      use, intrinsic :: iso_fortran_env
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use HORgrid_options
      use lam_options
      use step_options
      use glb_ld
      use gmm_geof
      use svri_mod
      use mem_nest
      use lun
      use rstr
      use omp_timing
      use tr3d
      use gmm_vt1
      implicit none

!object
!     Beginning of the integration. This subroutine
!     reads the data and performs initialization if required.
!     It then initiates the forward intergration of the model.

      character(len=16) :: Next_pilot_S
      logical :: rstrt_L= .false.
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: one=1.0d0, &
      sid=86400.0d0, rsid=one/sid
!     
!     ---------------------------------------------------------------
!
      call gemtime ( Lun_out, 'GEM_CTRL: START', .false. )

      call final_setup ()

      if ( .not. Rstri_rstn_L ) then
         call indata()
         if ( .not. Grd_yinyang_L .and. .not. Lam_ctebcs_L ) then
            dayfrac = Step_nesdt*rsid
            call incdatsd (Next_pilot_S, Step_runstrt_S, dayfrac)
            call itf_Iserv_request (Next_pilot_S,INs_list_S,1001,INs_nrequests)
         endif
      else
!!$omp parallel
         call vertical_metric ()
!!$omp end parallel
         if ( .not. Grd_yinyang_L .and. Lam_ctebcs_L ) then
            call nest_indata  (nest_u, nest_v , nest_w, nest_t   ,&
                               nest_q, nest_zd, nest_tr  ,&
                               nest_fullme,.false.,Step_runstrt_S,&
                               l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)
         endif
      endif
      call glbstat ( fis0,'ME',"indata",l_minx,l_maxx,l_miny,l_maxy,1,1,&
                     1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,1 )
      if (Schm_sleve_L) &
      call glbstat ( orols,'MELS',"indata",l_minx,l_maxx,l_miny,l_maxy, &
                 1,1, 1-G_halox,G_ni+G_halox,1-G_haloy,G_nj+G_haloy,1,1 )
                   
      call gemtime ( Lun_out, 'GEM_CTRL: INIT COMPLETED', .false. )
      call gtmg_stop ( 2 )

!      if (  Init_mode_L ) call initial (rstrt_L)

      if ( .not.rstrt_L ) call gem_run (rstrt_L)

      if (Lun_out > 0) write(Lun_out,3000) Lctl_step

 3000 format(/,'GEM_CTRL: END OF CURRENT TIME SLICE AT TIMESTEP',I8, &
             /,'===================================================')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine gem_ctrl
