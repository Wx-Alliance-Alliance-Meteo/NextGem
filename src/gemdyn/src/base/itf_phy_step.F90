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

!**s/r itf_phy_step - Apply the physical processes: CMC/RPN package

      subroutine itf_phy_step ( F_step_kount, F_lctl_step )
      use iso_c_binding
      use ctrl
      use glb_ld
      use phy_itf
      use itf_phy_cloud_objects, only: cldobj_displace,cldobj_expand,CLDOBJ_OK
      use itf_phy_filter, only: ipf_smooth_fld, sfcflxfilt_o, nsurfag
      use ens_options
      use init_options
      use lun
      use tr3d
      use rstr
      use path
      use wb_itf_mod
      use ptopo
      use omp_timing
      use HORgrid_options
      use step_options
      use gmm_pw
      use gmm_phy, only: phy_cplm, phy_cplt
      use dyn_fisl_options
      implicit none

      integer, intent(in) :: F_step_kount, F_lctl_step

!arguments
!  Name                    Description
!----------------------------------------------------------------
! F_step_kount             step count
! F_lctl_step              step number
!----------------------------------------------------------------

      integer,external :: itf_phy_prefold_opr

      integer :: err_geom, err_step, err_input, err_smooth, err, k
      logical :: cloudobj,must_stop_L
      real, dimension(:,:,:), pointer :: zud, zvd
!
!     ---------------------------------------------------------------
!
      if ( .not. Ctrl_phyms_L ) return
      if ( Rstri_user_busper_L .and. (F_step_kount == 0) ) return

      call gtmg_start ( 40, 'PHYSTEP', 1 )

!!$omp single
      if (Lun_out > 0) write (Lun_out,1001) F_lctl_step

      if (F_step_kount == 0) then
         call itf_phy_geom (err_geom)
         if (NTR_Tr3d_ntr > 0) then
            err = wb_put('itf_phy/READ_TRACERS', &
                                NTR_Tr3d_name_S(1:NTR_Tr3d_ntr))
         end if
         phy_cplm(:,:) = 1.
         phy_cplt(:,:) = 1.
      end if

      !call pw_glbstat('PW_BEF')
      !if (F_step_kount == 0) call itf_phy_glbstat('befinp')
!!$omp end single
      
!!$omp do
      do k= 1, G_nk
         pw_uu_copy(:,:,k) = pw_uu_plus(:,:,k)
         pw_vv_copy(:,:,k) = pw_vv_plus(:,:,k)
      end do
!!$omp end do
      
      call gtmg_start ( 45, 'PHY_input', 40 )
!!$omp single
      err_input = phy_input1 ( itf_phy_prefold_opr, F_step_kount, &
            Path_phyincfg_S, Path_phy_S, 'GEOPHY/Gem_geophy.fst' )
!!$omp end single copyprivate(err_input)

      !if (F_step_kount == 0) call itf_phy_glbstat('Aftinp')
      call gem_error_omp (err_input,'itf_phy_step','Problem with phy_input')
      call gtmg_stop  ( 45 )

      ! Smooth the thermodynamic state variables on request
!!$omp single
      err_smooth = min(&
           ipf_smooth_fld('rt',   'rt_smt',   3, 1), &
           ipf_smooth_fld('rdpr', 'rdpr_smt', 3, 1), &
           ipf_smooth_fld('rdqi', 'rdqi_smt', 3, 1), &
           ipf_smooth_fld('tcond','tcond_smt',3) &
           )

      ! Call digital filter to smooth Alfa, Beta surface fields
      if (sfcflxfilt_o > 1 .and. F_step_kount > 0) then
         call itf_dfilter('alfat','falfat',1)
         call itf_dfilter('alfaq','falfaq',1)
         call itf_dfilter('bt','fbt',nsurfag)
      endif
      
      ! Advect cloud objects
      if (.not.WB_IS_OK(wb_get('phy/deep_cloudobj',cloudobj))) cloudobj = .false.
!!$omp end single copyprivate(err_smooth,cloudobj)

      call gem_error_omp (err_smooth,'itf_phy_step','Problem with ipf_smooth_fld')

      if (cloudobj .and. F_lctl_step > 0) then
!!$omp single
         must_stop_L= (cldobj_displace() /= CLDOBJ_OK)
!!$omp end single copyprivate(must_stop_L)
         if (must_stop_L) call gem_error_omp (-1,'itf_phy_step','Problem with cloud object displacement')
      endif

      call gtmg_start ( 46, 'PHY_step', 40 )
      err_step = phy_step ( F_step_kount, F_lctl_step )
      call gtmg_stop  ( 46 )

      call gem_error_omp (err_step,'itf_phy_step','Problem with phy_step')
      !call itf_phy_glbstat('aftphy')

      call gtmg_start ( 47, 'PHY_update', 40 )

      call itf_phy_update3 ( F_step_kount > 0 )

      call gtmg_stop  ( 47 )
      !call pw_glbstat('PW_AFT')

      call gtmg_start ( 48, 'PHY_output', 40 )
!!$omp single
      call itf_phy_output ( F_lctl_step )
!!$omp end single
      call gtmg_stop  ( 48 )

!!$omp single
      call itf_phy_diag ()
      
!!$      if ( Init_mode_L ) then
!!$         if (F_step_kount ==   Init_halfspan) then
!!$            err = phy_snapshot('W')
!!$         else if (F_step_kount == 2*Init_halfspan) then
!!$            if (.not. ens_skeb_tndfix) then
!!$               nullify(zud, zvd)
!!$               err = phy_get(zud, 'phytd_udis', F_npath='V', F_bpath='P')
!!$               err = min(err, phy_get(zvd, 'phytd_vdis', F_npath='V', F_bpath='P'))
!!$               call gem_error (err, 'itf_phy_step', 'Cannot retrieve dissipations')
!!$            endif
!!$            err = phy_snapshot('R')
!!$            if (.not. ens_skeb_tndfix) then
!!$               err = phy_put(zud, 'phytd_udis', F_npath='V', F_bpath='P')
!!$               err = min(err, phy_put(zvd, 'phytd_vdis', F_npath='V', F_bpath='P'))
!!$               call gem_error (err, 'itf_phy_step', 'Cannot reset dissipations')
!!$               deallocate(zud, zvd)
!!$            endif
!!$         end if
!!$      end if
!!$omp end single
      
      call gtmg_stop ( 40 )

 1001 format(/,'PHYSICS : PERFORMING TIMESTEP #',I9, &
             /,'========================================')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine itf_phy_step
