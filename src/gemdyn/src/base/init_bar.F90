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

!**s/r init_bar - prepare data for autobarotropic runs (Williamson cases)

      subroutine init_bar ()
      use cstv
      use ctrl
      use dyn_fisl_options
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt0
      use gmm_vt1
      use gmm_vt2
      use inp_mod
      use mem_tracers
      use step_options
      use tdpack
      use tr3d
      use ver
      use wil_options
      implicit none

      !object
      !============================================================
      !Prepare data for autobarotropic runs (Williamson cases)
      !------------------------------------------------------------
      !NOTE: U,V output on Staggered grids
      !============================================================

      integer :: k,i0,in,j0,jn
      integer :: HLT_start, HLT_end, local_np
      real, dimension (:,:,:), pointer :: hu
      real, dimension (l_minx:l_maxx,l_miny:l_maxy,G_nk) :: gz_t
      real(kind=REAL64) :: FI_8
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0
!
!-------------------------------------------------------------------
!
      if (Vtopo_L)      call gem_error (-1,'INIT_BAR','Vtopo_L not available YET')

      if (Schm_sleve_L) call gem_error (-1,'INIT_BAR','  SLEVE not available YET')

      orols = 0. ; topo_low = 0. ; topo_high = 0.

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      !Setup Williamson Case 7: The 21 December 1978 Initial conditions are read
      !-------------------------------------------------------------------------
      if (Williamson_case==7) then

         call inp_data ( ut1(l_minx,l_miny,1), vt1(l_minx,l_miny,1), &
                         wt1(l_minx,l_miny,1), tt1(l_minx,l_miny,1), &
                         qt1(l_minx,l_miny,1),zdt1(l_minx,l_miny,1), &
                         p0,trt1,fis0, orols, .true., Step_runstrt_S,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)

         !Initialize PHI perturbation in q
         !--------------------------------
         FI_8= grav_8*Ver_z_8%m(1)
         do k=1,G_nk+1
            qt1(i0:in,j0:jn,k) = p0(i0:in,j0:jn) - FI_8
         end do

      end if

      !Initialize T/ZD/W/Q
      !-------------------
      tt1 = Cstv_Tstr_8 ; zdt1 = 0. ; wt1 = 0.

      !Initialize HU
      !-------------
      hu => tracers_P(Tr3d_hu)%pntr

      hu = 0.

      !Initialize d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Prepare initial conditions (staggered u-v,gz,s,topo) for Williamson cases
      !-------------------------------------------------------------------------
      call wil_init (ut1(l_minx,l_miny,1),vt1(l_minx,l_miny,1),gz_t,p0,fis0,qt1(l_minx,l_miny,1),l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

      !Required for CASE5 LAM version (NESTING)
      !----------------------------------------
      topo_high(i0:in,j0:jn,1) =      fis0(i0:in,j0:jn)
      topo_low (i0:in,j0:jn,1) = topo_high(i0:in,j0:jn,1)
      topo_high(i0:in,j0:jn,2) =     orols(i0:in,j0:jn)
      topo_low (i0:in,j0:jn,2) = topo_high(i0:in,j0:jn,2)

      !Estimate U-V and T on scalar grids
      !----------------------------------
!!$omp parallel private (local_np, HLT_start, HLT_end)
      call HLT_split (-2, G_nk+3, local_np, HLT_start, HLT_end)
      call gem_xch_halo (ut1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
      call gem_xch_halo (vt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
      call hwnd_stag2( pw_uu_plus,pw_vv_plus,ut1(l_minx,l_miny,1),vt1(l_minx,l_miny,1),&
                       l_minx,l_maxx,l_miny,l_maxy,G_nk   ,&
                       1-G_halox*west ,l_niu+G_halox*east ,&
                       1-G_haloy*south,l_njv+G_haloy*north, .false. )
!!$omp end parallel

      pw_tt_plus(:,:,1:G_nk) = tt1(:,:,1:G_nk) 
!
!-------------------------------------------------------------------
!
      return
      end subroutine init_bar
