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
module dyn_fisl_options
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   logical :: Euler_step_one = .true.
      
   !# Add souce level l_nk+1 in qt SL interpolation & trajectories
   logical :: SL_sfc  = .false.
   namelist /dyn_fisl  / SL_sfc
   namelist /dyn_fisl_p/ SL_sfc

   !# T* basic state temperature (K)
   real(kind=REAL64) :: Cstv_Tstr_8 = 240.0
   namelist /dyn_fisl  / Cstv_Tstr_8
   namelist /dyn_fisl_p/ Cstv_Tstr_8

   !# Inverse of mean height (m^-1) for Shallow-water
   real(kind=REAL64) :: Cstv_h0inv_8 

   !# Fraction of adjustment to be given to the ocean
   real(kind=REAL64) :: Cstv_psadj_8 = 1.d0
   namelist /dyn_fisl  / Cstv_psadj_8
   namelist /dyn_fisl_p/ Cstv_psadj_8

   !Schm

   !# Precicion order for 3D spacial operators
   integer :: Schm_POSO = 2
   namelist /dyn_fisl  / Schm_POSO
   namelist /dyn_fisl_p/ Schm_POSO

   !# Maximum iterations to solve trajectories
   integer :: Schm_itSL  = 6
   namelist /dyn_fisl  / Schm_itSL
   namelist /dyn_fisl_p/ Schm_itSL

   !# Trajectory precicion stopping criteria
   !# set to -1. constant number of iterations required(Schm_itSL above)
   real :: Schm_tolSL  = 1e-6
   namelist /dyn_fisl  / Schm_tolSL
   namelist /dyn_fisl_p/ Schm_tolSL
      
   !# Trajectory convergence rate stopping criteria
   real :: Schm_rateSL  = .1
   namelist /dyn_fisl  / Schm_rateSL
   namelist /dyn_fisl_p/ Schm_rateSL
      
   !# Maximum Picard iterations to solve non-linear Helmholtz problem
   integer :: Schm_itpc  = 10
   namelist /dyn_fisl  / Schm_itpc
   namelist /dyn_fisl_p/ Schm_itpc

   !# Picard position stopping criteria
   real :: Schm_tolpic  = 1e-6
   namelist /dyn_fisl  / Schm_tolpic
   namelist /dyn_fisl_p/ Schm_tolpic
      
   !# Picard convergence rate stopping criteria
   real :: Schm_ratepic  = .1
   namelist /dyn_fisl  / Schm_ratepic
   namelist /dyn_fisl_p/ Schm_ratepic
      
   !# *  -1: no blending between Yin and Yang
   !# *   0: blending at init only
   !# * > 0: blending at every nblendyy timestep
   integer :: Schm_nblendyy = -1
   namelist /dyn_fisl  / Schm_nblendyy
   namelist /dyn_fisl_p/ Schm_nblendyy

   !# * 0 -> No conservation of surface pressure
   !# * 1 -> Conservation of Total air mass Pressure
   !# * 2 -> Conservation of Dry air mass Pressure
   integer :: Schm_psadj = 0
   namelist /dyn_fisl  / Schm_psadj
   namelist /dyn_fisl_p/ Schm_psadj

   !# True-> print dry/wet air masses
   logical :: Schm_psadj_print_L = .false.
   namelist /dyn_fisl  / Schm_psadj_print_L
   namelist /dyn_fisl_p/ Schm_psadj_print_L

   !# True-> Tracers are mixing ratios with respect to dry air mass
   logical :: Schm_dry_mixing_ratio_L = .false.
   namelist /dyn_fisl  / Schm_dry_mixing_ratio_L
   namelist /dyn_fisl_P/ Schm_dry_mixing_ratio_L

   !# True-> use SLEVE vertical coordinate
   logical :: Schm_sleve_L = .false.
   !# True-> force the reading of MELS from geophysical file
   logical :: Schm_orols_fromgeophy_L = .false.
   namelist /dyn_fisl  / Schm_orols_fromgeophy_L
   namelist /dyn_fisl_P/ Schm_orols_fromgeophy_L
   !# True-> will retain 50% of Schm_orols_rc deltax waves
   real :: Schm_orols_rc = 10.
   namelist /dyn_fisl  / Schm_orols_rc
   namelist /dyn_fisl_P/ Schm_orols_rc
   !# True-> Number of pass in the filter
   integer :: Schm_orols_np = 20
   namelist /dyn_fisl  / Schm_orols_np
   namelist /dyn_fisl_P/ Schm_orols_np
   !# Filter type for oro filtering to obtain large scale oro
   character(len=16) :: Schm_orols_ftype_S = 'MC2'
!   namelist /dyn_fisl  / Schm_orols_ftype_S
!   namelist /dyn_fisl_p/ Schm_orols_ftype_S

   !# True-> to use topography
   logical :: Schm_Topo_L = .true.
   namelist /dyn_fisl  / Schm_Topo_L
   namelist /dyn_fisl_p/ Schm_Topo_L

   !# Interpolation type: TURBO(cubic), 3(cubic), 5(quintic)
   character(len=16) :: Schm_advec_type_S = 'TURBO'
   namelist /dyn_fisl  / Schm_advec_type_S
   namelist /dyn_fisl_p/ Schm_advec_type_S
   
   !# * 0   ->          NO advection
   !# * 1   -> traditional advection
   !# * 2   -> consistent advection with respect to off-centering
   !# * 3   -> reversed advection with respect to off-centering
   integer :: Schm_advec = 1
   namelist /dyn_fisl  / Schm_advec
   namelist /dyn_fisl_p/ Schm_advec

   !# True-> Modify slightly code behaviour to ensure bitpattern
   !# reproduction in restart mode using FST file
   logical :: Schm_bitpattern_L = .false.
   namelist /dyn_fisl/ Schm_bitpattern_L

   !# Apply water loading in the calculations
   logical :: Schm_wload_L = .false.
   namelist /dyn_fisl  / Schm_wload_L
   namelist /dyn_fisl_p/ Schm_wload_L

   logical Schm_opentop_L
   integer Schm_nith

contains

!**s/r dyn_fisl_nml - Read namelist dyn_fisl

      integer function dyn_fisl_nml (F_unf)
      use HORgrid_options
      use gem_options
!      use gem_options
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <clib_interface_mu.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer :: err!,G_halox,G_haloy
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         dyn_fisl_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=dyn_fisl_p)
            if ( F_unf == -2 ) write (Lun_out,nml=dyn_fisl)
         end if
         return
      end if

      dyn_fisl_nml= -1 ; nml_must= .false. ; nml_S= 'dyn_fisl'

      rewind(F_unf)
      read (F_unf, nml=dyn_fisl, end= 1001, err=1003)
      dyn_fisl_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         dyn_fisl_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      else
         if (Lun_out >= 0) write (Lun_out, 6009) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (dyn_fisl_nml < 0 ) return
      if ((Lun_out>=0).and.(dyn_fisl_nml==0)) write (Lun_out, 6004) trim(nml_S)
      dyn_fisl_nml= 1

      err = clib_toupper(Schm_orols_ftype_S)
      
!     if ( .not. Grd_yinyang_L) Schm_POSO = 2
      G_halox=4
      if (Schm_POSO==3) G_halox=5
      if (Schm_POSO==5) G_halox=6
      G_haloy=G_halox
      
 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
 6009 format (//,' NAMELIST ',A,' IS MANDATORY'//)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function dyn_fisl_nml

   function dyn_fisl_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function dyn_fisl_options_init

end module dyn_fisl_options
