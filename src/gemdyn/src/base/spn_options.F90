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
module spn_options
      use, intrinsic :: iso_fortran_env
      use :: transpose
      use :: gem_fft
      implicit none
   public
   save

   !# Nudging profile lower end in hyb level (eg. 1.0 or 0.8)
   !# profile will be set to 0. when hyb > 0.8
   real :: Spn_start_lev = 1.0
   namelist /spn  / Spn_start_lev
   namelist /spn_p/ Spn_start_lev

   !# Nudging profile upper end in hyb level (eg. 0.0 or 0.2)
   !# profile wll be set to 1.0 when hyb < 0.2
   real :: Spn_up_const_lev = 0.0
   namelist /spn  / Spn_up_const_lev
   namelist /spn_p/ Spn_up_const_lev

   !# Nudging profile transition shape('COS2' or 'LINEAR')
   !# Set the shape between Spn_start_lev and Spn_up_const_lev
   character(len=16) :: Spn_trans_shape_S = 'LINEAR'
   namelist /spn  / Spn_trans_shape_S
   namelist /spn_p/ Spn_trans_shape_S

   !# Nudging relaxation timescale - in hours
   real :: Spn_relax_hours = 10.
   namelist /spn  / Spn_relax_hours
   namelist /spn_p/ Spn_relax_hours

   !# The filter will be set to 0. for smaller scales (in km)
   real :: Spn_cutoff_scale_large = 300.
   namelist /spn  / Spn_cutoff_scale_large
   namelist /spn_p/ Spn_cutoff_scale_large

   !# The filter will be set to 1.0 for larger scales (in km)
   !# Transition between Spn_cutoff_scale_small and Spn_cutoff_scale_large
   !# will have a COS2 shape.
   real :: Spn_cutoff_scale_small = 100.
   namelist /spn  / Spn_cutoff_scale_small
   namelist /spn_p/ Spn_cutoff_scale_small

   !# Nudging interval - in sec
   !# Nudging is performed every every Spn_freq sec
   integer :: Spn_freq = -1
   namelist /spn  / Spn_freq
   namelist /spn_p/ Spn_freq

   !# Nudging weight in temporal space (.true. or .false.).
   !# If the driving fields are available every 6 hours and Spn_freq is
   !# set to 30 minutes then nudging will have more weight every six hours
   !# when the driving fields are available
   logical :: Spn_weight_L = .false.
   namelist /spn  / Spn_weight_L
   namelist /spn_p/ Spn_weight_L

   !# The weight factor when Spn_weight_L=.true.
   !# (The weigh factor is COS2**(Spn_wt_pwr), Spn_wt_pwr could  be set as
   !# 0, 2, 4, 6. If Spn_wt_pwr = 2, weight factor is COS2)
   integer :: Spn_wt_pwr = 2
   namelist /spn  / Spn_wt_pwr
   namelist /spn_p/ Spn_wt_pwr

   logical :: Spn_ON_L = .false.
   namelist /spn  / Spn_ON_L
   namelist /spn_p/ Spn_ON_L
      
   character(len=16) :: Spn_nudging_S = ' ' ! depricated

!    integer :: Spn_12smin, Spn_12smax ! Lower and upper bounds along x of the input array (before z->x transpose)
!    integer :: Spn_12sn               ! Number of local vertical levels (z) after the z->x transpose
!    integer :: Spn_22min , Spn_22max  ! Local lower and upper bounds of x after x->y transpose
!    integer :: Spn_22n , Spn_22n0     ! Number of local x-levels after x->y transpose and x halo size
!    integer :: Spn_22pil_w, Spn_22pil_e ! Local west/east pilot region size
   integer :: Spn_interval, Spn_ws   ! Spectral nudging application interval (#steps) and ??weight??
!    integer :: Spn_njnh               ! Number of local points in y that are outside the y-halo region, at start and after z->x transpose
!    integer :: Spn_nk12               ! Number of local points in z after transposes (both transposes))
!    integer :: Spn_ni22               ! Number of local points in x after x->y transpose

   !! Pilot regions applicable to the local grid, equal to pil_[nsew] on the appropriate boundary
   !! and zero otherwise
   integer :: Spn_pil_w, Spn_pil_e, Spn_pil_s, Spn_pil_n 

   !! Working grid parameters
   !! Number of points in i (x), j (y), and k (z) for the z-global working grid (xyz order)
   integer :: Spn_zgrid_lnx, Spn_zgrid_lny, Spn_zgrid_lnz
   !! Lower bounds _with respect to the global interior grid_ for the z-global working grid
   integer :: Spn_zgrid_llbx, Spn_zgrid_llby, Spn_zgrid_llbz
   !! Number of points in y, z, and x for the x-global working grid (yzx order)
   integer :: Spn_xgrid_lny, Spn_xgrid_lnz, Spn_xgrid_lnx 
   !! Lower bounds for the x-global working grid
   integer :: Spn_xgrid_llby, Spn_xgrid_llbz, Spn_xgrid_llbx
   !! Number of points in z, x, and y for the y-global working grid (zxy order)
   integer :: Spn_ygrid_lnz, Spn_ygrid_lnx, Spn_ygrid_lny
   !! Lower bounds for the y-global working grid
   integer :: Spn_ygrid_llbz, Spn_ygrid_llbx, Spn_ygrid_llby


   real :: Spn_weight ! Global weight factor for spectral nudging
   real, dimension(:), allocatable :: prof ! Profile of vertical weights (0 <= w <= 1)
   real(kind=REAL64) , dimension(:,:  ), allocatable :: Spn_flt ! Filter coefficients in wavenumber space
!    real(kind=REAL64) , dimension(:,:,:), allocatable :: Spn_fft,& 
!                                                 Spn_fdg!, Spn_wrk
   !! Spn_wrk is the z-contiguous grid, a subset of the tracer grid that excludes the global piloting region
   real(kind=REAL64), dimension(:,:,:), contiguous, pointer :: Spn_wrk
   !! Pointers for the grids used in the transpose operators, aliasing the respective parts of the transpose descriptors
   real(kind=REAL64), dimension(:,:,:), contiguous, pointer :: Spn_xgrid, & ! x-contiguous grid, after the first transpose
                                                               Spn_ygrid    ! y-contiguous grid, after the second transpose

   type(transpose_descriptor) :: zx_transpose, xy_transpose

   type(dft_descriptor) :: fft_x_forward, fft_x_reverse, fft_y_forward, fft_y_reverse

contains

!**s/r spn_nml - Read namelist spn

      integer function spn_nml (F_unf)
      use lun
      implicit none

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         spn_nml= 0
         if ( Lun_out >= 0) then
            if ( F_unf == -1 ) write (Lun_out,nml=spn_p)
            if ( F_unf == -2 ) write (Lun_out,nml=spn)
         end if
         return
      end if

      spn_nml= -1 ; nml_must= .false. ; nml_S= 'spn'

      rewind(F_unf)
      read (F_unf, nml=spn, end= 1001, err=1003)
      spn_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         spn_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (spn_nml < 0 ) return
      if ((Lun_out>=0).and.(spn_nml==0)) write (Lun_out, 6004) trim(nml_S)
      spn_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function spn_nml

end module spn_options
