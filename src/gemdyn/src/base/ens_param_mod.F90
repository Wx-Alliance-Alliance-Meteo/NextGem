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

module ens_param
   use wb_itf_mod
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      logical :: write_markov_l
      integer :: Ens_skeb_l,Ens_skeb_m
      integer :: itstep_s, iperiod_iau
      
      integer  :: lmin, lmax, nlon, nlat
      real :: fact3, cpdi, deltax
      real(kind=REAL64)  :: eps, dt
      real, dimension(:,:,:), allocatable  :: dsp_local, dsp_dif, dsp_gwd, psi_local
      real, dimension(:,:),   allocatable  :: fg1
      real, dimension(:),     allocatable  :: fg,fg2
      real, dimension(:),     allocatable  :: fact2
      real, dimension(:,:), pointer  :: f2
      real, dimension(:,:,:), pointer  :: f1,f1_str
      real, dimension(:,:), pointer  :: markov2
      real, dimension(:,:,:), pointer  :: markov1
      real(kind=REAL64), dimension(:,:), allocatable  :: pspec1
      real(kind=REAL64), dimension(:), allocatable  :: pspec2
      real(kind=REAL64), dimension(:), pointer  :: wrk2
      real(kind=REAL64), dimension(:,:), pointer  :: wrk1
      real(kind=REAL64), dimension(:,:,:), pointer  :: cc2
      real(kind=REAL64), dimension(:,:,:,:), pointer  :: cc1
     
      character(len=WB_MAXNAMELENGTH), dimension(:), allocatable :: vname
      integer, dimension(:), allocatable :: vlmin, vlmax, vnlon, vnlat
      real, dimension(:), allocatable :: vfmin, vfmax, vfstd, vfstr, vtau, veps
      real,    dimension(:,:), pointer, save :: tropwt=>null()

end module ens_param
