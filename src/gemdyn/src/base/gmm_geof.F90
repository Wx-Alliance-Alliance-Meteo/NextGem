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
module gmm_geof
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

!______________________________________________________________________
!                                                                      |
!  GMM VARIABLES ASSOCIATED WITH GEOPHYSICAL FIELDS (set_geof)         |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! fis0               | Phi srf at current timestep                     |
! topo_low           | Low resolution analysis orography before growth |
! topo_high          | High resolution target orography after growth   |
!-----------------------------------------------------------------------
!
!
      real, pointer, contiguous, dimension (:,:,:) :: orography => null()
      real, pointer, contiguous, dimension (:,:) :: fis0      => null()
      real, pointer, contiguous, dimension (:,:) :: fis0u     => null()
      real, pointer, contiguous, dimension (:,:) :: fis0v     => null()
      real, pointer, contiguous, dimension (:,:) :: orols     => null()
      real, pointer, contiguous, dimension (:,:) :: orolsu    => null()
      real, pointer, contiguous, dimension (:,:) :: orolsv    => null()
      real, pointer, contiguous, dimension (:,:,:) :: topo_low  => null()
      real, pointer, contiguous, dimension (:,:,:) :: topo_high => null()
      real, pointer, contiguous, dimension (:,:) :: me_full  => null()
      real, pointer, contiguous, dimension (:,:) :: me_large => null()
      real, pointer, contiguous, dimension (:,:,:) :: dgzm     => null()
      real, pointer, contiguous, dimension (:,:,:) :: dgzt     => null()
      real(kind=REAL64), allocatable, dimension (:,:,:) :: &
                            zthtu_8,zmomu_8, zthtv_8,zmomv_8
      
      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH), parameter :: gmmk_orography_s = 'OROGRAPHY'
      character(len=MAXNAMELENGTH), parameter :: gmmk_topo_low_s  = 'TOPOLOW'
      character(len=MAXNAMELENGTH), parameter :: gmmk_topo_high_s = 'TOPOHIGH'
      character(len=MAXNAMELENGTH), parameter :: gmmk_me_full_s   = 'MEFULL'
      character(len=MAXNAMELENGTH), parameter :: gmmk_me_large_s  = 'MELARGE'
      character(len=MAXNAMELENGTH), parameter :: gmmk_dgzm_s      = 'DGZM'
      character(len=MAXNAMELENGTH), parameter :: gmmk_dgzt_s      = 'DGZT'

end module gmm_geof
