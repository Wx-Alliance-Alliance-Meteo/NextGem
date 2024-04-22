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
module itf_phy_mem
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save
      real, dimension(:,:  ), pointer :: delq
      real, dimension(:,:,:), pointer :: &
                                   qt1i,pw_uu_plus0,pw_vv_plus0,&
                                   pw_tt_plus0,qw_phy,qw_dyn,tdu,tdv,tv
      real(kind=REAL64), dimension(:,:,:), pointer :: pm_dyn_8
      real(kind=REAL64), dimension(:,:  ), pointer :: p0_8

end module itf_phy_mem

module itf_phy_io
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save
   integer :: Phy_comm_id, Phy_npes
end module itf_phy_io
