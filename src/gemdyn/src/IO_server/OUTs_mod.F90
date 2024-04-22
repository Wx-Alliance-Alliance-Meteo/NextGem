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

module OUTs
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      type :: sfc_var
         sequence
         character(len=3) :: stag
         character(len=4) :: nv
         integer :: knd, nbits, indx, k0, skip
         real :: lvl
      end type sfc_var
      type :: fst
         character(len=32  ) :: type
         character(len=1024) :: name
         integer :: fnt
      end type fst
      type(fst) :: out(10)
      
      character(len=1024) :: Out_dirname_S
      character(len=32  ) :: Out_filenames_S(20)
      character(len=12  ) :: Out_etik_S
      character(len=4   ) :: Out_nomvar
      character(len=1   ) :: Out_typvar_S,chac1
      logical :: OUTs_server_L=.false.
      logical :: OUTs_1o1_L, Out_reduc_L, Out_rewrit_L
      integer :: OUTs_GEM_COMM,OUTs_me1o1,OUTs_gem1o1
      integer :: clients,data_wm,IOS_events=0
      integer :: nout_files, Out_kind, Out_nbit, OUTs_hostmyproc
      integer :: Out_i0,Out_in,Out_j0,Out_jn
      integer :: Out_ip2,Out_ip3,Out_npas,Out_dateo,Out_deet
      integer :: Out_ig1,Out_ig2,Out_ig3,Out_ig4,Out_unf
      integer :: SHARED_WIN
      integer, dimension (:    ), allocatable :: metaG
      real, dimension(:), pointer :: levels

end module OUTs
