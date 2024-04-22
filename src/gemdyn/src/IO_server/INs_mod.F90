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

module INs
   use iso_c_binding
   use, intrinsic :: iso_fortran_env
   use vGrid_Descriptors
   implicit none
   public
   save

      character(len=2048) :: Path_input_S
      character(len=64), dimension(:), allocatable :: INs_list_S
      integer :: INs_GEM_COMM,INs_me1o1,INs_gem1o1,INs_host
      integer :: GZ_win,NEST_win,INs_hostmyproc
      logical :: INs_server_L=.false. , INs_1o1_L
      integer :: INs_nia,INs_nja,INs_nka,INs_nid,INs_njd,INs_hord
       
      real, dimension(:,:,:), pointer :: GZ
      real, dimension(:    ), pointer :: ND
      real, dimension(:,:  ), allocatable :: GZbuf, NDbuf
      integer, dimension(:), allocatable :: DIP1, GZIP1, &
                             iBUF, INs_isend
      character(len=32), dimension(:), allocatable :: cBUF
      integer :: GZcaract(5),INs_maxreqs,INs_maxNKA 
      real(kind=REAL64), pointer :: vtbl_8(:,:,:), VGD_tbl_8(:)
      
      character(len=16) :: Inp_datev
      logical Inp_src_hauteur_L, SRC_GZ_L, INs_vgd_L
      
      integer Inp_nfiles, Inp_kind, Inp_vgdkind, Inp_handle, &
              Inp_cmcdate, INs_nreq, Inp_rtag, INs_nplans, INs_n123(3)
      integer, dimension(:), contiguous,pointer :: Inp_list_unf => null()
      type(vgrid_descriptor) :: Inp_vgd_src

      type :: REQ
         character(len=32) :: vname(2)
         character(len=4 ) :: stag
         character(len=1 ) :: src
         integer :: unf,nk,deb
      end type REQ
      type(REQ),dimension (:), allocatable :: SRL
      
      type :: STD
         character(len=1) :: type
         integer :: fnt
      end type STD
      type(STD) :: fst(10)
            
end module INs
