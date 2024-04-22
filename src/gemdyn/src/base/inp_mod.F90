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
module inp_mod
      use iso_c_binding
      use inp_options
      use vGrid_Descriptors
      use, intrinsic :: iso_fortran_env
      implicit none
      public
      save

      character(len=1 ) :: Inp_levtype_S
      character(len=16) :: Inp_datev
      logical :: Inp_src_PX_L=.false. , Inp_src_GZ_L=.false.
      logical Inp_src_hauteur_L, Inp_zd_L, Inp_w_L, Inp_qt_L
      integer Inp_nfiles , Inp_kind ,&
              Inp_version, Inp_handle , Inp_cmcdate
      integer :: Inp_iope, Inp_bcast, Inp_gtmg(2)
      integer, dimension(:), contiguous,pointer :: Inp_list_unf => null()
      integer :: Inp_PX_kind, Inp_GZ_kind, Inp_PX_nka, Inp_GZ_nka
      integer :: Inp_comm=0, Inp_window=0
      integer, parameter :: Inp_maxNKA=200
      real, dimension (:), pointer :: Inp_recv
      type(vgrid_descriptor) :: Inp_vgd_src
      real(kind=REAL64) Inp_pref_a_8

      type Inp_vrt3d
         logical me_L, mels_L, sleve_L
         integer nk,kind
         integer, dimension(:), pointer :: ip1
         real, dimension(:,:,:), pointer :: valq,valu,valv,sfc
      end type Inp_vrt3d
      type(Inp_vrt3d) :: GZ3d, PX3d
      
contains

!**s/r inp_init

      subroutine inp_init ()
      use glb_ld
      use ptopo
      use gem_options
      use, intrinsic :: iso_fortran_env
      implicit none

      include 'mpif.h'
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
 !     integer(C_INTPTR_T) :: WINSIZE
      integer :: dim, DISP_UNIT, INFO, err
      integer i,col_init,row_init,col,row,inc_row,inc_col
      integer :: color,ordinal_in_set,new_comm,me
      type(C_PTR), save :: basepntr
      integer ORIGIN_DATATYPE, TARGET_DATATYPE, TARGET_RANK
      integer ORIGIN_COUNT, TARGET_COUNT
      integer(C_INTPTR_T) :: TARGET_DISP
!     
!-------------------------------------------------------------------
!
      Inp_npes = max(1,min(Inp_npes ,Ptopo_npex*Ptopo_npey))
      col_init=0 ; row_init=0 ; inc_row=1 ; inc_col=1
      col= col_init ; row= row_init ; color= MPI_UNDEFINED
      do i=1, Inp_npes
         if (Ptopo_myproc == Ptopo_colrow(0,col,row)) then
            color=1
            ordinal_in_set = i
            exit
         endif
         col=col+1 ; row=row+1
         if (col>Ptopo_npex-1) then
            col= col_init ; row= row_init+inc_row
            inc_row= inc_row+1
         endif
         if (row>Ptopo_npey-1) then
            col= col_init+inc_col ; row= row_init
            inc_col= inc_col+1
         endif         
      enddo
      call MPI_COMM_SPLIT(COMM_grid,color,ordinal_in_set,Inp_COMM,err)
      me=-99
      if(Inp_COMM .ne. MPI_COMM_NULL) then 
         call mpi_comm_rank(Inp_COMM,me,err)
      endif
      Inp_iope= -1
      if (me == 0) Inp_bcast= Ptopo_myproc
      if (me >= 0) Inp_iope = me

      dim = (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*Inp_maxNKA
      WINSIZE   = dim * 4 ! will allocate dim reals of 4 bytes each
      DISP_UNIT = 4       ! window displacement unit = size of an integer
      INFO = 0            ! MPI_INFO_NULL

      call MPI_WIN_ALLOCATE (WINSIZE, DISP_UNIT, MPI_INFO_NULL,&
                             COMM_grid,basepntr, Inp_window, err)

      call C_F_POINTER ( basepntr, Inp_recv, [dim] )
      Inp_recv= 0.
!
!-------------------------------------------------------------------
!
      return
      end subroutine inp_init
      
end module inp_mod
