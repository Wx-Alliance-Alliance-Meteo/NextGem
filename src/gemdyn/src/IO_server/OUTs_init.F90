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

!**s/r OUTs_init

      subroutine OUTs_init (grid_info)
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use IOs
      use OUTs
      implicit none

      integer, intent(IN) :: grid_info(*)
      
      character(len=1) :: dumc
      logical :: signal
      integer :: dim,myhost,ierr
      integer :: status(MPI_STATUS_SIZE),DISP_UNIT
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      type(C_PTR), save :: basepntr
!     
!--------------------------------------------------------------------
!
! Allocate shared memory among server PEs of the same node
      call MPI_Comm_split_type(MY_WORLD_COMM, MPI_COMM_TYPE_SHARED,0,&
                               MPI_INFO_NULL, myhost, ierr)
      call MPI_COMM_rank (myhost,OUTs_hostmyproc ,ierr)

      dim = 0
      if (OUTs_hostmyproc == 0) dim = G_ni*G_nj*(G_nk+1)*2*IOS_ncolors
      disp_unit = 4
      WINSIZE = dim * disp_unit
      call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                   myhost, basepntr, SHARED_WIN, ierr)
      call MPI_Win_shared_query(SHARED_WIN, MPI_PROC_NULL, WINSIZE  ,&
                                disp_unit, basepntr,ierr)
      call C_F_POINTER ( basepntr, IOs_glbdata, [G_ni,G_nj,(G_nk+1)*2,IOS_ncolors] )
      call MPI_barrier (MY_WORLD_COMM,ierr)
      
      Out_rewrit_L= TRANSFER(grid_info(23), signal)
      Out_deet = grid_info(24)
      Out_etik_S(1:12) = TRANSFER(grid_info(25), dumc)//&
                         TRANSFER(grid_info(26), dumc)//&
                         TRANSFER(grid_info(27), dumc)         
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_init

