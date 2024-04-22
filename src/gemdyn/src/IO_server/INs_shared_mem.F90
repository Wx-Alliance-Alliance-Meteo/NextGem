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

!**s/r INs_shared_mem - Determine analysis dimensuins and allocate
!                       shared memory among server PEs (same node)

      logical function INs_shared_mem (F_unf, F_vname_S)
      use IOs
      use INs
      implicit none

      character(len=*), intent(IN) :: F_vname_S
      integer         , intent(IN) :: F_unf

#include <rmnlib_basics.hf>

      integer, parameter :: nlis = 1024
      integer :: nrec,n1,n2,n3,n4,lislon,liste(nlis)
      integer :: dim,DISP_UNIT,ierr
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      type(C_PTR), save :: basepntr
!
!---------------------------------------------------------------------
!
      INs_shared_mem= .false.
      nrec= fstinl (F_unf, n1,n2,n3, Inp_cmcdate,' ', &
                    -1,-1,-1,' ', F_vname_S,liste,lislon,nlis)
      
      if (lislon>0) then
         INs_nia= n1
         INs_nja= n2
         INs_nka= lislon
         INs_shared_mem= .true.
         if (INs_1o1_L) print*, 'Analysis at ',Inp_datev,&
         ' dimensions: ',INs_nia,INs_nja,INs_nka,&
         ' from TT in unit: ', F_unf
         n1= INs_nid
         n2= INs_njd
         n3= (INs_nka+1)*6
         dim = 0
         if (INs_hostmyproc == 0) dim= INs_hord*n3
         disp_unit = 4
         WINSIZE = dim * disp_unit
         call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                      INs_host, basepntr, GZ_win, ierr)
         call MPI_Win_shared_query(GZ_win, MPI_PROC_NULL, WINSIZE  ,&
                                   disp_unit, basepntr, ierr)
         call C_F_POINTER (basepntr,GZ,[n1,n2,n3] )
         call MPI_barrier (MY_WORLD_COMM,ierr)
         allocate (GZIP1(2*(INs_nka+2)))

         n3= (INs_nreq+2)*(INs_nka+2)
         if (dble(INs_hord*n3) > dble(huge(n3))*.98) then
            print*,'OUT of memory in INs_shared_mem'
            INs_shared_mem= .false.
            return
         endif

         dim = 0
         if (INs_hostmyproc == 0) dim= INs_hord*n3
         disp_unit = 4
         WINSIZE = dim * disp_unit
         call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                      INs_host, basepntr, NEST_win, ierr)
         call MPI_Win_shared_query(NEST_win, MPI_PROC_NULL, WINSIZE  ,&
                                   disp_unit, basepntr, ierr)
         dim= INs_hord*n3
         call C_F_POINTER (basepntr,ND,[dim] )
         call MPI_barrier (MY_WORLD_COMM,ierr)
         allocate (DIP1(n3*2))
      else
         print*,'could NOT determine Input dimensions with var= ',trim(F_vname_S) 
      endif
!
!---------------------------------------------------------------------
!
      return
      end function INs_shared_mem
      
