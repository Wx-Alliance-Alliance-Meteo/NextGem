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
      module IOserver
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      implicit none
      public
      save
      
contains

      subroutine SVR_init (F_server_S, F_COMM, F_server_L, F_comm_L, &
                           F_1on1, F_rank, F_pe, F_services)
      use gem_options
      use domains
      use glb_ld
      use glb_pil
      use MiMd
      use geomh
      use hgc
      use lun
      use out_mod
      use ptopo
      use levels
      use tr3d
      use ver
      use vGrid_Descriptors
      implicit none

      character(len=*), intent(IN) :: F_server_S
      integer         , intent(IN) :: F_COMM
      logical, intent(OUT) :: F_server_L,F_comm_L
      integer, intent(OUT) :: F_1on1, F_rank, F_pe
      integer, dimension (:,:), pointer, intent(INOUT) :: F_services
      
      include 'mpif.h'

      logical :: Server_OK_L
      integer :: grid_info(100), err
      integer :: i, tag, myproc, me_1on1
      integer :: status(mpi_status_size,1)
      real(kind=REAL64), pointer :: vtbl_8(:,:,:) => null()
!     
!-------------------------------------------------------------------
!
      F_server_L= .false.
      do i=1,MiMd_ncolors
         if ( trim(MiMd_world(i)%name_S)==trim(F_server_S) ) then
            F_server_L=.true.
            F_rank=MiMd_world(i)%rank ! rank of Server
            exit
         endif
      end do

      if (F_server_L) then
         
         F_comm_L = F_COMM .ne. MPI_COMM_NULL

         allocate (F_services(3,F_rank))
     
         err= vgd_get ( Ver_vgdobj, 'VTBL', vtbl_8, quiet=.true.)
         grid_info(1) = G_ni
         grid_info(2) = G_nj
         grid_info(3) = G_nk
         grid_info(4) = G_halox
         grid_info(5) = G_haloy
         grid_info(6) = l_minx 
         grid_info(7) = l_maxx 
         grid_info(8) = l_miny
         grid_info(9) = l_maxy 
         grid_info(10)= Tr3d_ntr
         grid_info(11)= TRANSFER(Tr3d_anydate_L,i)
         grid_info(12)= Glb_pil_w
         grid_info(13)= Glb_pil_e
         grid_info(14)= Glb_pil_s
         grid_info(15)= Glb_pil_n
         grid_info(16)= Out_Hmaxdim
         grid_info(17)= Domains_ngrids
         grid_info(18:20) = ubound(vtbl_8)
         grid_info(21)= Level_npres
         grid_info(22)= size(geomh_4output)
         grid_info(23)= TRANSFER(Out_rewrit_L, i)
         grid_info(24)= Out_deet
         grid_info(25)= TRANSFER(Out_etik_S(1:4 ), i)
         grid_info(26)= TRANSFER(Out_etik_S(5:8 ), i)
         grid_info(27)= TRANSFER(Out_etik_S(9:12), i)

         if (F_comm_L) then
            call MPI_COMM_rank (F_COMM,me_1on1,err)
            F_1on1 = -1*(me_1on1-1)
            call MPI_barrier (F_COMM, err)
            tag= 110
            call MPI_recv ( Server_OK_L, 1, MPI_LOGICAL, &
                            F_1on1, tag, F_COMM, status, err)
           
            if (.not. Server_OK_L) goto 999
            tag= 111
            call MPI_send ( grid_info,size(grid_info),MPI_INTEGER     ,&
                            F_1on1,tag  ,F_COMM,err )
            call MPI_send ( Ptopo_gindx,size(Ptopo_gindx),MPI_INTEGER ,&
                            F_1on1,tag+1,F_COMM,err )
            call MPI_send ( geomh_4output,size(geomh_4output),MPI_REAL,&
                            F_1on1,tag+2,F_COMM,err )
            call MPI_send ( Level_allpres,Level_npres,MPI_REAL        ,&
                            F_1on1,tag+3,F_COMM,err )
            call MPI_send ( vtbl_8,size(vtbl_8),MPI_DOUBLE_PRECISION  ,&
                            F_1on1,tag+4,F_COMM,err )
            tag= 311
            call MPI_recv ( F_services, size(F_services), MPI_INTEGER,&
                            F_1on1, tag, F_COMM, status, err)
            call MPI_barrier (F_COMM, err)
         endif

 999     call MPI_bcast ( Server_OK_L, 1, MPI_LOGICAL,0,&
                          COMM_multigrid, err )

         if (.not. Server_OK_L) then
            if (Lun_out>0) write(Lun_out,'(/10("#")/2x,2a,/,10("#"))')&
                           trim(F_server_S),' has disconnected'
            F_server_L = .false.
            F_comm_L = .false.
            return
         endif
                        
         call MPI_bcast ( F_services, size(F_services), MPI_INTEGER,0,&
                          COMM_multigrid, err )       

         call MPI_comm_rank (MPI_COMM_WORLD,myproc,err)
         do i=1,F_rank
            if ( (myproc>=F_services(2,i)) .and. &
                 (myproc<=F_services(3,i)) ) then
               F_pe = F_services(1,i)
               exit
            endif
         end do

      else
         if (Lun_out >= 0) write(6,'(/2x3a/)') 'SVR_init: SERVER ',trim(F_server_S),' unavailable'
      endif
!     
!-------------------------------------------------------------------
!
      return
      end subroutine SVR_init
      
end module IOserver
