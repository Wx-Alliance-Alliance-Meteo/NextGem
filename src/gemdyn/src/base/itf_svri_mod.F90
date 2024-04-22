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
      module svri_mod
      use iso_c_binding
      use omp_timing
      use, intrinsic :: iso_fortran_env
      use vGrid_Descriptors
      implicit none
      public
      save

      character(len=64), dimension(:), allocatable :: INs_list_S
      character(len=64), dimension(1000) :: INs_send_list_S
      logical :: INs_server_L=.false., INs_comm_L=.false.
      integer :: INs_COMM,INs_1on1,INs_rank,INs_pe,INs_nrequests
      integer, dimension (:,:), pointer :: INs_services => null()
      integer :: GZcaract(5),INs_gmyproc,INs_dimgzH, INs_dimgz3d
      integer :: INs_request(2),INs_sent_nreqs,INs_recv_nreqs,&
                 INs_nplans,INs_irecv(5),INs_maxNKA, INs_maxreqs
      integer, dimension(:), pointer :: INs_GZIP1,INs_DIP1
      real, dimension(:), pointer :: INs_GZ,INs_ND

      character(len=32), dimension(:), allocatable :: cBUF
      integer, dimension (:), pointer :: iBUF
      real   , dimension (:), allocatable :: GZbuf, NDbuf
      real(kind=REAL64), pointer :: VGD_tbl_8(:)
      
      real   , dimension (:), pointer :: Rcv_GZr,Rcv_NDr
      type :: REQ
         character(len=32) :: vname(2)
         integer :: nk,deb
      end type REQ
      type(REQ), dimension(:), allocatable :: SRL
      
contains
      subroutine itf_Iserv_init ()
      use gem_options
      use glb_ld
      use inp_options
      use path
      use ptopo
      use tr3d
      implicit none

      include 'mpif.h'
      character(len=7) vname
      integer:: i,n1,n2,dim,cnt,tag,err
!
!------------------------------------------------------------------
!
      call MPI_COMM_rank (MPI_COMM_WORLD,INs_gmyproc ,err)
      INs_dimgzH = (l_ni+2*G_halox)*(l_nj+2*G_haloy)
      allocate (INs_list_S(1000))
      INs_irecv= MPI_REQUEST_NULL
      
      if (INs_server_L) then
         INs_list_S(1)= ''
         INs_list_S(2)= 'TEMPERATURE:Q:R'
         INs_list_S(3)= 'SFCPRES:Q:R'
         INs_list_S(4)= 'QT1:Q:R'
         INs_list_S(5)= 'WT1:Q:R'
         INs_list_S(6)= 'ZDT1:Q:R'
         INs_list_S(7)= 'UV:UV:R'
         cnt=7
         do i=1,Tr3d_ntr
            vname= 'TR/'//trim(Tr3d_name_S(i))
            if ( any (Inp_blacklist_S(1:MAX_blacklist) == trim(vname)) ) cycle
            cnt= cnt+1
            INs_list_S(cnt)= trim(vname)//':Q:R'
         end do
         INs_nrequests= cnt
         INs_request= MPI_REQUEST_NULL
         if (INs_comm_L) then
            tag = 121
            call MPI_send (Path_input_S,len(Path_input_S),MPI_CHARACTER,&
                           INs_1on1,tag  ,INs_COMM,err )
            call MPI_send (INs_nrequests,1,MPI_INTEGER,&
                           INs_1on1,tag+1 ,INs_COMM,err )
            call MPI_recv (INs_maxNKA, 1, MPI_INTEGER,INs_1on1,tag+2,&
                           INs_COMM, MPI_STATUSES_IGNORE, err)
         endif
      endif
! In the dynamics the number of requests is constant (INs_nrequests)
! and therefore the I-svr is set with that value. Be aware of that
! fact if this server is used with physics input.

      call MPI_bcast (INs_maxNKA,1,MPI_INTEGER, 0,&
                      COMM_multigrid, err)               
      INs_maxreqs= INs_nrequests
! Memory allocation here is set to a constant maximum value
! which must exactly match the I-svr allocation. It cannot be
! dynamic because of the non-blocking nature of the isend/irecv
      dim= 2*INs_maxreqs + INs_maxNKA*INs_maxreqs
      allocate (iBUF(dim))
      allocate (VGD_tbl_8(10*INs_maxNKA))
      n1= Ptopo_gindx(2,1)-Ptopo_gindx(1,1)+1+2*G_haloy
      n2= Ptopo_gindx(4,1)-Ptopo_gindx(3,1)+1+2*G_haloy
      dim= n1*n2*(INs_maxNKA*6 + 1)
      allocate (Rcv_GZr(dim))
      dim= n1*n2*INs_maxNKA*INs_maxreqs
      allocate (Rcv_NDr(dim))
      allocate (cBUF(INs_maxreqs*2),SRL(INs_maxreqs))
!     
!------------------------------------------------------------------
!
      return         
      end subroutine itf_Iserv_init
      
!**s/r itf_Iserv_request - Send a request for dateV

      subroutine itf_Iserv_request ( F_dateV_S, F_list_S, F_tag, F_nr )
      implicit none

      character(len=16 ),intent(IN) :: F_dateV_S
      integer, intent(IN) :: F_tag,F_nr
      character(len=*), intent(INOUT) :: F_list_S(F_nr)
      
      include 'mpif.h'
      character(len=4) :: tag_S
      integer :: tag,err
      integer :: status(MPI_STATUS_SIZE)      
!
!------------------------------------------------------------------
!
      if (INs_comm_L) then
         call MPI_waitall (size(INs_request),INs_request,status,err)
         print*, 'SENDING a request to IN-server for datev: ', F_dateV_S
         tag= 801
         write(tag_S,'(i4)') F_tag
         F_list_S(1)= "VALID_DATE:"//trim(F_dateV_S)//":"//tag_S
         INs_sent_nreqs = F_nr
         INs_send_list_S(1:F_nr)= F_list_S(1:F_nr)
         call MPI_isend (INs_sent_nreqs,1,MPI_INTEGER,&
                         INs_1on1,tag,INs_COMM,INs_request(1),err )
         call MPI_isend (INs_send_list_S,F_nr*len(INs_send_list_S(1)),&
               MPI_CHARACTER,INs_1on1,tag+1,INs_COMM,INs_request(2),err)
      endif
! Immediately post the irecv            
      call itf_Iserv_recv ()
!     
!------------------------------------------------------------------
!
      return         
      end subroutine itf_Iserv_request
      
      subroutine itf_Iserv_recv ()
      use inp_mod
      use ptopo
      implicit none

      include 'mpif.h'
      integer:: tag,err
!     
!------------------------------------------------------------------
!
      if (INs_server_L) then
         if (INs_comm_L) then
            tag= 3001
            call MPI_irecv ( cBUF,size(cBUF)*len(cBUF(1)),MPI_CHARACTER,&
                          INs_1on1,tag,INs_COMM,INs_irecv(1),err)
            call MPI_irecv ( iBUF,size(iBUF),MPI_INTEGER, INs_1on1,&
                          tag+1,INs_COMM,INs_irecv(2),err)
            call MPI_irecv (VGD_tbl_8,size(VGD_tbl_8),MPI_DOUBLE_PRECISION,&
                            INs_1on1,tag+2,INs_COMM,INs_irecv(3),err)
         endif
         tag= 4001
         !print*, 'RCV_GZ: ',size(Rcv_GZr),tag+INs_gmyproc,INS_pe
         call MPI_irecv ( Rcv_GZr,size(Rcv_GZr),MPI_REAL, INS_pe,&
                          tag+INs_gmyproc,MPI_COMM_WORLD,INs_irecv(4),err)
         tag= 9001
         !print*, 'RCV_ND: ',size(Rcv_NDr),tag+INs_gmyproc,INS_pe
         call MPI_irecv ( Rcv_NDr,size(Rcv_NDr),MPI_REAL, INS_pe,&
                          tag+INs_gmyproc,MPI_COMM_WORLD,INs_irecv(5),err)
      endif
!     
!------------------------------------------------------------------
!
      return
      end subroutine itf_Iserv_recv

      subroutine itf_Iserv_GZ3d ()
      use gem_options
      use glb_ld
      use inp_mod
      use ptopo
      implicit none

      include 'mpif.h'
      integer :: n1,n2,n3,nm,nreqs,nplans,err,n
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
!     
!------------------------------------------------------------------
!
      Inp_src_GZ_L= .false.
      if (INs_server_L) then
         call MPI_bcast (cBUF,size(cBUF)*len(cBUF(1)),MPI_CHARACTER, 0,&
                         COMM_multigrid, err)         
         call MPI_bcast (iBUF,size(iBUF),MPI_INTEGER, 0,&
                         COMM_multigrid, err)         
         call MPI_bcast (VGD_tbl_8,size(VGD_tbl_8),MPI_DOUBLE_PRECISION, 0,&
                         COMM_multigrid, err)         
         GZ3d%nk     = iBUF(1)
         GZ3d%kind   = iBUF(2)
         GZ3d%me_L   = iBUF(3)==1
         GZ3d%mels_L = iBUF(4)==1
         Inp_pref_a_8= iBUF(5)
         GZ3d%ip1 => iBUF(6:5+GZ3d%nk) ; nm=5+GZ3d%nk
         INs_recv_nreqs= iBUF(nm+1)
         INs_nplans    = iBUF(nm+2) ; nm=nm+2
         INs_dimgz3d   = INs_dimgzH*GZ3d%nk
         nreqs= INs_recv_nreqs ; nplans= INs_nplans
         SRL(1:nreqs)%nk = iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         SRL(1:nreqs)%deb= iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         INs_DIP1 => iBUF(nm+1:INs_nplans)     ; nm=nm+INs_nplans
         n1= iBUF(nm+1) ; n2= iBUF(nm+2) ; n3= iBUF(nm+3)
         allocate ( vtbl_8(n1,n2,n3) )
         vtbl_8 = reshape(VGD_tbl_8,(/n1,n2,n3/))
         err= vgd_new ( Inp_vgd_src, vtbl_8 )
         deallocate (vtbl_8)

         SRL(1:nreqs)%vname(1)= cBUF(1:nreqs)                 
         SRL(1:nreqs)%vname(2)= cBUF(nreqs+1:2*nreqs)
         
         INS_GZ(1:) => Rcv_GZr(1:)
         INS_ND(1:) => Rcv_NDr(1:)
         
         if (associated(GZ3d%valq)) deallocate (GZ3d%valq)
         if (associated(GZ3d%valu)) deallocate (GZ3d%valu)
         if (associated(GZ3d%valv)) deallocate (GZ3d%valv)
         if (associated(GZ3d%valv)) deallocate (GZ3d%sfc )
         nullify (GZ3d%valq,GZ3d%valu,GZ3d%valv,GZ3d%sfc)

         allocate (GZ3d%valq(l_minx:l_maxx,l_miny:l_maxy,GZ3d%nk),&
                   GZ3d%valu(l_minx:l_maxx,l_miny:l_maxy,GZ3d%nk),&
                   GZ3d%valv(l_minx:l_maxx,l_miny:l_maxy,GZ3d%nk),&
                   GZ3d%sfc (l_minx:l_maxx,l_miny:l_maxy,2))
         GZ3d%valq= 0. ; GZ3d%sfc= 0. ! for -C to work properly
         call reshapeH ( INS_GZ                  ,GZ3d%valq  ,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,GZ3d%nk     )
         call reshapeH ( INS_GZ(  INs_dimgz3d+1:),GZ3d%valu  ,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,GZ3d%nk     )
         call reshapeH ( INS_GZ(2*INs_dimgz3d+1:),GZ3d%valv  ,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,GZ3d%nk     )
         GZ3d%sfc(:,:,1) = GZ3d%valq(:,:,GZ3d%nk)
         call reshapeH ( INS_GZ(3*INs_dimgz3d+1:),&
                         GZ3d%sfc(l_minx:,l_miny:,2),&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,1           )
         Inp_src_GZ_L= .true.
      endif
!     
!------------------------------------------------------------------
!
      return
      end subroutine itf_Iserv_GZ3d
            
      subroutine reshapeH (src,dst,s1,s2,s3,s4,d1,d2,d3,d4,nk)
      implicit none
      
      integer, intent(IN) :: s1,s2,s3,s4,d1,d2,d3,d4,nk
      real, intent(IN ) :: src(s1:s2,s3:s4,nk)
      real, intent(OUT) :: dst(d1:d2,d3:d4,nk)
      
      integer i,j
!     
!------------------------------------------------------------------
!
      do j=s3,s4
         do i=s1,s2
            dst(i,j,:) = src(i,j,:)
         end do
      end do
!     
!------------------------------------------------------------------
!
      return
      end subroutine reshapeH
      
end module svri_mod
