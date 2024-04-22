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
      module svro_mod
      use iso_c_binding
      use omp_timing
      use out_mod
      use, intrinsic :: iso_fortran_env
      implicit none
      public
      save

      logical :: OUTs_server_L=.false., OUTs_comm_L=.false.
      logical :: OUTs_this_step_L=.false.
      integer :: OUTs_COMM
      integer :: OUTs_window,OUTs_myHost,OUTs_blob
      integer :: me_1on1,OUTs_1on1, OUTs_cntwm, OUTs_nvar_indx
      integer :: OUTs_rank,OUTs_pe,OUTs_wm(4)
      integer :: OUTs_myHost_size,OUTs_myHost_rank,OUTs_request(3)
      integer, parameter :: OUTs_nplans =5000 ! maximum number of 2D plans
      integer, parameter :: OUTs_sorties=200  ! maximum number of sorties line in outcfg.out
      integer, parameter :: OUTs_nvar   =1000 ! maximum number of output variables
      integer :: max_meta=0, max_nplans=0, OUTs_err=0
      integer, dimension (:  ), pointer :: OUTs_meta
      integer, dimension (:,:), pointer :: OUTs_services => null()
      real   , dimension (:,:), pointer :: Shared_outbuf
      real   , dimension (:,:), pointer :: OUTs_data
      
contains
      subroutine itf_Oserv_init ()
      use glb_ld
      implicit none

      include 'mpif.h'
!
!------------------------------------------------------------------
!
      if (OUTs_server_L) then
         allocate (OUTs_meta(OUTs_sorties*7+OUTs_nvar*(6+G_nk*2)))
         allocate (OUTs_data(Out_Hmaxdim,OUTs_nplans))
         OUTs_request= MPI_REQUEST_NULL
      endif
!     
!------------------------------------------------------------------
!
      return         
      end subroutine itf_Oserv_init
      
!** s/r OUTs_start  - liste of output files to open

      subroutine OUTs_start ( )
      use dimout
      use ctrl
      use grdc_options
      use step_options
      use lun
      use out3
      use outd
      use outp
      use out_listes
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      include 'mpif.h'
      character(len=32) :: filenames_S(100)
      logical keep,ontimec,ontime_phy
      INTEGER(KIND=MPI_ADDRESS_KIND) TARGET_DISP
      integer ORIGIN_COUNT, TARGET_COUNT, TARGET_RANK,flag
      integer, dimension(0:MAXSET,Lctl_step:Lctl_step) :: casc_sorties
      integer, dimension(MAXSET) :: casc_levels
      integer :: i,nfiles,send(6),tag,n,wm,err,ireq
      integer :: status(mpi_status_size,Ptopo_numproc)
!     
!------------------------------------------------------------------
!
      if ( .not. OUTs_server_L) return

      ontimec= .false. ; ontime_phy= .false.
      if ((Grdc_ni>0).and.(Grdc_nj>0)) then
         if ( Lctl_step >= Grdc_start.and.Lctl_step <= Grdc_end ) then
            ontimec = (mod(Lctl_step+Grdc_start,Grdc_ndt) == 0)
         end if
      endif
      if ( Ctrl_phyms_L ) ontime_phy= outp_sorties(0,Lctl_step) > 0

      if ( (outd_sorties(0,Lctl_step) < 1) .and.&
           (.not. ontime_phy) .and. (.not. ontimec) ) return

      casc_sorties(0,Lctl_step)= 0
      if ( ontimec ) casc_sorties(0,Lctl_step)= 1
      
      call gtmg_start ( 77, 'OUTs_SYNC', 1)
      
      Out_nplans=0; Out_cntM= 0
      OUTs_cntwm = -1 ; OUTs_err= 0

      OUTs_this_step_L= .true.
      nfiles = 0
      
      call out_filenames (outd_sorties,Outd_lev,'d',&
                          lbound(outd_sorties,2)   ,&
                          ubound(outd_sorties,2),nfiles,filenames_S )
      if ( Ctrl_phyms_L ) &
         call out_filenames (outp_sorties,Outp_lev,'p',&
                             lbound(outp_sorties,2)   ,&
                             ubound(outp_sorties,2),nfiles,filenames_S )
      call out_filenames (casc_sorties,casc_levels,'casc',&
                          lbound(casc_sorties,2)   ,&
                          ubound(casc_sorties,2),nfiles,filenames_S )

      ORIGIN_COUNT = 1
      TARGET_COUNT = ORIGIN_COUNT
      TARGET_RANK  = OUTs_1on1
      TARGET_DISP  = 0

      call gemtime ( Lun_out, 'begin sync', .false. )
      ireq= 3
      call MPI_waitall (ireq,OUTs_request,status,err)

      if (OUTs_comm_L) then
         call MPI_barrier (OUTs_COMM, err)
         send(1) = nfiles
         send(2) = Out_ip2
         send(3) = Out_ip3
         send(4) = Out_npas
         send(5) = TRANSFER(Out_typvar_S(1:1), n)
         send(6) = Out_dateo
         n=nfiles*len(Out_filenames_S(1))
         tag= 601
         call MPI_send (send,size(send),MPI_INTEGER,&
                        OUTs_1on1,tag,OUTs_COMM,err)
         call MPI_send (Out_dirname_S,len(Out_dirname_S),&
                        MPI_CHARACTER,OUTs_1on1,tag+1,OUTs_COMM,err)
         call MPI_send (filenames_S,n,MPI_CHARACTER,&
                        OUTs_1on1,tag+2,OUTs_COMM,err)
      endif

      call gemtime ( Lun_out, 'Server-O open window', .false. )

      call gtmg_stop ( 77 )
!
!------------------------------------------------------------------
!
      return
      end subroutine OUTs_start

      subroutine OUTs_metaS ()
      use lun
      implicit none
      
      include 'mpif.h'
      integer :: n,err
!
!-------------------------------------------------------------------
!
      if (OUTs_comm_L) then
         OUTs_meta(Out_cntM+1) = TRANSFER(Out_prefix_S(1:4), n)
         OUTs_meta(Out_cntM+2) = Out_gridi0
         OUTs_meta(Out_cntM+3) = Out_gridin
         OUTs_meta(Out_cntM+4) = Out_gridj0
         OUTs_meta(Out_cntM+5) = Out_gridjn
         OUTs_meta(Out_cntM+6) = TRANSFER(Out_reduc_L, n)
         OUTs_meta(Out_cntM+7) = -1 ; OUTs_nvar_indx = Out_cntM+7
         Out_cntM= Out_cntM+7
      endif
!     
!-------------------------------------------------------------------
!
      return
      end subroutine OUTs_metaS
      
      subroutine OUTs_metaF (F_nfstecr, F_indx)
      implicit none
      integer, intent(IN) :: F_nfstecr, F_indx
!
!-------------------------------------------------------------------
!
      if (OUTs_comm_L) then
         if (F_nfstecr>0) then
            OUTs_meta(F_indx) = F_nfstecr
         else
            Out_cntM= Out_cntM-7
         endif
      endif
!     
!-------------------------------------------------------------------
!
      return
      end subroutine OUTs_metaF

      subroutine OUTs_end (F_sigready_L)
      use lun
      use ptopo
      implicit none

      logical, intent(IN) :: F_sigready_L
      include 'mpif.h'
      
      integer :: v1,v2,v3,err,tag,len,me,toto(100)
      integer :: status(mpi_status_size,Ptopo_numproc)
      integer, save :: cnt=0
!
!-------------------------------------------------------------------
!
      if ( .not. OUTs_server_L ) return
      
      if ( OUTs_this_step_L ) then

         call gemtime ( Lun_out, 'Server-O: SENDING1', .false. )
         call gtmg_start ( 79, 'OUTs_SENData', 1)
         if (OUTs_comm_L) then
         call MPI_barrier (OUTs_COMM, err)
         v1= 1 ; v2= Out_cntM ; v3= 0
         if (F_sigready_L) v3 = -1
         OUTs_wm(1:4) = (/v1,v2,v3,Out_nplans/)
         tag=901
         call MPI_iSend (OUTs_wm,size(OUTs_wm),MPI_INTEGER,OUTs_1on1,&
                         tag,OUTs_COMM,OUTs_request(1),err)
         if (Out_cntM>0) then
            call MPI_iSend (OUTs_meta,Out_cntM,MPI_INTEGER,OUTs_1on1,&
                           tag+1,OUTs_COMM,OUTs_request(2),err)
         endif
         endif
      
         call gemtime ( Lun_out, 'Server-O: SENDING2', .false. )
         tag= 1001 ; len= Out_Hmaxdim*Out_nplans
         if (Out_nplans>0) then
            call MPI_iSend (OUTs_data,len,MPI_REAL,OUTs_pe,&
                         tag,MPI_COMM_WORLD,OUTs_request(3),err )
         endif

         OUTs_this_step_L= .false.
         call gemtime ( Lun_out, 'Server-O: closed window', .false. )
         call gtmg_stop ( 79 )
         
      endif
!     
!-------------------------------------------------------------------
!
      return
      end subroutine OUTs_end

      subroutine itf_svro_fstecr ( fa,lminx,lmaxx,lminy,lmaxy,rf,nomvar,&
                                 mul,add,kind,lstep,nkfa,ind_o,nk_o,nbit)
      use lun
      use out_mod
      use out_meta
      use ptopo
      use glb_ld
      implicit none

      character(len=*), intent(IN) :: nomvar
      integer, intent(IN) :: lminx,lmaxx,lminy,lmaxy,&
                             nkfa,nbit,nk_o,kind,lstep
      integer, intent(IN) :: ind_o(nk_o)
      real   , intent(IN) :: fa(lminx:lmaxx,lminy:lmaxy,nkfa), &
                             rf(nkfa), mul,add
      include 'mpif.h'
      integer i,err
!
!----------------------------------------------------------------------
!
      if (OUTs_comm_L) then
         OUTs_meta(Out_cntM+1) = TRANSFER(nomvar(1:4)      , i)
         OUTs_meta(Out_cntM+2) = TRANSFER(Out_stag_S(1:3)  , i)
         OUTs_meta(Out_cntM+3) = kind
         OUTs_meta(Out_cntM+4) = nbit
         OUTs_meta(Out_cntM+5) = nk_o
         if (Out_stag_S(2:2)=='P') then
            call find_indx (OUTs_meta(Out_cntM+6:),rf,nk_o)
         else
            OUTs_meta(Out_cntM+6:Out_cntM+6+nk_o-1)= ind_o(1:nk_o)
         endif
         Out_cntM= Out_cntM+6+nk_o-1
         max_meta = max(max_meta,Out_cntM)
      endif

      max_nplans = max(max_nplans,Out_nplans+nk_o)
      if (Out_nplans+nk_o>OUTs_nplans) then
         if (Lun_out > 0) then
            OUTs_err=-1
            write(Lun_out,'("IOS_data buffer not large enough: ",2i6)' ) Out_nplans+1+nk_o,OUTs_nplans
         end if
         return
      endif
      
      call fa2rma (OUTs_data(1:,Out_nplans+1), fa, mul, add, ind_o, nk_o, &
                   lminx,lmaxx,lminy,lmaxy,nkfa)

      Out_nplans = Out_nplans + nk_o
      Out_nfstecr= Out_nfstecr+ 1
!
!--------------------------------------------------------------------
!
      return
      end subroutine itf_svro_fstecr

      subroutine fa2rma (dest, src, mul, add, ind_o, nk_o, &
                         lminx,lmaxx,lminy,lmaxy,nk)
      use glb_ld
      use out_mod
      use ptopo
      implicit none

      integer, intent(IN) :: lminx,lmaxx,lminy,lmaxy,nk,nk_o
      integer, intent(IN) :: ind_o(nk_o)
      real, intent(IN)  :: src(lminx:lmaxx,lminy:lmaxy,nk), mul,add
      real, intent(OUT) :: dest(*)

      integer i,j,k,cnt
    !  real w(lminx:lmaxx,lminy:lmaxy,nk),f2rc(G_ni,G_nj,nk)
!     
!--------------------------------------------------------------------
!
      do k= 1, nk_o
         cnt= (k-1)*Out_Hmaxdim
         do j= 1, l_nj
            do i= 1, l_ni
               cnt=cnt+1
               dest(cnt)= src(i,j,ind_o(k))*mul + add
            !   w(i,j,k) = dest(cnt)
            end do
         end do
      end do
      
!!$      call glbcolc2 ( f2rc, 1,G_ni,1,G_nj,1,nk,&
!!$                       w,lminx,lmaxx,lminy,lmaxy,1,nk)
!
!--------------------------------------------------------------------
!
      return
      end subroutine fa2rma

      subroutine find_indx (dest, pres, n)
      use out_mod
      use levels
      implicit none
      
      integer, intent(IN)  :: n
      real   , intent(IN)  :: pres(n)
      integer, intent(OUT) :: dest(n)

      integer i,k
!
!--------------------------------------------------------------------
!
      do k=1,n
         dest(k) = -99
         do i= 1, Level_npres
            if (Level_allpres(i)==pres(k)) then
               dest(k) = i
               exit
            endif
         end do
      end do
!
!--------------------------------------------------------------------
!
      return
      end subroutine find_indx
      subroutine create_on_node_shared_space ()
      use ptopo
      implicit none

      include 'mpif.h'
      integer :: i, cnt, tag
      integer :: ierr, gnumproc,gmyproc,me,DISP_UNIT, key, info
      integer, dimension (:,:), allocatable :: host_pe0s,G_host_pe0s
      integer(kind=INT64) :: hugedim,maxp
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      type(C_PTR), save :: basepntr
!     
!-------------------------------------------------------------------
!      
      call MPI_COMM_rank (COMM_WORLD,key ,ierr)
      CALL MPI_COMM_SPLIT_TYPE(COMM_world, MPI_COMM_TYPE_SHARED,&
                               key, info, OUTs_myHost, ierr)
      call MPI_COMM_size (OUTs_myHost, OUTs_myHost_size ,ierr)
      call MPI_COMM_rank (OUTs_myHost, OUTs_myHost_rank ,ierr)

      call MPI_COMM_size (COMM_world,gnumproc,ierr)
      call MPI_COMM_rank (MPI_COMM_WORLD,gmyproc ,ierr)
      call MPI_COMM_rank (COMM_WORLD,me ,ierr)
      
      allocate (host_pe0s(gnumproc,2),G_host_pe0s(gnumproc,2)) ; host_pe0s=0
      hugedim=0 ; maxp=OUTs_nplans
      if (OUTs_myHost_rank==0) then
         host_pe0s(me+1,1) = gmyproc+1
         host_pe0s(me+1,2) = OUTs_myHost_size
         hugedim= Out_Hmaxdim*maxp*OUTs_myHost_size
      endif

      call MPI_REDUCE (host_pe0s, G_host_pe0s, gnumproc*2, &
                       MPI_INTEGER,MPI_SUM,0,COMM_WORLD,ierr)
      if (OUTs_comm_L) then
         G_host_pe0s(:,1) = G_host_pe0s(:,1)-1
         host_pe0s= -1 ; cnt=0
         do i=1, gnumproc
            if (G_host_pe0s(i,1) >= 0) then
               cnt=cnt+1
               host_pe0s(cnt,1) = G_host_pe0s(i,1)
               host_pe0s(cnt,2) = G_host_pe0s(i,2)
            endif
         end do
         tag= 10
         call MPI_send ( host_pe0s,size(host_pe0s),MPI_INTEGER,&
                         OUTs_1on1,tag  ,OUTs_COMM,ierr )
      endif

      disp_unit = 4
      WINSIZE = hugedim * disp_unit
      call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                   OUTs_myHost, basepntr, OUTs_blob, ierr)
      call MPI_Win_shared_query(OUTs_blob, MPI_PROC_NULL, WINSIZE  ,&
                                disp_unit, basepntr,ierr)
      call C_F_POINTER ( basepntr, Shared_outbuf, [Out_Hmaxdim*OUTs_nplans,OUTs_myHost_size] )
!     
!-------------------------------------------------------------------
!
      return
      end subroutine create_on_node_shared_space
      
end module svro_mod
