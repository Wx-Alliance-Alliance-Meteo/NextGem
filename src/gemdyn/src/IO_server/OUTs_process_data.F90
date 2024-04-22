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
      subroutine OUTs_process_data (F_deb, F_fin, F_nplans)
      use iso_c_binding
      use vGrid_Descriptors
      use MiMd
      use IOs
      use OUTs
      use omp_timing
      implicit none

      integer, intent(IN) :: F_deb, F_fin, F_nplans

      character(len=1024) fn
      character(len=3) stag,dumc3
      character(len=4) prefix,dumc4
      logical bool
      integer ORIGIN_COUNT, TARGET_COUNT, TARGET_RANK
      integer :: i,k,n,nko,nvar,skip,ii,jj,kk,offs,sk1,k01
      integer :: mpx,nklocal,irest,kstart,kend,k0,dim
      integer :: indo(10000),err,nfstecr,len,tag,cnt,ns
      integer :: status(MPI_STATUS_SIZE),sfc_lcl,sfc_i0,sfc_in
      integer(C_INTPTR_T) :: TARGET_DISP
      real :: gem_data(client_Hplane*G_nk*10)
      real, dimension (:,:), allocatable :: dataH
      type (sfc_var), dimension(:), allocatable :: vsfc
!     
!--------------------------------------------------------------------
!
      if (OUTs_1o1_L) print*, 'Proceeding with marker:', F_deb, F_fin
      
      call gtmg_start ( 20, 'GET_data', 13)
      if (Lun_out>0) call clock ( Lun_out, 'Recieving data', .false. )
      dim= client_Hplane*F_nplans
      allocate (dataH(dim,client_pestart:client_peend))
      allocate (vsfc(5000))
      do n=client_pestart,client_peend
         tag= 1001
         call MPI_recv ( dataH(1,n), dim, MPI_REAL, n, &
                         tag, MPI_COMM_WORLD, status, err)
      end do
      call gtmg_stop ( 20 )
      if (Lun_out>0) call clock (Lun_out, 'Recieving ... DONE', .false.)

      skip=0 ; k0=0 ; ns=0
      call gtmg_start ( 21, 'Processing', 13)

 987  if (skip>=F_fin) goto 888
      ns= ns+1
      prefix= TRANSFER(metaG(skip+1), dumc4)
      call up2low (prefix, dumc4) ; prefix=dumc4
      call OUTs_whichFST (Out_unf, prefix)
      
      Out_i0  = metaG(skip+2)
      Out_in  = metaG(skip+3)
      Out_j0  = metaG(skip+4)
      Out_jn  = metaG(skip+5)
      Out_reduc_L = TRANSFER(metaG(skip+6), bool)
      nvar= metaG(skip+7)
      skip= skip+7 ; cnt=0

      do i= 1, nvar
         call gtmg_start ( 31, 'Assembling', 21)
         Out_nomvar= TRANSFER(metaG(skip+1), dumc4)
         stag      = TRANSFER(metaG(skip+2), dumc3)
         call low2up (stag ,dumc3)
         stag  = dumc3
         Out_kind= metaG(skip+3)
         Out_nbit= metaG(skip+4)
         nko     = metaG(skip+5)
         indo(1:nko) = metaG(skip+6:skip+6+nko-1)
         nullify (levels)
         if (stag(2:2) == 'S') levels => Ver_i
         if (stag(2:2) == 'M') levels => Ver_hybM
         if (stag(2:2) == 'T') levels => Ver_hybT
         if (stag(2:2) == 'P') levels => Level_allpres
         if (stag(3:3) == 'D') then
            if (stag(2:2) == 'M') levels => hybM_diag
            if (stag(2:2) == 'T') levels => hybT_diag
            nko= 1 ; indo(1) = G_nk+1
         endif
         if (nko == 1) then
            cnt=cnt+1
            vsfc(cnt)%nv   = Out_nomvar
            vsfc(cnt)%knd  = Out_kind
            vsfc(cnt)%nbits= Out_nbit
            vsfc(cnt)%indx = i
            vsfc(cnt)%k0   = k0
            vsfc(cnt)%stag = stag
            vsfc(cnt)%skip = skip
            if ( stag(2:2) == 'S' ) then
               vsfc(cnt)%lvl= 0.
            else
               vsfc(cnt)%lvl= levels(G_nk+1)
            endif
            skip= skip+6+nko-1
            k0  = k0+nko
            cycle
         endif
         skip  = skip+6+nko-1
         call OUTs_wrtref (stag)
         
         do n=client_pestart,client_peend
            offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
            call lcl2glb (dataH(1,n), model_gindx(1,offs), k0, nko)
         end do
         
         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 31 )
         call gtmg_start ( 32, 'FSTECR', 21)

         call splitW ( myproc_IOS, numproc_IOS,nko,1,&
                       sfc_lcl, kstart, kend )
         call OUTs_fstecr (levels,indo,nko,kstart,kend,k0,stag(2:2))

         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 32 )
         k0= k0+nko
      enddo
      
      if (cnt>=1) then
         call gtmg_start ( 33, 'OneLevel', 21)
         call splitW ( myproc_IOS, numproc_IOS,cnt,1,&
                       sfc_lcl, sfc_i0, sfc_in)
         do i=1,cnt
         do n=client_pestart,client_peend
            offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
            call lcl2glb_sfc (dataH(1,n), model_gindx(1,offs), vsfc(i)%k0, i)
         end do
         end do
         call MPI_barrier (MY_WORLD_COMM,err)
         call OUTs_fstecr_sfc (vsfc,sfc_i0, sfc_in)
         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 33 )
      endif
      goto 987
      
 888  deallocate (dataH,vsfc)
      call gtmg_stop ( 21 )

      if (OUTs_1o1_L) print*, 'Marker:', F_deb, F_fin, ' ...DONE'
      if (Lun_out>0) call clock ( Lun_out, 'Processing ...DONE', .false. )
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_process_data

      subroutine lcl2glb (src, lcl_indx, F_k0, F_nk)
      use IOs
      use OUTs
      implicit none
      
      integer, intent(IN) :: lcl_indx(4), F_k0, F_nk
      real, intent(IN) :: src(*)

      integer i,j,k,cnt
!
!--------------------------------------------------------------------
!
      if ( F_nk > ubound(IOs_glbdata,3) ) then
         if (Lun_out>0) &
         print*, 'Insufficient storage in array IOs_glbdata --ABORT',&
                  F_nk, ubound(IOs_glbdata,3)
         stop
      endif
      
      do k= 1, F_nk
         cnt= F_k0*client_Hplane + (k-1)*client_Hplane
         do j= lcl_indx(3), lcl_indx(4)
            do i= lcl_indx(1), lcl_indx(2)
               cnt=cnt+1
               IOs_glbdata(i,j,k,IOS_couleur+1)= src(cnt)
            end do
         end do
      end do
!
!--------------------------------------------------------------------
!
      return
      end subroutine lcl2glb
      
      subroutine lcl2glb_sfc (src, lcl_indx, F_k0, F_k)
      use IOs
      use OUTs
      implicit none
      
      integer, intent(IN) :: lcl_indx(4), F_k0, F_k
      real, intent(IN) :: src(*)

      integer i,j,cnt
!
!--------------------------------------------------------------------
!
      if ( F_k > ubound(IOs_glbdata,3) ) then
         if (Lun_out>0) &
         print*, 'Insufficient storage in array IOs_glbdata --ABORT',&
                  F_k, ubound(IOs_glbdata,3)
         stop
      endif
      
      cnt= F_k0*client_Hplane
      do j= lcl_indx(3), lcl_indx(4)
         do i= lcl_indx(1), lcl_indx(2)
            cnt=cnt+1
            IOs_glbdata(i,j,F_k,IOS_couleur+1)= src(cnt)
         end do
      end do
!
!--------------------------------------------------------------------
!
      return
      end subroutine lcl2glb_sfc
