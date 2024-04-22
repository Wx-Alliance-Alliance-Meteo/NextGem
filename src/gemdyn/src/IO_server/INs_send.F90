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

!**s/r INs_send

      subroutine INs_send (F_GZ,F_ND,Minx,Maxx,Miny,Maxy,F_err)
      use, intrinsic :: iso_fortran_env
      use iso_c_binding
      use omp_timing
      use IOs
      use INs
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,F_err
      real, intent(IN) :: F_GZ(Minx:Maxx,Miny:Maxy,*),&
                          F_ND(Minx:Maxx,Miny:Maxy,*)
      
      integer :: i,j,k,cnt,n,n1,n2,offs,tag,clients,dim,nm,err
         integer i0,in,j0,jn                 
!     
!--------------------------------------------------------------------
!
      if (Lun_out>0) call clock ( Lun_out, 'Waitall', .false. )
      call gtmg_start ( 20, 'waitall', 13)
      call MPI_waitall (size(INs_isend),INs_isend,&
                        MPI_STATUSES_IGNORE,err)
      call gtmg_stop ( 20 )
      if (F_err<0) GZcaract=-9

      n1= model_gindx(2,1)-model_gindx(1,1)+1+2*G_haloy
      n2= model_gindx(4,1)-model_gindx(3,1)+1+2*G_haloy
      dim= n1*n2*(GZcaract(1)*3 + 1)
      if (dim>ubound(GZbuf,1)) then
         print*, 'INS_send: NOT enough space in GZbuf - abort'
         print*, 'INS_send: required: ',dim,' allocated: ',ubound(GZbuf,1)
         stop
      endif
      
      if (INs_1o1_L) then
         VGD_tbl_8(:) = reshape(vtbl_8,(/size(vtbl_8)/))
         do n=1,INs_nreq
            cBUF(         n) = SRL(n)%vname(1)
            cBuf(INs_nreq+n) = SRL(n)%vname(2)
         end do
         n=0 ; tag= 3001
         iBUF(n+1:n+size(GZcaract)) = GZcaract ; n= n+size(GZcaract)
         iBUF(n+1:n+GZcaract(1)   ) = GZIP1(1:GZcaract(1)) ; n= n+GZcaract(1)
         iBUF(n+1) = INs_nreq ; iBUF(n+2) = INs_nplans ; n= n+2 ; nm=n
         iBUF(nm+1:nm+INs_nreq)= SRL(1:INs_nreq)%nk    ; nm=nm+INs_nreq
         iBUF(nm+1:nm+INs_nreq)= SRL(1:INs_nreq)%deb   ; nm=nm+INs_nreq
         iBUF(nm+1:nm+INs_nplans) = DIP1(1:INs_nplans) ; nm=nm+INs_nplans
         iBUF(nm+1:nm+3) = INs_n123 ; nm=nm+3
         call MPI_isend ( cBUF,size(cBUF)*len(cBUF(1)),MPI_CHARACTER,&
                          INs_gem1o1,tag,INs_GEM_COMM,INs_isend(1),err)
         call MPI_isend ( iBUF,size(iBUF),MPI_INTEGER, INs_gem1o1,&
                          tag+1,INs_GEM_COMM,INs_isend(2),err)
         call MPI_isend (VGD_tbl_8,size(VGD_tbl_8),MPI_DOUBLE_PRECISION,&
                         INs_gem1o1,tag+2,INs_GEM_COMM,INs_isend(3),err)
      endif
      if (Lun_out>0) call clock ( Lun_out, 'SEND_completed1', .false. )

      clients= 3
      call MPI_barrier (MY_WORLD_COMM,err)
      do n=client_pestart,client_peend
         tag= 4001
         offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
         i0= model_gindx(1,offs)-G_halox
         in= model_gindx(2,offs)+G_halox
         j0= model_gindx(3,offs)-G_haloy
         jn= model_gindx(4,offs)+G_haloy
         dim= (in-i0+1)*(jn-j0+1)*(GZcaract(1)*3 + 1)
         cnt=0
         do k=1, GZcaract(1)*3 + 1
         do j=model_gindx(3,offs)-G_haloy,model_gindx(4,offs)+G_haloy
         do i=model_gindx(1,offs)-G_halox,model_gindx(2,offs)+G_halox
            cnt= cnt+1
            GZbuf(cnt,n)= F_GZ(i,j,k)
         end do
         end do
         end do
         clients= clients+1
         !print*, 'SND_GZ: ',cnt,tag+n,n,clients,ubound(GZbuf,1)
         cnt= ubound(GZbuf,1)
         call MPI_isend ( GZbuf(1,n),cnt,MPI_REAL,n,tag+n,&
                          MPI_COMM_WORLD,INs_isend(clients),err )
         tag= 9001
         cnt=0
         do k=1,INs_nplans
         do j=model_gindx(3,offs)-G_haloy,model_gindx(4,offs)+G_haloy
         do i=model_gindx(1,offs)-G_halox,model_gindx(2,offs)+G_halox
            cnt= cnt+1
            NDbuf(cnt,n)= F_ND(i,j,k)
         end do
         end do
         end do
         clients= clients+1
         !print*, 'SND_ND: ',cnt,tag+n,n,clients,ubound(NDbuf,1)
         cnt= ubound(NDbuf,1)
         call MPI_isend ( NDbuf(1,n),cnt,MPI_REAL,n,tag+n,&
                          MPI_COMM_WORLD,INs_isend(clients),err)
      end do
      if (Lun_out>0) call clock ( Lun_out, 'SEND_completed2', .false. )
!     
!--------------------------------------------------------------------
!
      return
      end subroutine INs_send

