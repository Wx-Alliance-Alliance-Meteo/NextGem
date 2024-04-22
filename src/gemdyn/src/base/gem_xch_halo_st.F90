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

!**s/r gem_xch_halo_st - MPI halo exchange based on 2-sided MPI comm

      subroutine gem_xch_halo_st ( f, lminx,lmaxx,lminy,lmaxy,nk, F_halo )
      use ptopo
      use glb_ld
      implicit none

      integer, intent (IN) :: lminx,lmaxx,lminy,lmaxy,nk,F_halo
      real, intent (INOUT) :: f(lminx:lmaxx,lminy:lmaxy,nk)

      include 'mpif.h'
      integer i, j, k, lni,lnj, len,lenb, halo, tid
      integer i0,in,j0,jn,proc,request(4),ierr
      integer tag1,tag2,ireq
      real, dimension((G_maxldim+2*(1-lminx))*(1-lminx)*nk,2) :: buf,recv
!
!---------------------------------------------------------------------
!
      request= MPI_REQUEST_NULL
      halo= F_halo ; if (halo<0) halo= 1-lminx
      if ((nk<1).or.(halo<1)) return
      ireq= 0
      tag1= 8528 ; tag2= tag1*2
! Exchanging first on the west-east axis
      j0= 1-halo*south ; jn= l_nj+halo*north ; lnj=jn-j0+1
      lenb= halo*lnj*nk
      if (.not. l_east) then
         ireq= ireq+1
         proc= Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow)
         call MPI_IRecv (recv(1,1), lenb, MPI_REAL, proc,tag1+proc,&
                                 COMM_grid,request(ireq),ierr)
         ireq= ireq+1 ; len= 0
         do k = 1, nk
            do j = j0, jn
               do i = l_ni-halo+1, l_ni
                  len      = len + 1
                  buf(len,2) = f(i,j,k)
               end do
            end do
         end do
         call MPI_ISend (buf(1,2),len,MPI_REAL,&
                         proc,tag1+Ptopo_myproc,&
                         COMM_grid, request(ireq),ierr )         
      endif
      if (.not. l_west) then
         ireq= ireq+1
         proc= Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow)
         call MPI_IRecv (recv(1,2), lenb, MPI_REAL, proc,tag1+proc,&
                                 COMM_grid,request(ireq),ierr)
         ireq = ireq+1 ; len= 0
         do k = 1, nk
            do j = j0, jn
               do i = 1, halo
                  len      = len + 1
                  buf(len,1) = f(i,j,k)
               end do
            end do
         end do
         call MPI_ISend (buf(1,1),len,MPI_REAL,&
                         proc,tag1+Ptopo_myproc,&
                         COMM_grid, request(ireq),ierr )
      endif

      call MPI_waitall (ireq,request,MPI_STATUSES_IGNORE,ierr)

      if (.not. l_west) then
         len = 0
         do k = 1, nk
            do j = j0, jn
               do i = 1-halo, 0
                  len      = len + 1
                  f(i,j,k) = recv(len,2)
               end do
            end do
         end do
      endif

      if (.not. l_east) then
         len = 0
         do k = 1, nk
            do j = j0, jn
               do i = l_ni+1, l_ni+halo
                  len      = len + 1
                  f(i,j,k) = recv(len,1)
               end do
            end do
         end do
      endif
      
! Then exchanging on the south-north axis and include halos for corners
      i0= 1-halo ; in= l_ni+halo ; lni=in-i0+1
      lenb= halo*lni*nk ; ireq= 0
      if (.not. l_north) then
         ireq= ireq+1
         proc= Ptopo_colrow(Ptopo_couleur,Ptopo_mycol,Ptopo_myrow+1)
         call MPI_IRecv (recv(1,1), lenb, MPI_REAL, proc,tag2+proc,&
                                 COMM_grid,request(ireq),ierr)
         ireq= ireq+1 ; len= 0
         do k = 1, nk
            do j = l_nj-halo+1, l_nj
               do i = i0, in
                  len      = len + 1
                  buf(len,2) = f(i,j,k)
               end do
            end do
         end do
         call MPI_ISend (buf(1,2),len,MPI_REAL,&
                         proc,tag2+Ptopo_myproc,&
                         COMM_grid, request(ireq),ierr )         
      endif
      if (.not. l_south) then
         ireq= ireq+1
         proc= Ptopo_colrow(Ptopo_couleur,Ptopo_mycol,Ptopo_myrow-1)
         call MPI_IRecv (recv(1,2), lenb, MPI_REAL, proc,tag2+proc,&
                                 COMM_grid,request(ireq),ierr)
         ireq = ireq+1 ; len= 0
         do k = 1, nk
            do j = 1, halo
               do i = i0, in
                  len      = len + 1
                  buf(len,1) = f(i,j,k)
               end do
            end do
         end do
         call MPI_ISend (buf(1,1),len,MPI_REAL,&
                         proc,tag2+Ptopo_myproc,&
                         COMM_grid, request(ireq),ierr )
      endif

      call MPI_waitall (ireq,request,MPI_STATUSES_IGNORE,ierr)

      if (.not. l_north) then
         len = 0
         do k = 1, nk
            do j = l_nj+1, l_nj+halo
               do i = i0, in
                  len      = len + 1
                  f(i,j,k) = recv(len,1)
               end do
            end do
         end do
      endif

      if (.not. l_south) then
         len = 0
         do k = 1, nk
            do j = 1-halo, 0
               do i = i0, in
                  len      = len + 1
                  f(i,j,k) = recv(len,2)
               end do
            end do
         end do
      endif
!     
!---------------------------------------------------------------------
!
      return
      end subroutine gem_xch_halo_st
