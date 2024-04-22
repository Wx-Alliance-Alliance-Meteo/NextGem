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

!**s/r glb_orderd_sum 
!
      subroutine glb_orderd_sum (F_sum, F_src, Minx, Maxx, Miny, Maxy,&
                                 Nk,i0,in,j0,jn,k0,kn)
      use glb_ld
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(IN) :: Minx,Maxx,Miny,Maxy,Nk,i0,in,j0,jn,k0,kn
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(IN) :: F_src
      real(kind=REAL64), intent(OUT) :: F_sum

      include 'mpif.h'
      integer :: i,j,k,err,tag,n,nj,proc
      integer :: status(mpi_status_size,1)
      real(kind=REAL64) :: sum_lcl(l_nj,0:Ptopo_npey-1)
!
!----------------------------------------------------------------------
!
!!$omp single
      sum_lcl=0.d0 ; F_sum=0.d0
      if (Ptopo_mycol==0) then
         do k=k0,kn
         do j=j0,jn
         do i=i0,in
            sum_lcl(j,Ptopo_myrow) = sum_lcl(j,Ptopo_myrow) + F_src(i,j,k)
         end do
         end do
         end do
         if (me_in_row/=num_in_row-1) then
            tag=711+me_in_row
            call MPI_send ( sum_lcl(1,Ptopo_myrow), l_nj, MPI_DOUBLE_PRECISION, me_in_row+1, tag, COMM_row,err )
         endif
      else
         tag=711+me_in_row-1
         call MPI_recv ( sum_lcl(1,Ptopo_myrow), l_nj, MPI_DOUBLE_PRECISION,me_in_row-1,tag, COMM_row, status, err)
         do k=k0,kn
         do j=j0,jn
         do i=i0,in
            sum_lcl(j,Ptopo_myrow) = sum_lcl(j,Ptopo_myrow) + F_src(i,j,k)
         end do
         end do
         end do
         if (me_in_row/=num_in_row-1) then
            tag=711+me_in_row
            call MPI_send ( sum_lcl(1,Ptopo_myrow), l_nj, MPI_DOUBLE_PRECISION, me_in_row+1, tag, COMM_row,err )
         endif
      endif
      if (Ptopo_mycol==Ptopo_npex-1) then
         if (Ptopo_npey>1) then
            if (Ptopo_myrow>0) then
               tag=811+Ptopo_myrow
               call MPI_send ( sum_lcl(1,Ptopo_myrow), l_nj, MPI_DOUBLE_PRECISION, 0, tag, COMM_col,err )
            else
               do n=1,Ptopo_npey-1
                  tag=811+n
                  nj=Ptopo_gindx_alongY(2,n+1)-Ptopo_gindx_alongY(1,n+1)+1
                  call MPI_recv ( sum_lcl(1,n), nj, MPI_DOUBLE_PRECISION,n,tag, COMM_col, status, err)
               end do
            endif
         endif
         if ((Ptopo_mycol==Ptopo_npex-1).and.(Ptopo_myrow==0)) then
            do n=0,Ptopo_npey-1
               nj=Ptopo_gindx_alongY(2,n+1)-Ptopo_gindx_alongY(1,n+1)+1
               do j=1,nj
                  F_sum= F_sum + sum_lcl(j,n)
               end do
            end do
         endif
      endif
      proc= Ptopo_colrow(Ptopo_couleur,Ptopo_npex-1,0)
      call MPI_bcast (F_sum, 1, MPI_DOUBLE_PRECISION, proc, COMM_grid, err)
!!$omp end single
!
!----------------------------------------------------------------------
!
      return
      end subroutine glb_orderd_sum
