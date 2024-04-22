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

!**s/r periodicX

      subroutine periodicX ()
      use mem_tracers
      use gmm_vt0
      use glb_ld
      use tr3d
      use ptopo
      implicit none

      include 'mpif.h'
      integer i,j,k,n,cnt,tag,err,req
      real, dimension(max(pil_w,pil_e)*l_nj*l_nk*(6+Tr3d_ntr)) :: buf1,buf2
!     
!----------------------------------------------------------------------
!     
!!$omp single
      if (Ptopo_npex==1) then
         
         do n=1,Tr3d_ntr
            do k=1,G_nk
               do j=1,l_nj
                  do i=1,pil_w
                     tracers_M(n)%pntr(i,j,k)= tracers_M(n)%pntr(l_ni-2*pil_e+i,j,k)
                  end do
                  do i=l_ni-pil_e+1,l_ni
                     tracers_M(n)%pntr(i,j,k)= tracers_M(n)%pntr(i-l_ni+2*pil_w,j,k) 
                  end do
               end do
            end do
         end do
         do k=1,G_nk
            do j=1,l_nj
               do i=1,pil_w
                  ut0(i,j,k) = ut0(l_ni-2*pil_e+i-1,j,k)
                  vt0(i,j,k) = vt0(l_ni-2*pil_e+i,j,k)
                  tt0(i,j,k) = tt0(l_ni-2*pil_e+i,j,k)
                  wt0(i,j,k) = wt0(l_ni-2*pil_e+i,j,k)
                  zdt0(i,j,k)= zdt0(l_ni-2*pil_e+i,j,k)
                  qt0(i,j,k) = qt0(l_ni-2*pil_e+i,j,k)
               end do
               do i=l_ni-pil_e+1,l_ni
                  ut0(i-1,j,k) = ut0(i-l_ni+2*pil_w,j,k) 
                  vt0(i,j,k) = vt0(i-l_ni+2*pil_w,j,k) 
                  tt0(i,j,k) = tt0(i-l_ni+2*pil_w,j,k) 
                  wt0(i,j,k) = wt0(i-l_ni+2*pil_w,j,k) 
                  zdt0(i,j,k)= zdt0(i-l_ni+2*pil_w,j,k) 
                  qt0(i,j,k) = qt0(i-l_ni+2*pil_w,j,k) 
               end do
            end do
         end do

      else

         if (Ptopo_mycol==0) then
            cnt=0
            do n=1,Tr3d_ntr
               do k=1,G_nk
                  do j=1,l_nj
                     do i=1+pil_w,2*pil_w
                        cnt=cnt+1
                        buf1(cnt)= tracers_M(n)%pntr(i,j,k)
                     end do
                  end do
               end do
            end do
            do k=1,G_nk
               do j=1,l_nj
                  do i=1+pil_w,2*pil_w
                     cnt=cnt+1 ; buf1(cnt)= ut0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= vt0(i,j,k) 
                     cnt=cnt+1 ; buf1(cnt)= tt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= wt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)=zdt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= qt0(i,j,k)
                  end do
               end do
            end do
            tag= 711
            call MPI_isend ( buf1, size(buf1), MPI_REAL, Ptopo_npex-1, tag, COMM_row, req, err )
            tag= 712
            call MPI_recv ( buf2, size(buf2), MPI_REAL, Ptopo_npex-1, tag, COMM_row, MPI_STATUSES_IGNORE, err)

            cnt=0
            do n=1,Tr3d_ntr
               do k=1,G_nk
                  do j=1,l_nj
                     do i=1,pil_w
                        cnt=cnt+1
                        tracers_M(n)%pntr(i,j,k)= buf2(cnt)
                     end do
                  end do
               end do
            end do
            do k=1,G_nk
               do j=1,l_nj
                  do i=1,pil_w
                     cnt=cnt+1 ; ut0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; vt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; tt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; wt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ;zdt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; qt0(i,j,k)= buf2(cnt)
                  end do
               end do
            end do
            call MPI_waitall (1,req,MPI_STATUSES_IGNORE,err)
            
         else if (Ptopo_mycol==Ptopo_npex-1) then

            cnt=0
            do n=1,Tr3d_ntr
               do k=1,G_nk
                  do j=1,l_nj
                     do i=l_ni-2*pil_e+1,l_ni-pil_e
                        cnt=cnt+1
                        buf1(cnt)= tracers_M(n)%pntr(i,j,k)
                     end do
                  end do
               end do
            end do
            do k=1,G_nk
               do j=1,l_nj
                  do i=l_ni-2*pil_e+1,l_ni-pil_e
                     cnt=cnt+1 ; buf1(cnt)= ut0(i-1,j,k)
                     cnt=cnt+1 ; buf1(cnt)= vt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= tt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= wt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)=zdt0(i,j,k)
                     cnt=cnt+1 ; buf1(cnt)= qt0(i,j,k)
                  end do
               end do
            end do
            tag= 712
            call MPI_isend ( buf1, size(buf1), MPI_REAL, 0, tag, COMM_row, req, err )
            tag= 711
            call MPI_recv ( buf2, size(buf2), MPI_REAL, 0, tag, COMM_row, MPI_STATUSES_IGNORE, err)

            cnt=0
            do n=1,Tr3d_ntr
               do k=1,G_nk
                  do j=1,l_nj
                     do i=l_ni-pil_e+1,l_ni
                        cnt=cnt+1
                        tracers_M(n)%pntr(i,j,k)= buf2(cnt)
                     end do
                  end do
               end do
            end do
            do k=1,G_nk
               do j=1,l_nj
                  do i=l_ni-pil_e+1,l_ni
                     cnt=cnt+1 ; ut0(i-1,j,k)= buf2(cnt)
                     cnt=cnt+1 ; vt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; tt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; wt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ;zdt0(i,j,k)= buf2(cnt)
                     cnt=cnt+1 ; qt0(i,j,k)= buf2(cnt)
                  end do
               end do
            end do
            call MPI_waitall (1,req,MPI_STATUSES_IGNORE,err)
         endif
      endif
!!$omp end single
!     
!----------------------------------------------------------------------
!     
      return
      end subroutine periodicX
