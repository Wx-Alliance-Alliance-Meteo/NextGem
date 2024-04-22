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
!----------------------------------LICENCE END ---------------------------------
      subroutine gem_error (F_errorCode, F_FromSubName, F_Message)
      use iso_c_binding
      use lun
      use ptopo
      implicit none

      include 'mpif.h'
      
      integer :: F_errorCode
      character(len=*) :: F_FromSubName
      character(len=*) :: F_Message

      integer :: errcode, err
!
!     ---------------------------------------------------------------
!
      call MPI_allreduce (F_errorCode, errcode,1,MPI_INTEGER,&
                               MPI_MIN,COMM_MULTIGRID,err)
      if (errcode < 0) then
         if (Lun_out > 0) write(Lun_out,2000) F_FromSubName, F_Message
         call gem_stop ()
      endif

2000  format (/60('#')/2x,a': ',a/2x,'ABORT'/60('#')/)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine gem_error
      
      subroutine gem_error_omp (F_errorCode, F_FromSubName, F_Message)
      use iso_c_binding
      use lun
      use ptopo
      implicit none

      include 'mpif.h'
      
      integer :: F_errorCode
      character(len=*) :: F_FromSubName
      character(len=*) :: F_Message

      integer :: errcode, err
!
!     ---------------------------------------------------------------
!
!!$omp single
      call MPI_allreduce (F_errorCode, errcode,1,MPI_INTEGER,&
                               MPI_MIN,COMM_MULTIGRID,err)
      if (errcode < 0) then
         if (Lun_out > 0) write(Lun_out,2000) F_FromSubName, F_Message
      endif
!!$omp end single copyprivate(errcode)

      if (errcode < 0) call gem_stop ()

2000  format (/60('#')/2x,a': ',a/2x,'ABORT'/60('#')/)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine gem_error_omp
