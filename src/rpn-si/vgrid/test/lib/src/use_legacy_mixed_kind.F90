! * libdescrip - Vertical grid descriptor library for FORTRAN programming
! * Copyright (C) 2016  Direction du developpement des previsions nationales
! *                     Centre meteorologique canadien
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
program constructor

   ! Object : try to construct the descriptor but there are many different vertical structures
   !          in the link file, this should produce an error.

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: fnom,fstouv,fstlnk,i
  integer :: stat
  integer, parameter :: nfiles=3
  integer, dimension(nfiles) :: iun
  character(len=500),dimension(nfiles) :: files

  logical :: OK=.false.
  iun=0
  files(1)="data/dm_1001_from_model_run"
  files(2)="data/dm_1002_from_model_run"
  files(3)="data/dm_5001_from_model_run"
  !files(4)="data/dm_5002_from_model_run"

  do i=1,nfiles
     stat=fnom(iun(i),files(i),"RND",0)
     if(stat.lt.0)then
        print*,'ERROR with fnom on ',trim(files(i))
        error stop 1
     endif
     stat=fstouv(iun(i),"RND")
     if(stat.lt.0)then
        print*,'Error with fstouv on ',trim(files(i))
        error stop 1
     endif
  enddo

  stat=fstlnk(iun,nfiles)

  stat=vgd_new(d,iun(1))
  OK=stat.eq.VGD_ERROR

  call ut_report(OK,'Grid_Descriptors::vgd_vintage (5001)')

end program constructor
