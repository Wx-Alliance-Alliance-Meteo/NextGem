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
!================================================================================================
program tests
#include <rmn/msg.h>
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_free, vgd_new,vgd_levels, vgd_putopt
  use Unit_Testing, only: ut_report


  implicit none

  integer,external :: fnom,fstouv,fstfrm,fclos
  integer :: lu,ier,k
  integer, parameter :: nfiles=1
  ! Test files that generate lots of error messages
  character(len=200), dimension(nfiles) :: files
  logical :: ok = .false.
  type(vgrid_descriptor) :: vgd
  
  call msg_verbosity(MSG_DEBUG)
  ! Assign name of files here because in data statements
  ! gfortran will want all files to be the same length in compilation
  files(1)="data/dm_all_sorts_of_toc"
  !files(2)="data/dm_2001_5001_from_editfst"
  !files(3)="data/dm_2001_from_editfst"
  !files(4)="data/dm_1001_just_p0"
  !files(5)="data/dm_1001_from_model_run"

  ier = vgd_putopt("ALLOW_SIGMA",.true.)

  do k=1,nfiles

     lu=0
     ier=fnom(lu,trim(files(k)),"RND+R/O",0)
     print *,'opening file: ',trim(files(k))
     !feeding fnom with 0 forces lu to obtain a number
     if(ier.lt.0)then
        print*,'(Test) ERROR with fnom on file ',trim(files(k))
        error stop 1
     endif
     ier=fstouv(lu,'RND')
     if(ier.lt.0)then
        print*,'(Test) No record in RPN file ',trim(files(k))
        error stop 1
     endif
     
     ! Suppress (Cvgd) and (vgd) messages for vgd_new

     ier = vgd_new(vgd,unit=lu,format="fst",ip1=999,ip2=999, quiet=.true.)
     ier = vgd_new(vgd,unit=lu,format="fst",ip1=999,ip2=-1, quiet=.true.)
     ier = vgd_new(vgd,unit=lu,format="fst",ip1=-1,ip2=999, quiet=.true.)
     ier = vgd_new(vgd,unit=lu,format="fst",ip1=-1,ip2=-1, quiet=.true.)

     ier = vgd_free(vgd)

     ier = fstfrm(lu)
     
     ier=fclos(lu)

  end do

  ok = .true.
  
  call ut_report(ok,message='Grid_Descriptors::vgd_levels level calculation status')
     
end program tests
