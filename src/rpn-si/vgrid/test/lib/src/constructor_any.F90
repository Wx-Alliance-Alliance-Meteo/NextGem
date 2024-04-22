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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,operator(==),VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d,d2
  integer :: lu=0,lu2=0,lutxt=69
  integer :: stat,ip1,ip2
  integer :: fnom,fstouv,fstfrm,fclos

  flush(6)
 
  stat=fnom(lu,"data/dm_5002_from_model_run-interp","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  stat=fnom(lu2,"output.fst","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu2,'RND')
  open(unit=lutxt,file='data/dm_5002_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)
  stat = vgd_write(d,unit=lu2,format="fst")
  stat = vgd_new(d2,unit=lu2,format="fst",ip1=ip1,ip2=ip2)
  if(.not.d2==d)stat=vgd_ERROR

  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_get CA')

  stat=fstfrm(lu)
  stat=fstfrm(lu2)
  stat=fclos(lu)
  stat=fclos(lu2)
  call system("rm -f output.fst")

end program constructor
