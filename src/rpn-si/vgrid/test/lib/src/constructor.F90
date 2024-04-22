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

  ! Object: construct vgrid_descriptor

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new, VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: lu=10,lutxt=69
  integer :: stat,ip1,ip2
  integer :: fnom,fstouv,fstfrm
  logical :: ok

  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then 
     print*,'No record in RPN file'
     error stop 1
  endif
  open(unit=lutxt,file='data/dm_5002_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Construct a new set of 3D coordinate descriptors
  print*,'The following error on optional ip2 is expected'
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1)
  if(stat.eq.VGD_OK)ok=.false.

  print*,'The following error on optional ip1 is expected'
  stat = vgd_new(d,unit=lu,format="fst",ip2=ip2)
  if(stat.eq.VGD_OK)ok=.false.

  print*,'The following error on format is expected'
  stat = vgd_new(d,unit=lu,format="toto")
  if(stat.eq.VGD_OK)ok=.false.  

  print*,'The following error matching ip1 ip2 is expected'
  stat = vgd_new(d,unit=lu,format="fst",ip1=123,ip2=456)
  if(stat.eq.VGD_OK)ok=.false.  
  
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)
  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_get CA')

  stat=fstfrm(lu)

end program constructor
