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
program tests
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_put,VGD_ERROR, VGD_OK
  use Unit_Testing, only: ut_report
  
  
  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat,lu=10,lutxt=69,ip1,ip2
  logical :: OK = .true.
  
  integer, external :: fnom,fstouv,fclos,fstfrm

  stat=fnom(lu,"data/dm_5005_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  open(unit=lutxt,file='data/dm_5005_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  print*,'The following error vgrid descriptor not constructed is normal'
  stat = vgd_put(d,key='TOTO',value=1)
  if(stat /= VGD_ERROR) OK=.false.

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)
  if(stat /= VGD_OK) OK=.false.

  print*,'The following error "Cvgd_put_int, invalid vgrid" is normal'
  stat = vgd_put(d,key='TOTO',value=1)
  if(stat /= VGD_ERROR) OK=.false.

  print*,'The following error "Cvgd_put_char, invalid key" is normal'
  stat = vgd_put(d,key='TOTO',value="ABC")
  if(stat /= VGD_ERROR) OK=.false.

  call ut_report(OK,message='Grid_Descriptors::vgd_get get valid value')

  ! Close files
  stat = fstfrm(lu)
  stat = fclos(lu)

end program tests
