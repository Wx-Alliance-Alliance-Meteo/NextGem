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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write, VGD_OK, VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10,lu2=20,lutxt=69
  integer :: stat,ip1,ip2
  integer :: fnom,fstouv,fstfrm
  real(kind=8), dimension(:), pointer :: my_a
  logical :: ok=.true.

  nullify(my_a)

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

  stat=fnom(lu2,"data_out/toc_toc.rpn","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom toc toc'
     error stop 1
  endif
  stat=fstouv(lu2,'RND')

  ! Try to write unconstructed
  print*,'The following error on vgrid descriptor not constructed is normal '
  stat = vgd_write(d,unit=lu2,format="fst")
  if(stat /= VGD_ERROR) ok=.false.

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)
  if(stat /= VGD_OK)ok=.false.
  stat = vgd_write(d,unit=lu2,format="fst")
  if(stat /= VGD_OK)ok=.false.

  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_get CA')

  stat=fstfrm(lu)
  stat=fstfrm(lu2)

end program constructor
