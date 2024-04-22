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
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_new,Vgd_print,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: stat
  integer :: fnom,fstouv,fstfrm

  stat=fnom(lu,"data/dm_5002_toc_dyn_toc_phy","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  !if(stat.le.0)then
  !   print*,'No record in RPN file'
  !   error stop 1
  !endif

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_new'
     stat=fstfrm(lu)
     error stop 1
  endif
  stat = vgd_print(d)

  call ut_report(stat,'Grid_Descriptors, deal_with_two_toc_toc_same_ver_put_diff_ips vgd_new,  5002')

  stat=fstfrm(lu)

end program constructor
