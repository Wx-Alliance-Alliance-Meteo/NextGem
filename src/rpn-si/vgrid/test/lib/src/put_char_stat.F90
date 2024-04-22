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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_put,VGD_LEN_NAME,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat,fnom,fstouv,lu=10
  logical :: ok = .true.

  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,lu,format='fst')
  if(stat == VGD_ERROR)then
     ok = .false.
     print*,'Error with vgd_put'
  end if

  ! Change an element of the structure
  stat = vgd_put(d,key='ETIKE',value='123456789123')
  if(stat == VGD_ERROR)then
     ok = .false.
     print*,'Error with vgd_put'
  end if

  call ut_report(ok,message='Grid_Descriptors::vgd_put put character value')

end program tests
