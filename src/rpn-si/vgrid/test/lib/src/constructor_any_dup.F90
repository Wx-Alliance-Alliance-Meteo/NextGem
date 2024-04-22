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

  ! Objetc: The constructor vgd_new shoul produce an error since there 
  ! are more than one !! matching

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,VGD_ERROR,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: lu=0
  integer :: stat
  integer :: fnom,fstouv,fstfrm
  integer, parameter :: nmax=1000

  stat=fnom(lu,"data/dm_5002_from_model_run-interp_dup","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)

  call ut_report(stat==VGD_ERROR,'Grid_Descriptors, vgd_new, vgd_get CA')

  stat=fstfrm(lu)

end program constructor
