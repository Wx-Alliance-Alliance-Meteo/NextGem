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
program defaults

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_putopt,VGD_MISSING,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: vgrid
  integer :: lu=0
  integer :: stat
  integer :: fnom,fstouv,fstfrm
  real(kind=8) :: ptop
  
  stat=fnom(lu,"data/dm_1001_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_putopt("ALLOW_SIGMA",.true.)
  stat = vgd_new(vgrid,unit=lu,format="fst")
  if(stat.eq.VGD_ERROR)then
     print*,'In test error with vgd_new'
     error stop 1
  endif

  stat = vgd_get(vgrid,'PTOP',ptop)
  if (abs(ptop-dble(VGD_MISSING))<10.*epsilon(ptop)) then
     stat = 0
  else
     stat = -1
  endif
  call ut_report(stat,'Grid_Descriptors, rebuild from table')

  stat=fstfrm(lu)

end program defaults
