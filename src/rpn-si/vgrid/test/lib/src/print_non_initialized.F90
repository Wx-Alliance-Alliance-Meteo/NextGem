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
program main
  use vGrid_Descriptors, only : vgrid_descriptor, vgd_print,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  integer :: ier
  type(vgrid_descriptor) :: vgd
  logical :: ok=.false.

  print*,'The following error is expected since vgd is not initiaized'
  ier = vgd_print(vgd)
  if(ier.eq.VGD_ERROR)then
     ok=.true.
     print*,'Problem printing vgd'
  endif

  call ut_report(ok,'Test print non initialized')  

end program main
