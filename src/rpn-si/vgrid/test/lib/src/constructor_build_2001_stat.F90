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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,VGD_ERROR,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat
  real(kind=8), dimension(4) :: pres=(/1000.,925.,850.,700./)

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,kind=2,version=1,nk=size(pres),ip1_m=int(pres),a_m_8=pres,b_m_8=pres)
  if(stat.eq.VGD_ERROR)then    
     print*,'UNEXPECTED problem with vgd_new'
     error stop 1
  endif

  call ut_report(stat,'Grid_Descriptors::vgd_new vertical build initializer (2001)')

end program constructor
