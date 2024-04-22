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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer, parameter :: nk=27
  integer :: stat
  integer, dimension(:), pointer :: vipm,vipt,work_i
  real, dimension(nk) :: sigma
  real, dimension(:), pointer :: vcdm,vcdt,work
  real(kind=8), dimension(:), pointer :: b_m_8,a_m_8,b_t_8,a_t_8,work_8
  logical :: OK=.true.
  logical, parameter :: write_control_L=.false.

  nullify(vipm,vipt,work_i,vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8)

  sigma=(/0.011000, 0.027000, 0.051000, 0.075000, 0.101000, 0.127000, 0.155000,&
 0.185000, 0.219000, 0.258000, 0.302000, 0.351000, 0.405000, 0.460000, 0.516000,&
 0.574000, 0.631000, 0.688000, 0.744000, 0.796000, 0.842000, 0.884000, 0.922000,&
 0.955000, 0.980000, 0.993000, 1.000000/)

  stat = vgd_new(vgd,kind=1,version=1,hyb=sigma)
  
  stat = vgd_print(vgd)

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer (1001) value')

end program constructor
