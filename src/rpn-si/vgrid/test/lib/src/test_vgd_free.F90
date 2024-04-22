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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_getopt,vgd_putopt,vgd_print,vgd_free,vgd_get,VGD_ERROR,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer :: stat
  logical :: ok=.true.

  ! Structure elements
  real, dimension(3) :: hyb=(/0.1,  0.5, 0.9/)
  real :: rcoef1=0.,rcoef2=1.  
  real(kind=8) :: ptop_8=1000d0,pref_8=100000d0,value_8
  
  print*,'First build';flush(6)
  stat = vgd_new(vgd,kind=5,version=2,hyb=hyb,rcoef1=rcoef1,rcoef2=rcoef2,ptop_8=ptop_8,pref_8=pref_8)
  if(stat==VGD_ERROR)then
     print*,'This vgd_new error should not happen, please fixit'
     error stop 1
  end if  
  stat = vgd_get(vgd,'PREF - reference pressure',value_8)
  if(stat==VGD_ERROR)then
     print*,'This vgd_get error should not happen, please fixit'
     error stop 1
  end if  
  if(abs(value_8-pref_8)/pref_8 > epsilon(pref_8))then
     print*,'This error on value_8 should not happen, please fixit'
     error stop 1    
  endif

  stat = vgd_free(vgd)
  if(stat==VGD_ERROR)then
     print*,'This vgd_free error should not happen, please fixit'
     error stop 1
  end if  

  stat = vgd_get(vgd,'PREF - reference pressure',value_8)
  if(stat==VGD_ERROR)then
     print*,'The above error in normal'
  else
     print*,'The above call to vgd_get should produce an error and does not, please fixit'
     error stop 1    
  end if

  call ut_report(ok,message='Grid_Descriptors:: test allow reshape')

end program tests
