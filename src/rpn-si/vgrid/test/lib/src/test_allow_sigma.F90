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
program test_allow_sigma

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_putopt,vgd_getopt,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: lu=0
  integer :: stat
  integer :: fnom,fstouv,fstfrm
  logical :: ok=.true., my_bol_value

  stat=fnom(lu,"data/dm_1001_from_model_run","RND+R/O",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Try to get the value of ALLOW_SIGMA. it should be .false. by default
  stat = vgd_getopt("ALLOW_SIGMA",my_bol_value)
  print*,'==========================================================================='
  if(stat == VGD_ERROR)then
     print*, "ERROR in test, cannot get ALLOW_RESHAPE value"
     error stop 1
  endif
  print*,'my_bol_value=',my_bol_value
  if(my_bol_value)then
      print*, "ERROR in test, ALLOW_SIGMA should be .false."
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  print*,'==========================================================================='
  print*,'The following error message telling that "sigma coordinate construction is not ALLOWED" is normal'
  stat = vgd_new(d,unit=lu,format="fst")
  if(stat/=VGD_ERROR)then
     print*,'ERROR in test, should not be allowed to construct a new sigma coordinate since ALLOW_SIGMA should be false by default'
     ok=.false.
  endif
  print*,'============================================================================'
  stat = vgd_putopt("ALLOW_SIGMA",.true.)
  if(stat == VGD_ERROR)then
     print*, "ERROR in test, cannot put option ALLOW_SIGMA"
     error stop 1
  endif
  stat = vgd_new(d,unit=lu,format="fst")
  if(stat==VGD_ERROR)then
     print*,'ERROR in test, with vgd_new. Constructing sigma should be possible.'
     ok=.false.
  endif
  
  call ut_report(stat,'Grid_Descriptors, vgd_new, on 1001')

  stat=fstfrm(lu)

end program test_allow_sigma
