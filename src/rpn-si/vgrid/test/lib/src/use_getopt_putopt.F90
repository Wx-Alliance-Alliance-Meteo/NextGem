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
  use vGrid_Descriptors, only: vgd_getopt,vgd_putopt,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  integer :: stat
  logical :: ok=.true.,key_value_L

  !=======
  ! getopt
  !-------
  !
  ! Call with wrong key
  stat=vgd_getopt('THIS_KEY_WILL_NEVER_EXIS',key_value_L)
  print*,'The above error message is nornal'
  stat=vgd_getopt('THIS_KEY_WILL_NEVER_EXIS',key_value_L,quiet=.true.)
  ! The above produces no output

  stat=vgd_getopt('ALLOW_RESHAPE',key_value_L)
  if(stat==VGD_ERROR)ok=.false.
  print*,'ALLOW_RESHAPE=',key_value_L
  if(key_value_L)then
     print*,'Default ALLOW_RESHAPE must by false it is ',key_value_L
     ok=.false.
  endif

  stat=vgd_getopt('ALLOW_SIGMA',key_value_L)
  if(stat==VGD_ERROR)ok=.false.
  print*,'ALLOW_SIGMA=',key_value_L
  if(key_value_L)then
     print*,'Default ALLOW_SIGMA must by false it is ',key_value_L
     ok=.false.
  endif

  !=======
  ! putopt
  !-------
  ! Call with wrong key
  stat=vgd_putopt('THIS_KEY_WILL_NEVER_EXIS',.true.)
  print*,'The above error message is nornal'
  !
  stat=vgd_putopt('ALLOW_RESHAPE',.true.)
  if(stat==VGD_ERROR)then
     print*,'This should not happend please fixit'
     error stop 1
  endif
  stat=vgd_getopt('ALLOW_RESHAPE',key_value_L)
  if(.not.key_value_L)then
     print*,'Key ALLOW_RESHAPE should be .true. and it is not'
     ok=.false.
  endif

  stat=vgd_putopt('ALLOW_SIGMA',.true.)
  if(stat==VGD_ERROR)then
     print*,'This should not happend please fixit 2'
     error stop 1
  endif
  stat=vgd_getopt('ALLOW_SIGMA',key_value_L)
  if(.not.key_value_L)then
     print*,'Key ALLOW_SIGMA should be .true. and it is not'
     ok=.false.
  endif
  
  call ut_report(ok,message='Grid_Descriptors:: vgd_getopt vgd_putopt')

end program tests
