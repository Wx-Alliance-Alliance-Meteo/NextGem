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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_print,vgd_putopt,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: i, stat
  integer :: fnom,fstouv,fstfrm
  logical :: ok = .true.
  integer, parameter :: nfiles=14
  character(len=200), dimension(nfiles) :: files=(/&
       "data/gg_5001_from_model_run_plus_toc  ",&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_2001_from_editfst             ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_5999_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)

  stat = vgd_putopt("ALLOW_SIGMA",.true.)

  do i=1,nfiles
     if( fnom(lu+i,files(i),"RND",0) < 0)then
        print*,'ERROR with fnom, on file',trim(files(i))
        error stop 1
     end if
     if( fstouv(lu+i,'RND') <= 0 )then
        print*,'No record in RPN file ',trim(files(i))
        error stop 1
     endif
     ! Construct a new set of 3D coordinate descriptors
     if( vgd_new(d,unit=lu+i,format="fst") == VGD_ERROR )then
        print*,'ERROR: problem with vgd_new on file ',trim(files(i))
        error stop 1
     endif
     if( vgd_print(d) == VGD_ERROR) ok = .false.
     stat=fstfrm(lu+i)
  enddo
  call ut_report(ok,'Grid_Descriptors, vgd_new, vgd_print')

end program constructor
