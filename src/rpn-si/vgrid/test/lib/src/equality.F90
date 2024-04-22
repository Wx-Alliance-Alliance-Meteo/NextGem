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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_free,operator(==),vgd_putopt,VGD_OK,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: vgrid,vgrid_rebuilt
  integer :: lu=10,i,ier
  integer :: stat
  integer :: fnom,fstouv,fstfrm

  integer, parameter :: nfiles=14
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5003_from_model_run           ",&
       "data/dm_5004_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)

  stat = VGD_OK

  ier = vgd_putopt("ALLOW_SIGMA",.true.)

  do i=1,nfiles
     print*,'==============================================='
     print*,'FILE = ',trim(files(i))
     ier=fnom(lu+i,files(i),"RND",0)     
     if(ier.lt.0)then
        print*,'ERROR with fnom on file ',files(i)
        error stop 1
     endif
     ier=fstouv(lu+i,'RND')
     if(ier.le.0)then
        print*,'No record in RPN file'
        error stop 1
     endif
     ! Construct a new set of 3D coordinate descriptors
     ier = vgd_new(vgrid        ,unit=lu+i,format="fst")
     ier = vgd_new(vgrid_rebuilt,unit=lu+i,format="fst")
     if (.not. vgrid == vgrid_rebuilt) then
        stat = VGD_ERROR
     endif
     ier = vgd_free(vgrid)
     ier = vgd_free(vgrid_rebuilt)
     ier = fstfrm(lu+i)

  enddo

  print*,'**********************************************************'
  if(stat == VGD_OK)then
     print*,'ALL TESTS ON EQUALITY PASSED'
  else
     print*,'SOME OR ALL TESTS ON EQUALITY FAILED'
  endif
  print*,'**********************************************************'

  call ut_report(stat,'Grid_Descriptors, equality')

end program constructor
