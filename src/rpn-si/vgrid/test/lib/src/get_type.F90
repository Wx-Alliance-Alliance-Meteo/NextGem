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
program get_put_get_all

  ! Goal, test vgd_put on all Vcode
  use vGrid_Descriptors, only: vgd_putopt, VGD_OK, VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none
  integer :: i, test_it
  integer, parameter :: nfiles=15
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_2001_from_editfst             ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5003_from_model_run           ",&
       "data/dm_5004_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_5999_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)
  
  if( vgd_putopt("ALLOW_SIGMA",.true.) == VGD_ERROR )then
     print*,'Error with vgd_putopt on ALLOW_SIGM'
     error stop 1
  endif
  
  do i=1, nfiles
     if( test_it(files(i),i) == VGD_ERROR )then
        error stop 1
     endif
  enddo

  call ut_report(VGD_OK,'Grid_Descriptors, vgd_put')

end program get_put_get_all

!=======================================================================
!=======================================================================
!=======================================================================

integer function test_it(F_file, ind) result(status)
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_new, vgd_get, &
       vgd_free, VGD_OK, VGD_ERROR, VGD_PRES_TYPE, VGD_HEIGHT_TYPE
   implicit none
   integer, intent(in) :: ind
   character(len=*), intent(in) :: F_file
   ! Local variables
   integer :: lu, fnom, fstouv, fstfrm, fclos, ier
   integer :: vcode, type
   type (vgrid_descriptor) :: vgd

   status = VGD_ERROR
   
   lu = 10 + ind
   print*,'Treating file ',trim(F_file)
   if( fnom(lu, F_file, 'RND+R/O', 0) < 0 )then
      print*,'ERROR: with fnom on ', trim(F_file)
      return
   endif
   if( fstouv(lu, 'RND') == 0 )then
      print*,'ERROR: no record in file ', trim(F_file)
      return
   endif   
   if( vgd_new(vgd, lu, 'fst') == VGD_ERROR )then
      print*,'Error with vgd_new on file ', trim(F_file)
      error stop 1
   endif
   
   if( vgd_get(vgd,"TYPE",type) == VGD_ERROR )then
      print*,'TEST ERROR: cannot get key TYPE for file',trim(F_file)
      error stop 1
   endif

   ier = vgd_get(vgd,"VCOD",vcode)

   select case (vcode)
   case (1001,1002,1003,2001,5001,5002,5003,5004,5005,5100,5999)
      if(type /= VGD_PRES_TYPE)then
         print*,'TEST ERROR: vcode',vcode,' should be of type pressure and it is not'
         error stop 1
      endif
   case (1,4001,21001,21002)
      if(type /= VGD_HEIGHT_TYPE)then
         print*,'TEST ERROR: vcode',vcode,' should be of type height and it is not'
         error stop 1
      endif
   case DEFAULT
      print*,'TEST ERROR: vcode',vcode,' not supported in test get_type, please add in code'
      error stop 1
 end select
   
   ier = vgd_free(vgd)
   status = fstfrm(lu)
   status = fclos(lu)
   status = VGD_OK
   
 end function test_it
