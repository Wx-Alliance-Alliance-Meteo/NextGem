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
  use vGrid_Descriptors, only: VGD_ERROR,vgd_putopt
  use Unit_Testing, only: ut_report
  
  implicit none
  
  integer, parameter :: nfiles=14
  integer :: i,test_it
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
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
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE"&
       /)
  logical :: OK=.true.

  if(vgd_putopt("ALLOW_SIGMA",.true.) == VGD_ERROR)OK=.false.
  
  do i=1,nfiles
     if(test_it(files(i),i) == VGD_ERROR)OK=.false.
  end do

  call ut_report(OK,'Grid_Descriptors, vgd_new, vgd_get CA')
   
end program constructor

integer function test_it(F_file,F_index) result(istat)
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_new,vgd_print,vgd_write,VGD_OK,VGD_ERROR,operator(==),vgd_put,vgd_print
   implicit none

   integer :: F_index
   character(len=*) :: F_file
   
   ! Local variables
   
   type(vgrid_descriptor) :: vgd,vgd2
   integer :: stat, kind,version,lu1,lu2,lu3
   integer :: fnom,fstouv,fstfrm
   character(len=128) :: my_file

   istat = VGD_ERROR
   
   lu1=10+F_index*3
   lu2=lu1+1
   lu3=lu2+1
   
   print*,"=========================================================="
   print*,'TESTING ',trim(F_file)

   stat=fnom(lu1,F_file,"RND+R/O",0)
   if(stat.lt.0)then
      print*,'ERROR with fnom on ',trim(F_file)
      return
   endif
   stat=fstouv(lu1,'RND')
   if(stat.le.0)then
      print*,'No record in RPN file ',trim(F_file)
      error stop 1
   endif
   
   ! Construct a new set of 3D coordinate descriptors
   if(vgd_new(vgd,unit=lu1,format="fst") == VGD_ERROR)return
   stat=fstfrm(lu1)

   if(vgd_get(vgd,key="KIND",value=kind)    == VGD_ERROR)return
   if(vgd_get(vgd,key="VERS",value=version) == VGD_ERROR)return
   write(my_file,*)kind*1000+version   
   my_file=adjustl(my_file)
   my_file="data_out/toctoc_"//trim(my_file)//".fst"

   stat=fnom(lu2,my_file,"RND",0)
   if(stat.lt.0)then
      print*,'ERROR on fnom with ',my_file
      error stop 1
   endif
   stat=fstouv(lu2,'RND')   
   if(stat.lt.0)then
      print*,'ERROR on fstouv with ',my_file
      error stop 1
   endif
   if(vgd_write(vgd,unit=lu2,format="fst") == VGD_ERROR)return
   stat=fstfrm(lu2)
   
   stat=fnom(lu3,my_file,"RND+R/O",0)
   if(stat.lt.0)then
      print*,'ERROR on fnom with reopenning ',my_file
      error stop 1
   endif
   stat=fstouv(lu3,'RND')   
   if(stat.lt.0)then
      print*,'ERROR on fstouv with reopenning ',my_file
      error stop 1
   endif
    
   if(vgd_new(vgd2,unit=lu3,format="fst") == VGD_ERROR)return
   stat=fstfrm(lu3)
      
   if (.not. (vgd2 == vgd) )then
      print*,'In test read_write_read_all, !! do not match'
      stat=vgd_print(vgd2)
      return
   else
      print*,'!! matched'
   endif

   istat=VGD_OK
   
end function test_it

