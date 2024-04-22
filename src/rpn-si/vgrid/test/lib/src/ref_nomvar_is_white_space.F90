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
program ref_nomvar_is_white_space

   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,VGD_LEN_RFLD,VGD_LEN_RFLS, VGD_OK, VGD_ERROR
   use Unit_Testing, only: ut_report

   implicit none
   
   type(vgrid_descriptor) :: vgd
   integer :: lu=10, stat
   integer :: fnom,fstouv,fstfrm
   character(len=VGD_LEN_RFLD) :: rfld
   character(len=VGD_LEN_RFLS) :: rfls
   logical :: ok=.true.

   if( fnom(lu,"data/dm_2001_from_editfst","RND+R/O",0) < 0 )then
      print*,'ERROR with fnom'
      error stop 1
   endif
   if( fstouv(lu,'RND') < 0)then 
      print*,'No record in RPN file'
      error stop 1
   endif
   if( vgd_new(vgd,unit=lu,format="fst") == VGD_ERROR) error stop 1

   if( vgd_get(vgd,"RFLD", rfld, .true.) == VGD_OK )then
      print*,"RFLD='",rfls,"'"
      print*,'In test, problem with vgd_get on "RFLD", should have returned an error but returned VGD_OK'
      ok = .false.
   endif

   if( vgd_get(vgd,"RFLS", rfls, .true.) == VGD_OK )then
      print*,"RFLS='",rfls,"'"
      print*,'In test, problem with vgd_get on "RFLS", should have returned an error but returned VGD_OK'
      ok = .false.
   endif
   
   stat = fstfrm(lu)

   call ut_report(ok,'For Vcode without ref fields nomvar must be "    "')

end program ref_nomvar_is_white_space
