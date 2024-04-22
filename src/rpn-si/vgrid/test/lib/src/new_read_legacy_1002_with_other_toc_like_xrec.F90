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
!================================================================================================
program tests
#include <rmn/msg.h>
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_free, vgd_new,vgd_levels,operator(==), VGD_OK
  use Unit_Testing, only: ut_report


  implicit none

  integer :: fnom,fstouv,fclos,fstlnk,ier,i
  logical :: ok = .false.
  integer, dimension(3) :: list=(/10,11,12/)
  type(vgrid_descriptor) :: vgd1, vgd2, vgd3
  character(len=72), dimension(3) :: file_S=(/"data/dm_1002_from_model_run       ","data/dm_5005_from_model_run       ",&
       "data/dm_21001_from_model_run_SLEVE"/)

  call msg_verbosity(MSG_DEBUG)

  ier=fnom(list(1),file_S(1),"RND+R/O",0)
  if(ier.lt.0)then
     print*,'(Test) ERROR with fnom on file ',trim(file_S(1))
     error stop 1
  endif
  ier=fstouv(list(1),'RND')
  if(ier.lt.0)then
     print*,'(Test) No record in RPN file ',trim(file_S(1))
     error stop 1
  endif

  ! First construct the Vcode 1002 (eta) !! from legacy encoding from
  ! a pure file having only P0, PT and TT, UU ... in it of Vcode 1002
  if(VGD_OK /= vgd_new(vgd1,unit=list(1),kind=1))then
     print*,"Error could not read vgrid descriptor 1002"
     error stop 1
  endif
  
  ! Now simulater what xrec or other applycation may be doing be joining other files
  ! with other !! with Vcode different from 1002 and try to vgd_new on kind 2 like above
  do i=2,3
     ier=fnom(list(i),file_S(i),"RND+R/O",0)
     if(ier.lt.0)then
        print*,'(Test) ERROR with fnom on file ',trim(file_S(i))
        error stop 1
     endif
     ier=fstouv(list(i),'RND')
     if(ier.lt.0)then
        print*,'(Test) No record in RPN file ',trim(file_S(i))
        error stop 1
     endif
  end do

  ier = fstlnk(list,3)
  ! Prove that we can read a kind 5 toc
  if(VGD_OK /= vgd_new(vgd2,unit=list(1),kind=5))then
     print*,"Error could not read vgrid descriptor of kind 5"
     error stop 1
  endif  
  ! Prove that we can still reconstruc 1002 from legacy even with other !! and records of diferent kind
   ! Prove that we can read a kind 5 toc
  if(VGD_OK /= vgd_new(vgd3,unit=list(1),kind=1))then
     print*,"Error could not read vgrid descriptor 1002 on file with other !!"
     error stop 1
  endif
  if (.not. vgd1 == vgd3) then
     print*,"Error vgrid from legacy constructed from file with 1002 and other vcode is not equal as pure legacy construction"
     error stop 1
  endif
  
  do i=1,3
     ier=fclos(list(i))
  end do

  ok = .true.
  
  call ut_report(ok,message='Grid_Descriptors::new read legace 1002 with other !!')
     
end program tests
