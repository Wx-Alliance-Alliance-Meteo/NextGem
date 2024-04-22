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
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_new,Vgd_print,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: stat
  integer :: fnom,fstouv,fstfrm
  logical :: ok=.true.

  stat=fnom(lu,"data/dm_all_sorts_of_toc","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  print*,'=========================================================================='
  print*,'0) The folling should produce an error because there are many !! in file' 
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)
  if(stat.eq.VGD_OK)then
     ok=.false.
     print*,'Test failed'
  endif

  print*,'=========================================================================='
  print*,'1) The folling should produce an error since there are more than one !!'
  print*,'satisfying the selection'
  stat = vgd_new(d,unit=lu,format="fst",ip1=2009,ip2=1000)
  if(stat.eq.VGD_OK)then
     ok=.false.
     print*,'Test failed'
  endif

  print*,'=========================================================================='
  print*,'2) The folling should produce an error since there more than one !! with '
  print*,'kind 1 but having a different version'
  stat = vgd_new(d,unit=lu,format="fst",ip1=2009,ip2=1000,kind=1)
  if(stat.eq.VGD_OK)then
     ok=.false.
     print*,'Test failed'
  endif     

  print*,'=========================================================================='
  print*,'3) The folling should NOT produce an error since there is only one !!'
  stat = vgd_new(d,unit=lu,format="fst",ip1=2009,ip2=1000,kind=1,version=1)
  if(stat.ne.VGD_OK)then
     ok=.false.
     print*,'Test failed'
  endif

  print*,'=========================================================================='
  print*,'4) The folling should NOT produce an error since there is only one !!'
  stat = vgd_new(d,unit=lu,format="fst",kind=2)
  if(stat.ne.VGD_OK)then
     ok=.false.  
     print*,'Test failed'
  endif

  print*,'=========================================================================='
  print*,'5) The folling should produce an error since option kind must be used with option version'
  stat = vgd_new(d,unit=lu,format="fst",version=1)
  if(stat.eq.VGD_OK)then
     ok=.false.  
     print*,'Test failed'
  endif

  print*,'=========================================================================='
  print*,'6) The folling should produce an error since vcode doesnt exist'
  stat = vgd_new(d,unit=lu,format="fst",kind=5,version=10)
  if(stat.eq.VGD_OK)then
     ok=.false.  
     print*,'Test failed'
  endif

  call ut_report(ok,'Grid_Descriptors, deal_with_two_toc_toc_same_ver_put_diff_ips vgd_new,  5002')

  stat=fstfrm(lu)

end program constructor
