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
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_new,Vgd_get,Vgd_print,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: stat
  integer :: fnom,fstouv,fstfrm,fclos,kind,version
  logical :: ok=.true.

  ! 1002
  stat=fnom(lu,"data/dm_1002_from_model_run_plus_toc","RND",0)
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
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_new'
     stat=fstfrm(lu)
     error stop 1
  endif
  stat = vgd_print(d)
  stat = vgd_get(d,key='KIND',value=kind)
  if(stat.ne.VGD_OK)ok=.false.
  stat = vgd_get(d,key='VERS',value=version)
  if(stat.ne.VGD_OK)ok=.false.
  print*,'1002 if of kind ',kind,' and version', version
    if(kind*1000+version.ne.1002)ok=.false.

  stat=fstfrm(lu)  
  stat=fclos(lu)

  call ut_report(ok,'Grid_Descriptors, vgd_get, LOGP 5002')  

  stat=fstfrm(lu)

end program constructor
