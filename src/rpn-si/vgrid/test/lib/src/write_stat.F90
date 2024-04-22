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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat,lui=0,luo=0,lutxt=69,ip1,ip2

  integer, external :: fnom,fstouv,fclos,fstfrm

  ! Open test file
  stat=fnom(lui,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lui,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  stat=fnom(luo,"data_out/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(luo,'RND')
  open(unit=lutxt,file='data/dm_5002_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lui,format="fst",ip1=ip1,ip2=ip2)

  ! Print information about the instance
  stat = vgd_write(d,unit=luo,format='fst')
  call close_fst(lui,luo)

  call ut_report(stat,message='Grid_Descriptors::vgd_write write descriptors')

end program tests

subroutine close_fst(lui,luo)

   implicit none
   integer lui,luo,fstfrm,stat

   stat=fstfrm(lui)
   stat=fstfrm(luo)

   print*,'fst files closed'

end subroutine close_fst
