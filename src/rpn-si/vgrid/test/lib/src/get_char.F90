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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,VGD_LEN_NAME
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat,lu=0,fnom,fstouv,fstfrm,fclos,lutxt=69,ip1,ip2
  character(len=VGD_LEN_NAME) :: rfld
  logical :: ok=.true.

  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  open(unit=lutxt,file='data/dm_5002_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)

  ! Get information about the coordinate
  stat = vgd_get(d,key='RFLD - reference field name',value=rfld)
  if(rfld/='P0')ok=.false.
  stat = vgd_get(d,key='RFLD - reference field name',value=rfld,quiet=.true.)
  if(rfld/='P0')ok=.false.
  stat = vgd_get(d,key='RFLD - reference field name',value=rfld,quiet=.false.)
  if(rfld/='P0')ok=.false.

  call ut_report(ok,message='Grid_Descriptors::vgd_get get character value')

  stat = fstfrm(lu)
  stat = fclos(lu)

end program tests
