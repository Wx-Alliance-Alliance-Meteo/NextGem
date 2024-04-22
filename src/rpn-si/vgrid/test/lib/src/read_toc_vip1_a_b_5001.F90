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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_print,vgd_get,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: stat,i
  integer :: fnom,fstouv,fstfrm
  integer, dimension(:), pointer :: vip1
  real(kind=8), dimension(:), pointer :: cam,cbm
  real, dimension(:), pointer :: vcrd

  nullify(vip1,cam,cbm,vcrd)

  stat=fnom(lu,"data/dm_5001_from_model_run_plus_toc","RND",0)
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
  stat = vgd_new(d,unit=lu,format="fst",ip1=2009,ip2=1000)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_new'
     stat=fstfrm(lu)
     error stop 1
  endif
  stat = vgd_print(d)

  stat = vgd_get(d,key='VIPM - level ip1 list',value=vip1)
  print*,'vip1'
  print*,vip1
  stat = vgd_get(d,key='CA_M - vertical A coefficient (m)',value=cam)
  print*,'cam'
  print*,cam
  stat = vgd_get(d,key='CB_M - vertical B coefficient (m)',value=cbm)
  print*,'cbm'
  print*,cbm

  ! Print a and b for ip1=93423264
  do i=1,size(vip1)
     if(vip1(i).eq.93423264)then
        print*,'ip1,a,b=',vip1(i),cam(i),cbm(i)
        exit
     endif
  enddo

  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_print 5001')

  stat=fstfrm(lu)

end program constructor
