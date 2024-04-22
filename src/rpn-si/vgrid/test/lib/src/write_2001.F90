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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_write,vgd_print,vgd_free,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d,d_rtn
  integer :: stat,lu=0,fnom,fstouv,fstfrm,fclos,k,ip1,ip2
  integer, parameter :: LEVS=4
  integer, dimension(LEVS) :: ip1s
  real :: pres
  real(kind=8), dimension(LEVS) :: pres_8=(/1000.,925.,850.,700./),bb=(/0.,0.,0.,0./)
  real(kind=8), dimension(:), pointer :: pres_rtn
  character(len=1) :: dummy_S=" "

  nullify(pres_rtn)

  call system('rm -f test.fst')
  stat=fnom(lu,"test.fst","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  do k=1,LEVS
     pres = pres_8(k)
     call convip(ip1s(k),pres,2,2,dummy_S,.false.)
     print*,'ip1s(k),pres',ip1s(k),pres
  enddo

  pres_8 = pres_8*100. !convert mb to Pa
  stat = vgd_new(d,kind=2,version=1,nk=size(pres_8),ip1_m=ip1s,a_m_8=pres_8,b_m_8=bb,ip1=111,ip2=222)

  stat = vgd_print(d)
  if(stat==VGD_ERROR)then
     call system('rm -f test.fst')
     error stop 1
  endif
  stat = vgd_write(d,unit=lu,format='fst')
  stat = vgd_get(d,'IP_1 - record ip1',ip1)
  stat = vgd_get(d,'IP_2 - record ip2',ip2)
  stat = vgd_new(d_rtn,unit=lu,format='fst',ip1=ip1,ip2=ip2)
  if(stat==VGD_ERROR)then
     call system('rm -f test.fst')
     error stop 1
  endif

  stat = vgd_get(d_rtn,key='COFA - vertical A coefficient',value=pres_rtn)

  stat = vgd_free(d)
  stat = vgd_free(d_rtn)

  call ut_report(abs(pres_rtn(1)-pres_8(1)) < epsilon(pres_8),'Grid_Descriptors::vgd_new vertical build initializer (2001) value')

  stat=fstfrm(lu)
  stat=fclos(lu)

  ! Cleanup
  call system('rm -f test.fst')

end program constructor
