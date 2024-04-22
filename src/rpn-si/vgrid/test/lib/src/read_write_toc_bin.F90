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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_write,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d,d_rtn
  integer :: stat,lu=0,k
  integer, parameter :: LEVS=4
  integer, dimension(LEVS) :: ip1s
  real(kind=8), dimension(LEVS) :: pres=(/1000.,925.,850.,700./),b=(/0.,0.,0.,0./)
  real(kind=8), dimension(:), pointer :: pres_rtn

  nullify(pres_rtn)

  lu = 10
  open(unit=lu,file='test.bin',form='unformatted')

  ! Construct a new set of 3D coordinate descriptors
  do k=1,LEVS
     call convip(ip1s(k),pres(k),2,2,'',.false.)
  enddo
  pres = pres*100. !convert mb to Pa
  stat = vgd_new(d,kind=2,version=1,nk=size(pres),ip1_m=ip1s,a_m_8=pres,b_m_8=b)
  stat = vgd_write(d,unit=lu,format='bin')
  rewind(lu)
  stat = vgd_new(d_rtn,unit=lu,format='bin',ip1=0,ip2=0)
  stat = vgd_get(d_rtn,key='COFA - vertical A coefficient',value=pres_rtn)
  call ut_report(abs(pres_rtn(1)-pres(1)) < epsilon(pres),'Grid_Descriptors::vgd_new vertical build initializer (2001) value')

  close(10)

  ! Cleanup
  call system('rm -f test.fst')

end program constructor
