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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels
  use Unit_Testing, only: ut_report

  implicit none

  integer :: stat,lu=0,ni,nj,nk,fnom,fstouv,fstfrm,fclos,fstkey,fstlir,fstinf,lutxt=69,ip1,ip2
  real :: epsilon=0.01
  real, dimension(:,:), pointer :: p0,px
  real, dimension(:,:,:), pointer :: lev
  type(vgrid_descriptor) :: d

  nullify(p0,px,lev)

  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  open(unit=lutxt,file='data/dm_5002_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Get physical levelling information
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)
  fstkey = fstinf(lu,ni,nj,nk,-1,'',93423264,-1,-1,'','TT')
  allocate(p0(ni,nj))
  stat = fstlir(p0,lu,ni,nj,nk,-1,'',-1,-1,-1,'','P0')
  p0 = p0*100. !mb to Pa
  allocate(px(ni,nj))
  stat = fstlir(px,lu,ni,nj,nk,-1,'',93423264,-1,-1,'','PX')
  stat = vgd_levels(d,sfc_field=p0,ip1_list=(/93423264/),levels=lev,in_log=.true.)
  call ut_report(abs(lev(2,5,1)-log(px(2,5)*100.))<epsilon,message='Grid_Descriptors::vgd_levels level calculation status')

  stat=fstfrm(lu)
  stat=fclos(lu)

end program tests
