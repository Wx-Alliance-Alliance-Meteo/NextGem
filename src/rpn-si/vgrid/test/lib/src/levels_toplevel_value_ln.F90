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
  use vGrid_Descriptors, only: vgd_levels
  use Unit_Testing, only: ut_report

  implicit none

  integer :: stat,lu=0,ni,nj,nk,fnom,fstouv,fstfrm,fclos,fstkey,fstlir,fstinf
  real :: epsilon=0.01
  real, dimension(:,:), pointer :: px
  real, dimension(:,:,:), pointer :: lev

  nullify(px,lev)

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

  ! Get physical levelling information
  fstkey = fstinf(lu,ni,nj,nk,-1,'',93423264,-1,-1,'','TT')

  stat = vgd_levels(unit=lu,fstkeys=(/fstkey/),levels=lev,in_log=.true.)

  allocate(px(ni,nj))
  stat = fstlir(px,lu,ni,nj,nk,-1,'',93423264,-1,-1,'','PX')
  call ut_report(abs(lev(6,8,1)-log(px(6,8)*100.))<epsilon,message='Grid_Descriptors::vgd_levels level calculation status')

  stat=fstfrm(lu)
  stat=fclos(lu)

end program tests
