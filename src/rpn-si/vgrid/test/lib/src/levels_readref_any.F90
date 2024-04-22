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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels, VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  integer :: stat,lu=0,i,ni,nj,nk,fnom,fstouv,fstfrm,fclos,fstlir,fstinl,count,fstprm, &
       dateo,deet,npas,nbits,datyp,ip2,ip3,ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,ex1,ex2,ex3
  integer :: testLev=10
  integer, dimension(100) :: fstkeys
  integer, dimension(:), allocatable :: ip1s
  real :: epsilon=0.01
  real, dimension(:,:), allocatable :: px
  real, dimension(:,:,:), pointer :: lev
  character(len=12) :: typvar,nomvar,etiket,grtyp
  type(vgrid_descriptor) :: d

  nullify(lev)
  
  stat=fnom(lu,"data/dm_5002_from_model_run-interp","RND",0)
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
  stat = vgd_new(d,unit=lu,format="fst",ip1=-1,ip2=-1)
  stat = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',fstkeys,count,size(fstkeys))
  allocate(ip1s(count))
  do i=1,count
     stat = fstprm(fstkeys(i),dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1s(i), &
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng,dltf,  &
          ubc,ex1,ex2,ex3)
  enddo
  allocate(px(ni,nj))
  stat = fstlir(px,lu,ni,nj,nk,-1,'',ip1s(testLev),-1,-1,'','PX')
  stat = vgd_levels(d,unit=lu,fstkeys=fstkeys(1:count),levels=lev)
  if(stat == VGD_ERROR)then
     print*,'ERROR with vgd_levels'
     error stop 1
  end if
  call ut_report(abs(lev(10,5,testLev)/100.-px(10,5))<epsilon,message='Grid_Descriptors::vgd_levels level calculation status')

  stat=fstfrm(lu)
  stat=fclos(lu)

end program tests
