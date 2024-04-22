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
!================================================================================================
program tests
#include <rmn/msg.h>
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_free, vgd_new,vgd_levels,operator(==), VGD_OK
  use Unit_Testing, only: ut_report


  implicit none

  integer :: lu=10,fnom,fstouv,fstinf,fstprm,fclos,ier,key,kind
  real :: pres
  logical :: ok = .false.
  type(vgrid_descriptor) :: vgd, vgd2
  character(len=1) :: blk_S
  character(len=72) :: file_S="data/dm_5005_from_model_run"

  ! For fstprm
  integer :: ig1,ig2,ig3,ig4,dateo,deet,npas,datyp,nbits
  integer :: ni,nj,nk
  integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=2) :: typvar
  character(len=4) :: nomvar
  character(len=12) :: etiket  

  call msg_verbosity(MSG_DEBUG)

  ier=fnom(lu,file_S,"RND+R/O",0)
  if(ier.lt.0)then
     print*,'(Test) ERROR with fnom on file ',trim(file_S)
     error stop 1
  endif
  ier=fstouv(lu,'RND')
  if(ier.lt.0)then
     print*,'(Test) No record in RPN file ',trim(file_S)
     error stop 1
  endif
  key = fstinf(lu,ni,nj,nk,-1," ",-1,-1,-1,' ',"TT")
  if(key < 0 )then
     print*,'Error cannot find a record TT in file ',trim(file_S)
     error stop 1
  endif
  ier=fstprm(key,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,&
       ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
       ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
  call convip(ip1,pres,kind,-1,blk_S,.false.)
  
  key = fstinf(lu,ni,nj,nk,-1," ",-1,-1,-1,' ',"!!")
  if(key < 0 )then
     print*,'Error cannot find !! in file ',trim(file_S)
     error stop 1
  endif
  ier=fstprm(key,dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,&
       ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
       ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)

  if(VGD_OK /= vgd_new(vgd,unit=lu,format="fst",ip1=ip1,ip2=ip2,kind=kind))then
     print*,"Error could not read vgrid descriptor 1"
     error stop 1
  endif
  if(VGD_OK /= vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket=etiket,ip1=ip1,ip2=ip2,ip3=0,kind=kind))then
     print*,"Error could not read vgrid descriptor 2"
     error stop 1
  endif
  if (.not. vgd == vgd2) then
     print*,'Error vgd and vgd2 must be equal but are not'
     error stop 1
  endif

  ! Test wrong datev
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=354514400,etiket=etiket,ip1=ip1,ip2=ip2,ip3=0,kind=kind,quiet=.true.))then
     print*,'Error should have return an error and did not with wrong datev'     
     error stop 1
  endif
  ! Test wrong etiket
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket="WRONG ETIKET",ip1=ip1,ip2=ip2,ip3=0,kind=kind,quiet=.true.) )then
     print*,'Error should have return an error and did not with wrong etiket'     
     error stop 1
  endif
  ! Test wrong ip1
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket=etiket,ip1=ip1+1,ip2=ip2,ip3=0,kind=kind,quiet=.true.) )then
     print*,'Error should have return an error and did not with wrong ip1'     
     error stop 1
  endif
  ! Test wrong ip2
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket=etiket,ip1=ip1,ip2=ip2+1,ip3=0,kind=kind,quiet=.true.) )then
     print*,'Error should have return an error and did not with wrong ip2'     
     error stop 1
  endif
  ! Test wrong ip3
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket=etiket,ip1=ip1,ip2=ip2,ip3=1,kind=kind,quiet=.true.) )then
     print*,'Error should have return an error and did not with wrong ip3'     
     error stop 1
  endif
  ! Test wrong kind
  print*,'The following error message on legacy encoding is normal for test wrong kind'
  if ( VGD_OK == vgd_new(vgd2,unit=lu,format="fst",datev=0,etiket=etiket,ip1=ip1,ip2=ip2,ip3=0,kind=1,quiet=.true.) )then
     print*,'Error should have return an error and did not with wrong kind'     
     error stop 1
  endif
  
  ier=fclos(lu)

  ok = .true.
  
  call ut_report(ok,message='Grid_Descriptors::vgd_levels level calculation status')
     
end program tests
