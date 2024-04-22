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
program constructor_build_5999
   
   ! Note
   ! This program was used to create the dm_5999_from_model_run file from data/dm_5002_from_model_run

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,vgd_write,vgd_put,VGD_ERROR,VGD_OK
  use Unit_Testing, only: ut_report
  

  !
  implicit none
  !
  type(vgrid_descriptor) :: vgd,vgd2
  integer :: stat,lui=10,luo=11,fnom,fstouv,ier,fstfrm,i,j,fstinl,fstluk,fstprm,fstecr
  integer, dimension(:), pointer :: ip1_m
  real(kind=8), dimension(:), pointer :: a_m_8,b_m_8
  real, dimension(:,:), pointer :: ff
  real :: dummy
  logical :: OK=.true.
  character(len=128) :: file
  character(len=4), dimension(6) :: nomvars=(/'>>','^^','P0','TT','UU','PX'/)
  ! Variable for fstprm
  integer ::dateo, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
       ig2, ig3, ig4, ip1, ip2, ip3, lng, nbits,&
       ni,  nj, nk, npas, swa, ubc
  character(len=12) :: etiket
  character(len=4)  :: nomvar
  character(len=2)  :: typvar
  character(len=1)  :: grtyp
  integer, parameter :: nmax=1000
  integer, dimension(nmax) :: liste
  integer :: infon
 
  nullify(ff,ip1_m,a_m_8,b_m_8)

  file='data/dm_5001_from_model_run'
  ier=fnom(lui,file,"RND",0)
  if(ier.lt.0)then
     print*,'ERROR with fnom on file ',file
     error stop 1
  endif
  ier=fstouv(lui,'RND')
  if(ier.lt.0)then
     print*,'No record in RPN file ',file
     error stop 1
  endif
  if(vgd_new(vgd,unit=lui,format="fst") == VGD_ERROR)OK=.false.  

  ! Retrieve A, B and ip1
  if(vgd_get(vgd,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8) == VGD_ERROR)OK=.false.
  if(vgd_get(vgd,key='CB_M - vertical A coefficient (m)'   ,value=b_m_8) == VGD_ERROR)OK=.false.
  if(vgd_get(vgd,key='VIPM - ip1 momentun'                 ,value=ip1_m) == VGD_ERROR)OK=.false.

  ! Try create with ip1 not of kind 5
  print*,'TEST build with wrong ip1 kind'
  print*,'   The following error message on ip1 kind is expected'
  do i=1,size(ip1_m)
     ip1_m(i)=1000+i
  enddo
  stat=vgd_new(vgd2,kind=5,version=999,nk=size(a_m_8),ip1=1,ip2=2,&
       a_m_8=a_m_8,b_m_8=b_m_8,ip1_m=ip1_m)
  if(stat==VGD_OK)then
     print*,'Error with test on wrong ip1 kind, vgd_new should have error stop 1 error but did not'
     OK=.false.
  endif
  ! Try create with ip1 repetition
  print*,'TEST build with wrong ip1 list'
  print*,'   The following error message on wronf ip1 list is expected'
  ip1_m=93423264
  stat=vgd_new(vgd2,kind=5,version=999,nk=size(a_m_8),ip1=1,ip2=2,&
       a_m_8=a_m_8,b_m_8=b_m_8,ip1_m=ip1_m)
  if(stat==VGD_OK)then
     print*,'Error with test on wrong ip1 list, vgd_new should have return error but did not'
     OK=.false.
  endif   
  ! Build 5999  
  print*,'TEST build with correct parameter'
  if(vgd_get(vgd,key='VIPM - ip1 momentun'                 ,value=ip1_m) == VGD_ERROR)OK=.false.
  stat=vgd_new(vgd2,kind=5,version=999,nk=size(a_m_8),ip1=1,ip2=2,&
       a_m_8=a_m_8,b_m_8=b_m_8,ip1_m=ip1_m)
  if(stat==VGD_ERROR)OK=.false.
  !
  ! Run this section to create dm_5999_from_model_run file
  if(.false.)then
     file='dm_5999_from_model_run'
     ier=fnom(luo,file,"RND",0)
     if(ier.lt.0)then
        print*,'ERROR with fnom on file ',file
        error stop 1
     endif
     ier=vgd_print(vgd2)     
     ier=fstouv(luo,'RND')
     if(ier.lt.0)then
        print*,'No record in RPN file ',file
        error stop 1
     endif     
     ! Copy fields
     LOOP_ON_NOMVAR: do i=1,size(nomvars)
        print*,nomvars(i)
        ier = fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ',nomvars(i),liste,infon,nmax)
        if(infon.gt.0)then
           LOOP_ON_LEVELS: do j=1,infon
              if(associated(ff))then
                 if(size(ff,dim=1).ne.ni .or. size(ff,dim=2).ne.nj)then
                    deallocate(ff)
                    allocate(ff(ni,nj))
                 endif
              else
                 allocate(ff(ni,nj))
              endif
              print*,'liste(j)',liste(j)
              ier = fstluk(ff,liste(j),ni,nj,nk)
              ier = fstprm(liste(j),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
                   ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
                   dltf,ubc,extra1,extra2,extra3)
              ier=fstecr(ff,dummy,-nbits,luo,dateo,deet,npas, &
                   ni,nj,nk,ip1,ip2,ip3,typvar,nomvar,etiket,grtyp, &
                   ig1,ig2,ig3,ig4,datyp,.true.)
           enddo LOOP_ON_LEVELS
        endif
     enddo LOOP_ON_NOMVAR
     if(vgd_put(vgd2,key="IP_1",value=ig1) == VGD_ERROR)OK=.false.
     if(vgd_put(vgd2,key="IP_2",value=ig2) == VGD_ERROR)OK=.false.  
     if(vgd_write(vgd2,luo) == VGD_ERROR)OK=.false.
     ier=fstfrm(luo)
     
     error stop 1
     
  endif

  call ut_report(OK,'Grid_Descriptors:: vgd_new new_gen')
  !
end program constructor_build_5999
