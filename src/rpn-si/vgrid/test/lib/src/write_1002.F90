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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lui=10,luo=20,nmax=1000
  integer, dimension(nmax) :: liste
  integer :: stat,k,infon,kind,keyhy,keypt,nii,njj,nkk,ig1tt,ig2tt
  integer :: fnom,fstouv,fstinl,fstprm,fstinf,fstluk
  real, dimension(:),allocatable :: eta
  real, dimension(:,:,:),allocatable :: ptop
  real(kind=8), dimension(:),allocatable :: a_m_8,b_m_8
  real(kind=8) :: ptop_8
  real :: pp,rcoef1
  integer, dimension(:),allocatable :: ip1_m
  character(len=10) :: blk_S

  ! For fstprm
  integer :: ig1,ig2,ig3,ig4,dateo,deet,npas,datyp,nbits
  integer :: ni,nj,nk
  integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=2) :: typvar
  character(len=4) :: nomvar
  character(len=12) :: etiket

  stat=fnom(lui,"data/dm_1002_from_model_run","RND+OLD",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom 1'
     error stop 1
  endif
  stat=fstouv(lui,'RND')

  stat=fnom(luo,"data_out/dm_1002","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom 2'
     error stop 1
  endif
  stat=fstouv(luo,'RND')

  stat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',liste,infon,nmax)
  if(stat.lt.0)then
     print*,'ERROR with fstinl'
     call close_fst(lui,luo)
     error stop 1
  endif
  
  nk=infon

  if(nk.le.0)then
     print*,'ERROR no TT records'
     call close_fst(lui,luo)
     error stop 1
  endif

  allocate(ip1_m(nk),a_m_8(nk),b_m_8(nk),eta(nk))

  do k=1,nk
     stat=fstprm(liste(k),dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1tt,ig2tt,ig3, &
          ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
     if(stat.lt.0)then
        print*,'ERROR with fstprm'
        call close_fst(lui,luo)
        error stop 1
     endif
     call convip(ip1,eta(k),kind,-1,blk_S,.false.)
     ip1_m(k)=ip1
     print*,eta(k)
     if(kind.ne.1)then
        print*,'ERROR wrong kind:',kind
        call close_fst(lui,luo)
        error stop 1
     endif
  enddo

  keypt=fstinf(lui,nii,njj,nkk,-1,' ',-1,-1,-1,' ','PT')
  if(keypt.lt.0)then
     print*,'WARNING cannot find PT. Will look for HY'
     keyhy=fstinf(lui,nii,njj,nkk,-1,' ',-1,-1,-1,' ','HY')
     if(keyhy.lt.0)then
        print*,'ERROR with fstinf on HY. Will exit(1)'
        call close_fst(lui,luo)
        error stop 1
     endif
  endif

  if (keypt.ge.0) then
     stat=fstprm(keypt,dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
          ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
     allocate(ptop(nii,njj,nkk))
     stat=fstluk(ptop,keypt,nii,njj,nkk)
     ptop_8=ptop(1,1,1)*100.d0
  endif
  if (keypt.lt.0 .and. keyhy.ge.0) then
     stat=fstprm(keyhy,dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
          ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
     rcoef1=ig2/1000.
     if (ig2.ne.1000) then
        print*,'ERROR with rcoef not equal 1.0 from HY. rcoef= ',rcoef1
        call close_fst(lui,luo)
        error stop 1
     endif
     call convip(ip1,pp,kind,-1,blk_S,.false.)
     ptop_8=pp*100.d0
  endif
  print*,'nk,ptop_8',nk,ptop_8

  do k=1,nk
     b_m_8(k)= eta(k)
     a_m_8(k)= ptop_8*(1.0-eta(k))
     print*,'k,b_m_8(k),a_m_8(k)=',k,b_m_8(k),a_m_8(k)
     flush(6)
  enddo

  stat=vgd_new(&
       d,                  &
       kind    = 1,        &
       version = 2,        &
       nk      = nk,       &
       ip1    = ig1tt,     &
       ip2    = ig2tt,     &
       ptop_8  = ptop_8,   &
       a_m_8   = a_m_8,    &
       b_m_8   = b_m_8,    &
       ip1_m   = ip1_m)

  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_new'
     call close_fst(lui,luo)
     error stop 1
  endif

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_write(d,luo,format='fst')
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_write'
     call close_fst(lui,luo)
     error stop 1
  endif

  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_write')

  call close_fst(lui,luo)

end program constructor

subroutine close_fst(lui,luo)

   implicit none
   integer lui,luo,fstfrm,stat

   stat=fstfrm(lui)
   stat=fstfrm(luo)

   print*,'fst files closed'

end subroutine close_fst
