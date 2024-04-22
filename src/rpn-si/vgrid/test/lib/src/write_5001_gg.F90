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
  integer :: stat,k,infon,kind,keyhy,nii,njj,nkk,ig1tt,ig2tt
  integer :: fnom,fstouv,fstinl,fstprm,fstinf
  real, dimension(:),allocatable :: hyb
  real(kind=8), dimension(:),allocatable :: a_m_8,b_m_8
  real(kind=8) :: pref_8,ptop_8,hybtop_8
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

  stat=fnom(lui,"data/gg_5001_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom 1'
     error stop 1
  endif
  stat=fstouv(lui,'RND')

  stat=fnom(luo,"data_out/gg_5001","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom 2'
     error stop 1
  endif
  stat=fstouv(luo,'RND')

  stat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',liste,infon,nmax)
  if(stat.lt.0)then
     print*,'ERROR with fstinl'
     error stop 1
  endif
  
  nk=infon

  if(nk.le.0)then
     print*,'ERROR no TT records'
     call close_fst(lui,luo)
     error stop 1
  endif

  allocate(ip1_m(nk),a_m_8(nk),b_m_8(nk),hyb(nk))

  do k=1,nk
     stat=fstprm(liste(k),dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1tt,ig2tt,ig3, &
          ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
     if(stat.lt.0)then
        print*,'ERROR with fstprm'
        call close_fst(lui,luo)
        error stop 1
     endif
     call convip(ip1,hyb(k),kind,-1,blk_S,.false.)
     ip1_m(k)=ip1
     print*,hyb(k)
     if(kind.ne.5)then
        print*,'ERROR wrong kind:',kind
        call close_fst(lui,luo)
        error stop 1
     endif
  enddo

  keyhy=fstinf(lui,nii,njj,nkk,-1,' ',-1,-1,-1,' ','HY')
  if(keyhy.lt.0)then
     print*,'ERROR with fstinf on HY'
     call close_fst(lui,luo)
     error stop 1
  endif
  stat=fstprm(keyhy,dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
       ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
       ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
  rcoef1=ig2/1000.
  pref_8=ig1*100.d0
  call convip(ip1,pp,kind,-1,blk_S,.false.)
  ptop_8=pp*100.d0

  print*,'nk,ptop_8,pref_8,rcoef1=',nk,ptop_8,pref_8,rcoef1

  hybtop_8=ptop_8/pref_8

  do k=1,nk
     b_m_8(k)=((hyb(k)-hybtop_8)/(1.-hybtop_8))**rcoef1
     a_m_8(k)=(hyb(k)-b_m_8(k))*pref_8
     print*,'k,b_m_8(k),a_m_8(k)=',k,b_m_8(k),a_m_8(k)
     flush(6)
  enddo

  stat=vgd_new(&
       d,                  &
       kind    = 5,        &
       version = 1,        &
       nk      = nk,       &
       ip1    = 0 ,        &
       ip2    = 0 ,        &
       ptop_8  = ptop_8,   &
       pref_8  = pref_8,   &
       rcoef1  = rcoef1,   &
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
