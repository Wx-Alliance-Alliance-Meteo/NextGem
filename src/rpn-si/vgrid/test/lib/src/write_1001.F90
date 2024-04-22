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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,vgd_print,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lui=10,luo=20,nmax=1000
  integer, dimension(nmax) :: liste
  integer :: stat,k,infon,kind,nii,njj,nkk,ig1tt,ig2tt
  integer :: fnom,fstouv,fstinl,fstprm
  real, dimension(:),allocatable :: sig
  real(kind=8), dimension(:),allocatable :: a_m_8,b_m_8
  integer, dimension(:),allocatable :: ip1_m
  character(len=10) :: blk_S

  ! For fstprm
  integer :: ig3,ig4,dateo,deet,npas,datyp,nbits
  integer :: ni,nj,nk
  integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=1) :: grtyp
  character(len=2) :: typvar
  character(len=4) :: nomvar
  character(len=12) :: etiket

  stat=fnom(lui,"data/dm_1001_from_model_run","RND+OLD",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom 1'
     error stop 1
  endif
  stat=fstouv(lui,'RND')

  stat=fnom(luo,"data_out/dm_1001","RND",0)
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

  allocate(ip1_m(nk),a_m_8(nk),b_m_8(nk),sig(nk))

  do k=1,nk
     stat=fstprm(liste(k),dateo,deet,npas,nii,njj,nkk,nbits,datyp,ip1,&
          ip2,ip3,typvar,nomvar,etiket,grtyp,ig1tt,ig2tt,ig3, &
          ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
     if(stat.lt.0)then
        print*,'ERROR with fstprm'
        call close_fst(lui,luo)
        error stop 1
     endif
     call convip(ip1,sig(k),kind,-1,blk_S,.false.)
     ip1_m(k)=ip1
     print*,sig(k)
     if(kind.ne.1)then
        print*,'ERROR wrong kind:',kind
        call close_fst(lui,luo)
        error stop 1
     endif
  enddo



  do k=1,nk
     b_m_8(k)= sig(k)
     a_m_8(k)= 0.0
     print*,'k,b_m_8(k),a_m_8(k)=',k,b_m_8(k),a_m_8(k)
     flush(6)
  enddo

  stat=vgd_new(&
       d,                  &
       kind    = 1,        &
       version = 1,        &
       nk      = nk,       &
       ip1    = ig1tt,     &
       ip2    = ig2tt,     &
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
  stat = vgd_print(d)
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
   integer lui,luo,stat
   integer :: fstfrm

   stat=fstfrm(lui)
   stat=fstfrm(luo)

   print*,'fst files closed'

end subroutine close_fst
