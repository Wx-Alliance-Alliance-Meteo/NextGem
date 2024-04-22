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
module mod_constructor_table
  
  type FSTD_ext
     integer :: ig1,ig2,ig3,ig4,dateo,deet,npas,datyp,nbits,ni,nj,nk
     integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3,datev
     character(len=1) :: grtyp
     character(len=2) :: typvar
     character(len=4) :: nomvar
     character(len=12) :: etiket
  end type FSTD_ext
contains
   integer function my_fstprm(fstkey,record) result(status)
      ! Use fstprm function to get information about the record
      integer, intent(in) :: fstkey               !Key from FST file record
      type(FSTD_ext) :: record                    !Record information
      integer :: error
      integer, external :: fstprm,fstinf
      real(kind=8) :: nhours
      status = -1
      error=fstprm(fstkey,record%dateo,record%deet,record%npas, &
           record%ni,record%nj,record%nk,record%nbits,record%datyp,record%ip1,record%ip2, &
           record%ip3,record%typvar,record%nomvar,record%etiket,record%grtyp, &
           record%ig1,record%ig2,record%ig3,record%ig4,record%swa, &
           record%lng,record%dltf,record%ubc,record%extra1,record%extra2, &
           record%extra3)
      if (error < 0) then
         write(6,*) 'cannot fstprm for fstkey ',fstkey
         return
      end if
      nhours=record%deet*record%npas/3600.d0
      call incdatr(record%datev,record%dateo,nhours)
      status = 0
   end function my_fstprm
end module mod_constructor_table

program constructor
  use mod_constructor_table, only: FSTD_ext,my_fstprm
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_new,Vgd_get,Vgd_print,VGD_OK,VGD_ERROR
  use Unit_Testing, only: ut_report
  !
  implicit none
  !
  type(vgrid_descriptor) :: d
  type(FSTD_ext) :: prm
  integer, parameter :: nliste=2
  integer,dimension(nliste) :: lu=(/10,11/), &
       ip1,ip2,vcode=(/5001,5003/)
  integer :: stat,i,vers,key,ni,nj,nk
  integer :: fnom,fstouv,fstinf,fstfrm,fclos,kind,fstlnk
  logical :: ok=.true.
  character(len=1000),dimension(nliste) :: files=(/ &
       "data/dm_5001_from_model_run_plus_toc    ", &
       "data/dm_5003_from_model_run             "&
       /)
  !
  do i=1,nliste
     stat=fnom(lu(i),files(i),"RND",0)
     if(stat.lt.0)then
        print*,'ERROR with fnom in file',trim(files(i))
        error stop 1
     endif
     stat=fstouv(lu(i),'RND')
     if(stat.le.0)then
        print*,'No record in RPN file',trim(files(i))
        error stop 1
     endif
     key=fstinf(lu(i),ni,nj,nk,-1,' ',-1,-1,-1,' ','TT')
     if(key.lt.0)then
        print*,'No TT record in RPN file',files(i)
        error stop 1
     endif
     stat=my_fstprm(key,prm)
     if(i.eq.3)then
        ! ig4 -> ip1 link
        ip1(i)=prm%ig4
        ip2(i)=-1
     else
        ! (ig1,ig2) -> (ip1,ip2) link
        ip1(i)=prm%ig1
        ip2(i)=prm%ig2
     endif
  enddo
  !
  stat = fstlnk(lu,nliste)
  !
  ! Test that an error is return if no !! are found with ips
  !stat = vgd_new(d,unit=lu(1),format="fst",ip1=1234,ip2=5678)
  !print*,'The above error is normal since there is not !! with these ips'
  !if(stat/=VGD_ERROR)ok=.false.
  !
  ! Test that an error is return if wild card are passed and there are more than on !! 
  !stat = vgd_new(d,unit=lu(1),format="fst",ip1=-1,ip2=-1)
  !print*,'The above error is normal with ips=-1 since there are 2 differents !! in linked files'
  !if(stat/=VGD_ERROR)ok=.false.
  !
  ! Test that an error is return if no ips is passed and there are more than on !!
  !stat = vgd_new(d,unit=lu(1),format="fst")
  !print*,'The above error is normal with no ips passed since there are 2 differents !! in linked files'
  !if(stat/=VGD_ERROR)ok=.false.
  !
  do i=1,nliste
     ! Construct a new set of 3D coordinate descriptors
     print*,'i, ip1(i),ip2(i)',i, ip1(i),ip2(i)
     stat = vgd_new(d,unit=lu(1),format="fst",ip1=ip1(i),ip2=ip2(i))
     if(stat.ne.VGD_OK)then
        print*,'ERROR: problem with vgd_new'
        error stop 1
     endif
     ! Get vcode to check if whats is read is ok
     stat = vgd_get(d,'KIND - vertical coordinate ip1 kind',kind)
     stat = vgd_get(d,'VERS - vertical coordinate version',vers)
     if(kind*1000+vers/=vcode(i))ok=.false.
  enddo
  !
  do i=1,nliste
     stat=fstfrm(lu(i))  
     stat=fclos(lu(i))
  enddo
  !
  call ut_report(ok,'Grid_Descriptors, vgd_get, LOGP 5002')  
  !
end program constructor
