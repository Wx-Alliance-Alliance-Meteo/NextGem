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
!================================================================================================
program tests
#include <rmn/msg.h>
  use vGrid_Descriptors, only: vgd_levels,vgd_putopt,VGD_ERROR
  use Unit_Testing, only: ut_report


  implicit none

  integer :: stat,lu=10,do_it,fstfrm,fclos,i,ier
  logical :: ok=.true.
  integer, parameter :: nfiles=9
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run",&
       "data/dm_1002_from_model_run",&
       "data/dm_5001_from_model_run",&
       "data/dm_5002_from_model_run",&
       "data/dm_5003_from_model_run",&
       "data/dm_5004_from_model_run",&
       "data/dm_5005_from_model_run",&
       "data/dm_5100_from_model_run",&
       "data/dm_5999_from_model_run"&
       /)
  
  call msg_verbosity(MSG_DEBUG)

  stat = vgd_putopt("ALLOW_SIGMA",.true.)

  do i=1,nfiles
     stat = do_it(lu,files(i))
     ier=fstfrm(lu)
     ier=fclos(lu)
     if(stat==VGD_ERROR)then
        ok=.false.
        exit
     endif
  enddo

  call ut_report(ok,message='Grid_Descriptors::vgd_levels level calculation status')

end program tests
!=========================================================================
!=========================================================================
!=========================================================================
integer function do_it(lu,file) result(status)

  use vGrid_Descriptors, only: vgd_levels,VGD_ERROR,VGD_OK
  use mod_constructor_table, only: FSTD_ext,my_fstprm 

  implicit none

  integer :: lu
  character(len=*) :: file

  !Local variables

  type(FSTD_ext) :: prm
  integer, parameter :: nmax=1000
  integer, dimension(nmax) :: liste
  integer :: ni,nj,nk,fnom,fstouv,fstkey,fstlir,fstinf,ier,kind,fstinl,infon
  real :: epsilon=0.01,pres
  real, dimension(:,:), pointer :: px
  real, dimension(:,:,:), pointer :: lev
  character(len=1) :: string

  nullify(lev,px)

  status=VGD_ERROR

  print*,'===================================================='
  print*,trim(file)

  ier=fnom(lu,file,"RND+R/O",0)
  if(ier.lt.0)then
     print*,'(Test) ERROR with fnom on file ',file
     return
  endif
  ier=fstouv(lu,'RND')
  if(ier.lt.0)then
     print*,'(Test) No record in RPN file ',file
     return
  endif

  ! Get physical levelling information
  ier = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ' ,'TT',liste,infon,nmax)
  if(infon.le.0)then
     print*,'(Test) ERROR cannot find any TT in input file ',file
     return
  endif    
  ! treat middle level record
  ier=my_fstprm(liste(infon/2),prm)
  call convip(prm%ip1,pres,kind,-1,string,.false.)
  print*,'(Test) Treating TT at ip1=',prm%ip1,' ->',pres

  fstkey = fstinf(lu,ni,nj,nk,-1,'',prm%ip1,-1,-1,'','TT')
  if(fstkey.le.0)then
     print*,'(Test) ERROR cannot find any TT in input file ',file
     return
  endif
  ier = vgd_levels(unit=lu,fstkeys=(/fstkey/),levels=lev)
  if(ier==VGD_ERROR)then
     print*,'(Test) Problem with vgd_levels'
     return
  endif

  allocate(px(ni,nj))  
  ier = fstlir(px,lu,ni,nj,nk,-1,'',prm%ip1,-1,-1,'','PX')
  if(ier.lt.0)then
     print*,'(Test) ERROR with fstlir on PX for ip1,file ',prm%ip1,file
     return
  endif

  if(abs(lev(10,10,1)/100.-px(10,10))>epsilon)then
     print*,'(Test) ERROR difference in pressure to high',lev(10,10,1)/100.,' VS',px(10,10)
     return
  else
     print*,'(Test) Difference in pressure are OK',lev(10,10,1)/100.,' VS',px(10,10)
  endif

  ! Test log option
  ier = vgd_levels(unit=lu,fstkeys=(/fstkey/),levels=lev,in_log=.true.)
  if(ier==VGD_ERROR)then
     print*,'(Test) Problem with vgd_levels in_log=.true.'
     return
  endif
  if(abs(lev(10,10,1)-log(px(10,10)*100.))>epsilon)then
     print*,'(Test) ERROR difference in log pressure to high ',lev(10,10,1),' VS',log(px(10,10)*100.)
     return
  else
     print*,'(Test) Difference in log pressure are OK ',lev(10,10,1),' VS',log(px(10,10)*100.)
  endif
  
  deallocate(px,lev)

  status=VGD_OK

end function do_it
