! * Libdescrip - Vertical grid descriptor library for FORTRAN programming
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

! NOTE : the control files for comparison can be produce by the test c_standard_atmophere_all.c

program stda76
#include <rmn/msg.h>
  use Unit_Testing, only: ut_report
  use vGrid_Descriptors, only: vgd_putopt, VGD_STDA76_SFC_T, VGD_STDA76_SFC_P, &
       VGD_ERROR
  implicit none
  integer :: stat,lu=10,i,stda_do_it
  logical :: ok=.true.
  integer, parameter :: nfiles=12
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_5999_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)
  call msg_verbosity(MSG_DEBUG)
  stat = vgd_putopt("ALLOW_SIGMA",.true.)

  if(abs(VGD_STDA76_SFC_T - 288.15)/288.15 > 1.e-5)then
     error stop 1
  endif
  if(abs(VGD_STDA76_SFC_P - 101325.)/101325. > 1.e-5)then
     error stop 1
  endif
 
  do i=1,nfiles
     stat = stda_do_it(lu,files(i))
     if(stat==VGD_ERROR)then
        error stop 1
     endif
  enddo
  call ut_report(ok,message='Grid_Descriptors::vgd_levels level calculation status')
end program stda76
!=========================================================================
!=========================================================================
!=========================================================================
integer function stda_do_it(lu,file) result(status)
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_stda76,vgd_standard_atmosphere_1976,VGD_OK,VGD_ERROR
  implicit none
  integer :: lu
  character(len=*) :: file
  !Local variables
  integer, dimension(:), pointer :: ip1s
  integer :: ier,fnom,fstouv,fstfrm,fclos,compare
  real, dimension(:), pointer :: temp, pres
  character(len=4) :: nomvar
  type(vgrid_descriptor) :: vgd
  status=VGD_ERROR
  nullify(ip1s,temp,pres)
  print*,'===================================================='
  print*,'Testing ',trim(file)
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
  ier = vgd_new(vgd,lu)
  if(ier==VGD_ERROR)then
     print*,'(Test) Problem with vgd_new'
     return
  endif
  if( vgd_get(vgd, "VIPT", ip1s) ==  VGD_ERROR )then
     print*,"ERROR with vgd_get for VIPT"
     return
  end if
  !==============================
  print*,"   Testing temperature"
  !------------------------------
  if( vgd_stda76(vgd, ip1s, temp, 'TEMPERATURE') == VGD_ERROR )then
     print*,"In test : ERROR with vgd_stda76_temp on TEMPERATURE"
     return
  endif
  if( compare(file, "_stda76_temp.txt", ip1s, temp, size(ip1s)) == VGD_ERROR )then
     print*,"In test : ERROR on values with vgd_stda76_temp on TEMPERATURE"
     return
  endif
  ! Testing compatibility wrapper
  if( vgd_standard_atmosphere_1976(vgd, ip1s, temp, 'TEMPERATURE') == VGD_ERROR )then
     print*,"In test : ERROR with vgd_standard_atmosphere_1976 on TEMPERATURE"
     return
  endif
  if( compare(file, "_stda76_temp.txt", ip1s, temp, size(ip1s)) == VGD_ERROR )then
     print*,'ERROR on values with vgd_standard_atmosphere_1976 on TEMPERATURE'
     return
  endif

  ier = vgd_get(vgd, "RFLD", nomvar)
  if( trim(nomvar) == "ME" ) then
     !===========================
     print*,"   Testing pressure"
     !---------------------------
     if( vgd_stda76(vgd, ip1s, pres, 'PRESSURE') == VGD_ERROR )then
        print*,"In test : ERROR with vgd_stda76 on PRESSURE with hgts"
        return
     endif
     if( compare(file, "_stda76_pres.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR on values with vgd_stda76 on PRESSURE with hgts"
        return
     endif
     ! Testing compatibility wrapper
     if( vgd_standard_atmosphere_1976(vgd, ip1s, pres, 'PRESSURE') == VGD_ERROR )then
        print*,"In test : ERROR with vgd_standard_atmosphere_1976 on PRESSURE with hgts"
        return
     endif
     if( compare(file, "_stda76_pres.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR on values with vgd_standard_atmosphere_1976 on PRESSURE with hgts"
        return
     endif
     
     !============================================
     print*,"   Testing pressure, option sfc_pres"
     !--------------------------------------------
     if( vgd_stda76(vgd, ip1s, pres, 'PRESSURE', sfc_pres=100000.) == VGD_ERROR )then
        print*,"In test : ERROR with vgd_stda76 on PRESSURE with hgts with sfc_pres"
        return
     endif
     if( compare(file, "_stda76_pres_sfc_pres_100000.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR values with vgd_stda76 on PRESSURE with hgts with sfc_pres"
        return
     endif
     ! Testing compatibility wrapper
      if( vgd_standard_atmosphere_1976(vgd, ip1s, pres, 'PRESSURE', sfc_pres=100000.) == VGD_ERROR )then
        print*,"In test : ERROR with vgd_standard_atmosphere_1976 on PRESSURE with hgts with sfc_pres"
        return
     endif
     if( compare(file, "_stda76_pres_sfc_pres_100000.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR on values with vgd_standard_atmosphere_1976 on PRESSURE with hgts with sfc_pres"
        return
     endif
     !============================================
     print*,"   Testing pressure, option sfc_temp"
     !--------------------------------------------
     if( vgd_stda76(vgd, ip1s, pres, 'PRESSURE', sfc_temp=273.) == VGD_ERROR )then
        print*,"In test : ERROR with vgd_stda76 on PRESSURE with hgts with sfc_temp"
        return
     endif
     if( compare(file, "_stda76_pres_sfc_temp_273.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR on values with vgd_stda76 on PRESSURE with hgts with sfc_temp"
        return
     endif
     ! Testing compatibility wrapper
     if( vgd_standard_atmosphere_1976(vgd, ip1s, pres, 'PRESSURE', sfc_temp=273.) == VGD_ERROR )then
        print*,"In test : ERROR with vgd_standard_atmosphere_1976 on PRESSURE with hgts with sfc_temp"
        return
     endif
     if( compare(file, "_stda76_pres_sfc_temp_273.txt", ip1s, pres, size(ip1s)) == VGD_ERROR )then
        print*,"In test : ERROR on values with vgd_standard_atmosphere_1976 on PRESSURE with hgts with sfc_temp"
        return
     endif
     
     deallocate(pres)
  endif
  deallocate(ip1s)
  deallocate(temp)
  ier=fstfrm(lu)
  ier=fclos(lu)
  status=VGD_OK
end function stda_do_it

integer function compare(F_filename, F_filetype, F_ip1s, F_temp, F_nl) result(status)
  
  use vgrid_descriptors, only : VGD_OK, VGD_ERROR

  implicit none
  integer :: F_nl
  character(len=*) :: F_filename, F_filetype
  integer, dimension(F_nl) :: F_ip1s
  integer :: k, nl
  real, dimension(F_nl) :: F_temp
  ! Local variables
  integer :: ip1
  real :: temp

  status = VGD_ERROR

  open(unit=62, file=trim(F_filename)//trim(F_filetype), status='OLD')

print*,trim(F_filename)//trim(F_filetype)

  read(62,'(5x,i8)')nl
  if( nl /= F_nl) then
     print*,"ERROR in tests, size problem with validation file"
     close(62)
     return
  endif
  do k=1,nl
     !ip1   95101840 value 269.503143
     !1234567890123456789012345678901234567890
     !         1         2         3
     read(62,'(4x,i10,7x,f10.6)')ip1, temp
     if( ip1 /= F_ip1s(k) .or. abs(temp - F_temp(k) ) > .01 )then
        print*,"ERROR differences found, expecting:"
        print*,ip1, temp
        print*,"got:"
        print*,F_ip1s(k),F_temp(k)
        print*,"ERROR TEST failled for file ",F_filename
        close(62)
        return
     endif
  enddo
  print*,"TEST on file %s is OK",F_filename
  close(62)
status = VGD_OK
 
end function compare
