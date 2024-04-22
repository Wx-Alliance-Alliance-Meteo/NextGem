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
program constructor_build_all
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_free,Vgd_new,VGD_ERROR,vgd_putopt
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: vgd
  integer, parameter :: lu=10
  integer :: stat,i,check_build,check_get_ABC
  integer :: fnom,fstouv,fstfrm,fclos
  logical :: ok

  integer, parameter :: nfiles=15
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_2001_from_editfst             ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5003_from_model_run           ",&
       "data/dm_5004_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_5999_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)

  stat = vgd_putopt("ALLOW_SIGMA",.true.)

  ok=.true.

  do i=1,nfiles
     print*,'======================================================'
     print*,'file = ',trim(files(i))
     stat=fnom(lu+i,files(i),"RND",0)
     if(stat.lt.0)then
        print*,'ERROR with fnom on file ',trim(files(i))
        error stop 1
     endif
     stat=fstouv(lu+i,'RND')
     if(stat.le.0)then
        print*,'No record in RPN file ',trim(files(i))
        error stop 1
     endif     
     if( vgd_new(vgd,unit=lu+i,format="fst",ip1=-1,ip2=-1) == VGD_ERROR)then
        print*,'ERROR: problem with vgd_new on file ',trim(files(i))
        stat=fstfrm(lu+i)
        error stop 1
     endif
     stat=fstfrm(lu+i)  
     stat=fclos(lu+i)
     if( check_get_ABC(vgd) == VGD_ERROR) ok=.false.
     if( check_build(vgd) == VGD_ERROR) ok=.false.
     if( vgd_free(vgd) == VGD_ERROR)then
        print*,'Error with vgd_free(vgd) on file ',trim(files(i))
        error stop 1
     endif
  end do

  call ut_report(ok,'Grid_Descriptors, vgd_get, LOGP 5002')  

end program constructor_build_all

!=============================================================
!=============================================================
!=============================================================
!=============================================================
integer function check_get_ABC(F_vgd) result(istat)
  ! The goal is to check if all coef A, B and C can be retreived.
  use, intrinsic :: iso_fortran_env
  use Vgrid_Descriptors, only: Vgrid_descriptor, vgd_get, VGD_ERROR, VGD_OK,VGD_NO_REF_NOMVAR
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  ! Local varibales
  integer :: i,k, vcode, ier,nk
  integer, parameter :: nkeys=9
  character(len=4), dimension(9) :: keys_S = (/"CA_M","CB_M","CC_M","CA_T","CB_T","CC_T","CA_W","CB_W","CC_W"/)
  character(len=4) :: rfls_S, key_nk
  real(kind=REAL64), dimension(:), pointer :: value_1d_8  
  istat = VGD_ERROR
  nullify(value_1d_8)
  ! Get Vcode
  if( vgd_get(F_vgd,'vcode',vcode) == VGD_ERROR )error stop 1
  write(6,'("=============== Testing check_get_ABC for Vcode ",i10,"=======================")')vcode
  do i=1,nkeys
     print*,'Checking key ',keys_S(i)
     flush(6)
     if(vgd_get(F_vgd,keys_S(i),value_1d_8) == VGD_ERROR)then
        print*,'Problem in test constructor, s/r check_get_ABC, cannot get key ',keys_S(i)       
        error stop
        return
     endif
     if(.not. associated(value_1d_8))then
        print*,'Problem in test constructor, s/r check_get_ABC, value_1d_8 is not assiciated for ',keys_S(i)
        error stop
        return
     endif
     ! Test that C coefs are zero if non SLEVE
     ier = vgd_get(F_vgd, "RFLS", rfls_S,quiet=.true.)
     if( rfls_S == VGD_NO_REF_NOMVAR .and. keys_S(i)(2:2) == "C")then
        print*,'Testing ',keys_S(i),' for zero value'
        !NL_M, NL_T or NL_W
        key_nk="NL_"//keys_S(i)(4:4)
        ier = vgd_get(F_vgd,key_nk,nk)
        do k=1,nk
           if(value_1d_8(k) /= 0.d0)then
              print*,'Problem in test constructor, s/r check_get_ABC, some value_1d_8 are not zero and should for key ',keys_S(i)
              print*,value_1d_8(k)
              error stop
              return
           endif
        end do
     endif
     deallocate(value_1d_8)
  end do  
  istat = VGD_OK
  return
end function check_get_ABC
!=============================================================
integer function check_build(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor, vgd_get, VGD_ERROR, VGD_OK
  
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: check_build_1001_2001_5999, check_build_1002, check_build_4001, check_build_5001, check_build_5002, check_build_5005, &
       check_build_5100, check_build_21001, check_build_21002
  integer :: vcode

  istat = VGD_ERROR

  ! Get Vcode
  if( vgd_get(F_vgd,'vcode',vcode) == VGD_ERROR )error stop 1
  print*,'Testing build for Vcode ', vcode

  select case (vcode)
  case (1001,2001)
     if(check_build_1001_2001_5999(F_vgd) == VGD_ERROR) error stop 1
  case (1002)
     if(check_build_1002(F_vgd) == VGD_ERROR) error stop 1
  case (4001)
     if(check_build_4001(F_vgd) == VGD_ERROR) error stop 1
  case (5001)
     if(check_build_5001(F_vgd) == VGD_ERROR) error stop 1
  case (5002)
     if(check_build_5002(F_vgd) == VGD_ERROR) error stop 1
  case (5003)
     if(check_build_5002(F_vgd) == VGD_ERROR) error stop 1
  case (5004)
     if(check_build_5002(F_vgd) == VGD_ERROR) error stop 1
  case (5005)
     if(check_build_5005(F_vgd) == VGD_ERROR) error stop 1
  case (5100)
     if(check_build_5100(F_vgd) == VGD_ERROR) error stop 1
  case (5999)
     if(check_build_1001_2001_5999(F_vgd) == VGD_ERROR) error stop 1
  case (21001)
     if(check_build_21001(F_vgd) == VGD_ERROR) error stop 1
  case (21002)
     if(check_build_21002(F_vgd) == VGD_ERROR) error stop 1
  case DEFAULT
     print*,'ERROR, no function for checking Vcode ',vcode
     print*,'       please add this function in test'
     error stop 1
  end select

  istat = VGD_OK

  return
  
end function check_build

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function get_hyb_from_ip1(F_hyb,F_ip1,F_kind) result(istat)

  use Vgrid_Descriptors, only: VGD_OK

  implicit none
  
  integer, intent(in) :: F_kind
  integer, dimension(:), intent(in) :: F_ip1
  real, dimension(:),intent(inout) :: F_hyb
  
  ! Local varibales
  integer :: k
  
  do k=1,size(F_ip1)
     call CONVIP( F_ip1(k), F_hyb(k), F_kind, -1, "", .false.)
  enddo

  istat = VGD_OK

  return

end function get_hyb_from_ip1
!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_1001_2001_5999(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
 
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m), &
       ip1=-1,           &
       ip2=-1,           &
       a_m_8=a_m_8,         &
       b_m_8=b_m_8,         &
       ip1_m=ip1_m)         &
       == VGD_ERROR) return

  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif

  deallocate(ip1_m, a_m_8, b_m_8)

  istat = VGD_OK

  return

end function check_build_1001_2001_5999

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_1002(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
 
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m
  real(kind=8) :: ptop_8
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PTOP - top level pressure"        , ptop_8) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m), &
              ip1=-1,           &
              ip2=-1,           &
              ptop_8=ptop_8,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              ip1_m=ip1_m)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_1002

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_4001(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
 
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m), &
              ip1=-1,           &
              ip2=-1,           &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              ip1_m=ip1_m)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_4001

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_5001(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
  
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m
  real :: rcoef1
  real(kind=8) :: ptop_8, pref_8
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PTOP"                             , ptop_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PREF"                             , pref_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m), &
              ip1=-1,           &
              ip2=-1,           &
              ptop_8=ptop_8,     &
              pref_8=pref_8,       &
              rcoef1=rcoef1,       &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              ip1_m=ip1_m)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, a_m_8, b_m_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_5001

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_5002(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==),VGD_LEN_NAME, VGD_ERROR, VGD_OK
  
   
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version
  integer, dimension(:), pointer :: ip1_m, ip1_t
  real :: rcoef1, rcoef2
  real(kind=8) :: ptop_8, pref_8
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8, a_t_8, b_t_8
  real(kind=8), dimension(:,:,:), pointer :: table
  type(vgrid_descriptor) :: vgd
  character(len=VGD_LEN_NAME) :: rfld_S

  nullify(ip1_m, ip1_t, a_m_8, b_m_8, a_t_8, b_t_8, table)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_T - vertical A coefficient (t)", a_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_T - vertical B coefficient (t)", b_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPT - level ip1 list (t)"        , ip1_t) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PTOP"                             , ptop_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PREF"                             , pref_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_2"                             , rcoef2) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m)-1, &
              ip1=-1,            &
              ip2=-1,            &
              ptop_8=ptop_8,     &
              pref_8=pref_8,     &
              rcoef1=rcoef1,     &
              rcoef2=rcoef2,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              a_t_8=a_t_8,       &
              b_t_8=b_t_8,       &
              ip1_m=ip1_m,       &
              ip1_t=ip1_t)       &
              == VGD_ERROR) return

  if( vgd_get(F_vgd,"VTBL", table) == VGD_ERROR ) return
  print*,'table(2,3,1)',table(2,3,1)
  deallocate(table)
  if( vgd_get(vgd,"VTBL", table) == VGD_ERROR ) return
  print*,'table(2,3,1)',table(2,3,1)

  if( vgd_get(F_vgd,"RFLD",rfld_S ) == VGD_ERROR ) return
  print*,'rfld_S<',rfld_S,'>'
  if( vgd_get(vgd,"RFLD",rfld_S ) == VGD_ERROR ) return
  print*,'rfld_S<',rfld_S,'>'

  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     !ier = vgd_print(F_vgd)
     !ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, ip1_t, a_m_8, b_m_8, a_t_8, b_t_8, table)
  
  istat = VGD_OK
  
  return
  
end function check_build_5002
!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_5005(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
   
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m, ip1_t
  real :: rcoef1, rcoef2
  real(kind=8) :: pref_8
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8, a_t_8, b_t_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, ip1_t, a_m_8, b_m_8, a_t_8, b_t_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_T - vertical A coefficient (t)", a_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_T - vertical B coefficient (t)", b_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPT - level ip1 list (t)"        , ip1_t) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PREF"                             , pref_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_2"                             , rcoef2) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m)-2, &
              ip1=-1,            &
              ip2=-1,            &
              pref_8=pref_8,     &
              rcoef1=rcoef1,     &
              rcoef2=rcoef2,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              a_t_8=a_t_8,       &
              b_t_8=b_t_8,       &
              ip1_m=ip1_m,       &
              ip1_t=ip1_t)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, ip1_t, a_m_8, b_m_8, a_t_8, b_t_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_5005

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_5100(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
  
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local varibales
  integer :: kind, version
  integer, dimension(:), pointer :: ip1_m, ip1_t
  real :: rcoef1, rcoef2, rcoef3, rcoef4
  real(kind=8) :: pref_8
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, ip1_t, a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_M - vertical C coefficient (m)", c_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_T - vertical A coefficient (t)", a_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_T - vertical B coefficient (t)", b_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_T - vertical C coefficient (t)", c_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPT - level ip1 list (t)"        , ip1_t) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"PREF"                             , pref_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_2"                             , rcoef2) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_3"                             , rcoef3) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_4"                             , rcoef4) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m)-2, &
              ip1=-1,            &
              ip2=-1,            &
              pref_8=pref_8,     &
              rcoef1=rcoef1,     &
              rcoef2=rcoef2,     &
              rcoef3=rcoef3,     &
              rcoef4=rcoef4,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              c_m_8=c_m_8,       &
              a_t_8=a_t_8,       &
              b_t_8=b_t_8,       &
              c_t_8=c_t_8,       &
              ip1_m=ip1_m,       &
              ip1_t=ip1_t)       &
              == VGD_ERROR) return

  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     !ier = vgd_print(F_vgd)
     !ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, ip1_t, a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_5100
!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_21001(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
  
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local variables
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m, ip1_t
  real :: rcoef1, rcoef2, rcoef3, rcoef4
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, ip1_t, a_m_8, b_m_8, c_m_8, b_m_8, a_t_8, b_t_8, c_t_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_M - vertical C coefficient (m)", c_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_T - vertical A coefficient (t)", a_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_T - vertical B coefficient (t)", b_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_T - vertical C coefficient (t)", c_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPT - level ip1 list (t)"        , ip1_t) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_2"                             , rcoef2) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_3"                             , rcoef3) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_4"                             , rcoef4) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m)-2, &
              ip1=-1,            &
              ip2=-1,            &
              rcoef1=rcoef1,     &
              rcoef2=rcoef2,     &
              rcoef3=rcoef3,     &
              rcoef4=rcoef4,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              c_m_8=c_m_8,       &
              a_t_8=a_t_8,       &
              b_t_8=b_t_8,       &
              c_t_8=c_t_8,       &
              ip1_m=ip1_m,       &
              ip1_t=ip1_t)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, ip1_t, a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_21001

!=============================================================
!=============================================================
!=============================================================
!=============================================================

integer function check_build_21002(F_vgd) result(istat)

  use Vgrid_Descriptors, only: Vgrid_descriptor,vgd_new, vgd_get, vgd_print, operator(==), VGD_ERROR, VGD_OK
  
   
  implicit none
  type(vgrid_descriptor), intent(in) :: F_vgd
  
  ! Local variables
  integer :: kind, version, ier
  integer, dimension(:), pointer :: ip1_m, ip1_t, ip1_w
  real :: rcoef1, rcoef2, rcoef3, rcoef4
  real(kind=8), dimension(:), pointer :: a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8, a_w_8, b_w_8, c_w_8
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1_m, ip1_t, ip1_w, a_m_8, b_m_8, c_m_8, b_m_8, a_t_8, b_t_8, c_t_8, a_w_8, b_w_8, c_w_8)
  
  istat = VGD_ERROR
  
  if( vgd_get(F_vgd,"KIND", kind) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VERS", version) == VGD_ERROR ) return  
  if( vgd_get(F_vgd,"CA_M - vertical A coefficient (m)", a_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_M - vertical B coefficient (m)", b_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_M - vertical C coefficient (m)", c_m_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_T - vertical A coefficient (t)", a_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_T - vertical B coefficient (t)", b_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_T - vertical C coefficient (t)", c_t_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CA_W - vertical A coefficient (w)", a_w_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CB_W - vertical B coefficient (w)", b_w_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"CC_W - vertical C coefficient (w)", c_w_8) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPM - level ip1 list (m)"        , ip1_m) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPT - level ip1 list (t)"        , ip1_t) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"VIPW - level ip1 list (w)"        , ip1_w) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_1"                             , rcoef1) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_2"                             , rcoef2) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_3"                             , rcoef3) == VGD_ERROR ) return
  if( vgd_get(F_vgd,"RC_4"                             , rcoef4) == VGD_ERROR ) return

  if( vgd_new(vgd,kind,version,size(ip1_m)-2, &
              ip1=-1,            &
              ip2=-1,            &
              rcoef1=rcoef1,     &
              rcoef2=rcoef2,     &
              rcoef3=rcoef3,     &
              rcoef4=rcoef4,     &
              a_m_8=a_m_8,       &
              b_m_8=b_m_8,       &
              c_m_8=c_m_8,       &
              a_t_8=a_t_8,       &
              b_t_8=b_t_8,       &
              c_t_8=c_t_8,       &
              a_w_8=a_w_8,       &
              b_w_8=b_w_8,       &
              c_w_8=c_w_8,       &
              ip1_m=ip1_m,       &
              ip1_t=ip1_t,       &
              ip1_w=ip1_w)       &
              == VGD_ERROR) return
  
  if(.not. vgd == F_vgd)then
     print*,'ERROR : build descriptor not equal has the one read from the file see vgd_print output below'
     ier = vgd_print(F_vgd)
     ier = vgd_print(vgd)
     return
  endif
  
  deallocate(ip1_m, ip1_t, ip1_w, a_m_8, b_m_8, c_m_8, a_t_8, b_t_8, c_t_8, a_w_8, b_w_8, c_w_8)
  
  istat = VGD_OK
  
  return
  
end function check_build_21002

