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
  
  use, intrinsic :: iso_fortran_env
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,vgd_levels,VGD_ERROR, VGD_NO_REF_NOMVAR
  use Unit_Testing, only: ut_report
  
  implicit none

  type(vgrid_descriptor) :: vgd
  integer :: i,k,nk_abc,ier
  integer, dimension(:), pointer :: ip1_list
  integer, parameter :: nk = 50
  integer, dimension(nk) :: ip1_m, ip1_w
  real :: work  
  real, dimension(8,4) :: just_for_hor_size
  real, dimension(:,:,:), pointer :: levels
  real(kind=8), dimension(nk) :: hgt_w= (/&
       0.0000,     1.0113,     2.0857,     3.2230,     4.4372,     5.7451,     7.1677,     8.7300,    10.4630,    12.4044,&
      14.5999,    17.1056,    19.9897,    23.3350,    27.2426,    31.8354,    37.2625,    43.7049,    51.3819,    60.5588,&
      71.5555,    84.7573,   100.6260,   119.7140,   142.6780,   170.2940,   203.4740,   243.2780,   290.9300,   347.8220,&
     415.5120,   495.7180,   590.2900,   701.1660,   830.3120,   979.6430,  1150.9300,  1345.6900,  1565.0900,  1809.8800,&
    2080.3100,  2376.1200,  2696.5800,  3040.5200,  3406.4000,  3792.4700,  4196.8200,  4617.4600,  5052.4700,  5500.0000/)
  real(kind=8), dimension(nk) :: hgt_m= (/&
          0.4940,     1.5414,     2.6457,     3.8195,     5.0782,     6.4406,     7.9296,     9.5730,    11.4050,    13.4671, &
         15.8101,    18.4956,    21.5988,    25.2114,    29.4447,    34.4342,    40.3441,    47.3737,    55.7643,    65.8073,&
         77.8539,    92.3261,   109.7290,   130.6660,   155.8510,   186.1260,   222.4750,   266.0400,   318.1270,   380.2130,&
        453.9380,   541.0890,   643.5670,   763.3330,   902.3390,  1062.4400,  1245.2900,  1452.2500,  1684.2800,  1941.8900,&
       2225.0800,  2533.3400,  2865.7000,  3220.8200,  3597.0300,  3992.4800,  4405.2200,  4833.2900,  5274.7800,  5727.9200/)
  real(kind=8), dimension(nk) :: b_m_8, b_w_8
  real(kind=REAL64), dimension(:), pointer :: value_1d_8
  logical, parameter :: write_control_L=.false.
  logical :: OK
  character(len=1) :: dum_S
  integer, parameter :: nkeys=9
  character(len=4), dimension(9) :: keys_S = (/"CA_M","CB_M","CC_M","CA_T","CB_T","CC_T","CA_W","CB_W","CC_W"/)
  character(len=4) :: rfls_S, key_nk
  
  nullify(ip1_list,levels,value_1d_8)
  
  OK = .true.
  just_for_hor_size = 0.
  b_m_8=0.
  b_w_8=0.
  
  do k=1,nk
     work=hgt_w(k)
     call convip( ip1_w(k), work, 0, 2,dum_S, .false.)     
     work=hgt_m(k)
     call convip( ip1_m(k), work, 0, 2,dum_S, .false.)
  enddo

  ! Build a new set of vertical coordinate descriptors Vcode 0001
  if( vgd_new(vgd,0,1,nk,0,0,a_m_8=hgt_m,a_w_8=hgt_w,b_m_8=b_m_8,b_w_8=b_w_8,ip1_m=ip1_m,ip1_w=ip1_w) == VGD_ERROR )then
     print*,'ERROR in tests with vgd_new'
     error stop 1
  endif

  k = vgd_print(vgd)

  print*,'Testing A M'
  if( vgd_get(vgd,key="VIPM",value=ip1_list) == VGD_ERROR )then
     print*,'ERROR in test with vgd_get on key VIPM'
     error stop 1
  endif  
  if( vgd_levels(vgd,ip1_list,levels) == VGD_ERROR )then
     print*,'ERROR in test with vgd_levels'
     error stop 1
  endif
  do k=1,nk
     if(abs(levels(1,1,k)) < epsilon(1.))then
        if( abs(levels(1,1,k)-hgt_m(k)) > 100.*epsilon(1.) )then
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,levels(1,1,k),'vs',hgt_m(k)
        endif
     else
        if( abs(levels(1,1,k)-hgt_m(k)) / abs(levels(1,1,k)) > 100.*epsilon(1.) )then
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,levels(1,1,k),'vs',hgt_m(k)
        endif
     endif
  enddo
  
  print*,'Testing A W'
  if( vgd_get(vgd,key="VIPW",value=ip1_list) == VGD_ERROR )then
     print*,'ERROR in test with vgd_get on key VIPW'
     error stop 1
  endif  
  if( vgd_levels(vgd,ip1_list,levels) == VGD_ERROR )then
     print*,'ERROR in test with vgd_levels'
     error stop 1
  endif
  do k=1,nk
     if(abs(levels(1,1,k)) < epsilon(1.))then
        if( abs(levels(1,1,k)-hgt_w(k)) > 100.*epsilon(1.) )then
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,levels(1,1,k),'vs',hgt_w(k)
        endif
     else
        if( abs(levels(1,1,k)-hgt_w(k)) / abs(levels(1,1,k)) > 100.*epsilon(1.) )then
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,levels(1,1,k),'vs',hgt_w(k)
        endif
     endif
  enddo

  ! Try to retreive all A, B, C
  do i=1,nkeys
     print*,'Checking key ',keys_S(i)
     flush(6)
     if(vgd_get(vgd,keys_S(i),value_1d_8) == VGD_ERROR)then
        print*,'Problem in test constructor, s/r check_get_ABC, cannot get key ',keys_S(i)       
        error stop
     endif
     if(.not. associated(value_1d_8))then
        print*,'Problem in test constructor, s/r check_get_ABC, value_1d_8 is not assiciated for ',keys_S(i)
        error stop
     endif
     ! Test that B and C coefs are zero
     if( keys_S(i)(2:2) == "B" .or. keys_S(i)(2:2) == "C")then
        print*,'Testing ',keys_S(i),' for zero value'
        !NL_M, NL_T or NL_W
        key_nk="NL_"//keys_S(i)(4:4)
        ier = vgd_get(vgd,key_nk,nk_abc)
        do k=1,nk_abc
           if(value_1d_8(k) /= 0.d0)then             
              print*,'Problem in test constructor, s/r check_get_ABC, some value_1d_8 are not zero and should for key ',keys_S(i)
              print*,value_1d_8(k)
              error stop
           endif
        end do
     endif
     deallocate(value_1d_8)
  end do

  call ut_report(OK,'Grid_Descriptors::vgd_levels height computation for ocean coordinate')
end program constructor
