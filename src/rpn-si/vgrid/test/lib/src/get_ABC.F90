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
program get_invalid

  ! Goal, try to get invalid key and test if status is VGD_ERROR
  ! Use file sigma which is the one having the least valid keys
  
  use, intrinsic :: iso_fortran_env
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_new, vgd_get, &
       vgd_free, VGD_OK, vgd_putopt, VGD_ERROR
  use Unit_Testing, only: ut_report
  implicit none
  integer :: lu=10, fnom, fstouv, fstfrm, fclos
  integer :: i,nk=-99,ier
  type (vgrid_descriptor) :: vgd
  real(kind=REAL64), dimension(:), pointer :: value_8_1D
  character(len=4):: key_S
  character(len=4), dimension(10) :: key_int_S = (/'DATE','IG_1','IG_2','IG_3','IG_4','IP_1','IP_2','IP_3','DIPM','DIPT'/)
  character(len=4), dimension(6) :: key_double_1D_S = (/"CA_T","CB_T","CC_T","CA_W","CB_W","CC_W"/)
  character(len=128) :: file_S
  logical :: OK_L = .true.
  
  nullify(value_8_1D)

  ier = vgd_putopt("ALLOW_SIGMA",.true.)

  file_S="data/dm_1001_from_model_run"
  if( fnom(lu, file_S, 'RND+R/O', 0) < 0 )then
     print*,'ERROR: with fnom on ', trim(file_S)
     error stop
  endif
  if( fstouv(lu, 'RND') == 0 )then
     print*,'ERROR: no record in file ', trim(file_S)
     error stop
  endif
  if( vgd_new(vgd, lu, 'fst') == VGD_ERROR )then
     print*,'Error with vgd_new on file ', trim(file_S)
     error stop
  endif

  print*,'associated(value_8_1D)',associated(value_8_1D)
  
  ! Try to get invalid double 1D
  do i=1,size(key_double_1D_S)
     key_S=key_double_1D_S(i)
     print*
     print*,'========= TESTING ',trim(key_S),'==============='
     ier = vgd_get(vgd, key_double_1D_S(i), value_8_1D)
     if( ier == VGD_ERROR )then
        print*,'The above CRITICAL ERROR is expected'
     else
        OK_L = .false.
     endif
     if( associated(value_8_1D) )then
        print*,'ERROR: vector value_8_1D was allocated and should have not'
        OK_L = .false.
     endif
     if(nk /= -99)then
        print*,'ERROR: value of nk was changed and should have not'
        OK_L = .false.
     endif
  end do
  
  ier = vgd_free(vgd)
  ier = fstfrm(lu)
  ier = fclos(lu)

  if(OK_L)then
      print*,'====================================='
     print*,'ALL TESTS OK'
  endif
  
  call ut_report(VGD_OK,'Try to get get invalid keys')
  
end program get_invalid
