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
!===============================================================================

program stda76_hgts_from_pres
#include <rmn/msg.h>
  use Unit_Testing, only: ut_report
  use vGrid_Descriptors, only: vgrid_descriptor, vgd_new, vgd_get, vgd_levels, &
       vgd_stda76_hgts_from_pres_list, VGD_ERROR
  implicit none
  integer :: fnom, fstouv, nl, k
  integer :: lu = 10
  integer, dimension(:), pointer :: ip1s
  real :: p0, ff
  real, dimension(:), pointer :: pres, hgts
  character(len=100) :: file
  type(vgrid_descriptor) :: vgd
  
  nullify(ip1s, pres, hgts)

  file="data/dm_5001_from_model_run"
  
  !Get any pres vertical descriptor
  if(fnom(lu,file,"RND+R/O",0) < 0)then
     print*,'(Test) ERROR with fnom on file ',file
     error stop 1
  endif
  if(fstouv(lu,'RND') < 0)then
     print*,'(Test) No record in RPN file ',file
     error stop 1
  endif
  if(vgd_new(vgd,lu) == VGD_ERROR)then
     print*,'(Test) Problem with vgd_new'
     error stop 1
  endif
  if( vgd_get(vgd, "VIPM", ip1s) ==  VGD_ERROR )then
     print*,"ERROR with Cvgd_get_int for VIPT"
     error stop 1
  end if
  p0 = 105000.
  if(vgd_levels(vgd,sfc_field=p0,ip1_list=ip1s,levels=pres) ==&
       VGD_ERROR)then
     print*,"ERROR with vgd_levels"
     error stop 1
  endif
  allocate(hgts(size(pres)))
  if(vgd_stda76_hgts_from_pres_list(hgts, pres, size(pres)) &
       == VGD_ERROR) error stop 1
  
  ! Data for control is produce by tests c_stda76_hgts_from_pres_list
  open(unit=11, file="data/c_stda76_hgts_from_pres.txt", &
       status='OLD')
  read(11,'(4x,i8)')nl
  if(nl /= size(pres))then
     print*,'In tests, size problem'
     error stop 1
  end if
  do k=1, nl
     read(11,*)ff
     if(abs(ff - hgts(k))/hgts(k) > 100.*epsilon(ff))then
        print*,'OUPS, got ',ff,' expected ',hgts(k)
        error stop 1
     endif
  end do

  call ut_report(.true.,'stda76 hgts from pres list')

end program stda76_hgts_from_pres
