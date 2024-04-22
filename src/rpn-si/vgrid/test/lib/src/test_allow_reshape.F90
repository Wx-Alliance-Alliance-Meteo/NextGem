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
program tests
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_getopt,vgd_putopt,vgd_print,vgd_levels,vgd_dpidpis,VGD_ERROR,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer :: stat,fnom,fstouv,lu=12,fstkey,fstluk,fstinf,ni,nj,nk,test_functions
  real, dimension(:,:), pointer :: p0
  logical :: ok=.true.

  ! Structure elements
  real, dimension(3) :: hyb=(/0.1,  0.5, 0.9/)
  real, dimension(4) :: hyb2=(/0.1,0.33, 0.66, 0.9/)
  real :: rcoef1=0.,rcoef2=1.
  real(kind=8) :: ptop=1000d0,pref=100000d0
    
  nullify(p0)

  !===================================
  !===================================
  ! Reshaping structure
  !===================================
  !===================================

  print*,'===================================================================================================='
  print*,'Reshaping of structure should not be a problem since it is not subject to the ALLOW_RESHAPE flag'
  print*,'First build';flush(6)
  stat = vgd_new(vgd,kind=5,version=2,hyb=hyb,rcoef1=rcoef1,rcoef2=rcoef2,ptop_8=ptop,pref_8=pref)
  if(stat==VGD_ERROR)then
     print*,'This error should not happen, please fixit'
     error stop 1
  end if
  print*,'Second build';flush(6)
  stat = vgd_new(vgd,kind=5,version=2,hyb=hyb2,rcoef1=rcoef1,rcoef2=rcoef2,ptop_8=ptop,pref_8=pref)
  if(stat==VGD_ERROR)then
     print*,'This error should not happen, since reshaping memory pointed by structure vgrid_descriptor is always allowed'
     error stop 1
  end if

  !===================================
  !===================================
  ! Open lu and get P0 for tests below
  !===================================
  !===================================  
 
  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom on lu'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file data/dm_5002_from_model_run'
     error stop 1
  endif  
  fstkey = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','P0') 
  allocate(p0(ni,nj))
  stat = fstluk(p0,fstkey,ni,nj,nk)
 
  !===================================
  !===================================
  ! Reshaping should produce error
  !===================================
  !===================================

  print*,'===================================================================================================='
  print*,' Reshaping on users pointers is not allowed by default'
  print*,' Test this'

  stat = test_functions(lu,p0,ni,nj,.false.)

  !===================================
  !===================================
  !Reshaping should not produce error
  !===================================
  !===================================

  stat=vgd_putopt('ALLOW_RESHAPE',.true.)
  if(stat==VGD_ERROR)then
     print*,'ERROR on putopt, this should not happend please fixit'
     error stop 1
  endif

  print*,'===================================================================================================='
  print*,' Reshaping on users pointers is now allowed'
  print*,' Test this'

  stat = test_functions(lu,p0,ni,nj,.true.)  

  call ut_report(ok,message='Grid_Descriptors:: test allow reshape')

end program tests
!========================================================================================================
!========================================================================================================
!========================================================================================================
!========================================================================================================
integer function test_functions(lu,p0,ni,nj,reshape_L) result(status)
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_getopt,vgd_putopt,vgd_levels,vgd_dpidpis,vgd_get,VGD_ERROR,VGD_OK
   

   implicit none
   
   integer :: lu,ni,nj
   real, dimension(ni,nj) :: p0
   logical :: reshape_L

   ! Local variables
   integer :: stat,stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size
   integer, dimension(:), pointer :: ip1_list,i1d
   real, dimension(:,:,:), pointer :: lev
   real, dimension(:), pointer :: lev_1d,r1d
   real(kind=8), dimension(:), pointer :: r81d
   real(kind=8), dimension(:,:,:), pointer :: table_8
   type(vgrid_descriptor) :: vgd

   nullify(ip1_list,i1d,lev,lev_1d,r1d,r81d,table_8)

   ! Allocate so they can be dealloca below
   allocate(ip1_list(1),lev(1,1,1),lev_1d(1),r1d(1),r81d(1),table_8(1,1,1))

   status=VGD_ERROR

   ip1_list(1)=93423264
   stat = vgd_new(vgd,lu,'fst')

   print*,'-----------------------------';flush(6)
   print*,'Testing levels_withref';flush(6)
   deallocate(lev)
   stat_not_alloc = vgd_levels(vgd,ip1_list,lev,p0)
   stat_alloc_correct_size = vgd_levels(vgd,ip1_list,lev,p0)
   deallocate(lev);allocate(lev(1,1,1))
   stat_alloc_wrong_size = vgd_levels(vgd,ip1_list,lev,p0)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)

   print*,'----------------------------------';flush(6)
   print*,'Testing levels_withref_prof';flush(6)
   deallocate(lev_1d)
   stat_not_alloc = vgd_levels(vgd,ip1_list,lev_1d,100000.)
   stat_alloc_correct_size = vgd_levels(vgd,ip1_list,lev_1d,100000.)
   deallocate(lev_1d);allocate(lev_1d(2))
   stat_alloc_wrong_size = vgd_levels(vgd,ip1_list,lev_1d,100000.)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)

   print*,'-----------------------------';flush(6)
   print*,'Testing dpidpis_withref';flush(6)
   deallocate(lev)
   stat_not_alloc=vgd_dpidpis(vgd,ip1_list,lev,p0)
   stat_alloc_correct_size=vgd_dpidpis(vgd,ip1_list,lev,p0)
   deallocate(lev);allocate(lev(1,1,1))
   stat_alloc_wrong_size=vgd_dpidpis(vgd,ip1_list,lev,p0)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)
   
   print*,'----------------------------------';flush(6)
   print*,'Testing dpidpis_withref_prof';flush(6)
   deallocate(lev_1d)
   stat_not_alloc=vgd_dpidpis(vgd,ip1_list,lev_1d,100000.)
   stat_alloc_correct_size=vgd_dpidpis(vgd,ip1_list,lev_1d,100000.)
   deallocate(lev_1d);allocate(lev_1d(2)) 
   stat_alloc_wrong_size=vgd_dpidpis(vgd,ip1_list,lev_1d,100000.)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)
   
   print*,'----------------------------------';flush(6)
   print*,'Testing get_allocate_i1d ip1m';flush(6)
   stat_not_alloc=vgd_get(vgd,'VIPM - level ip1 list (m)',i1d)
   stat_alloc_correct_size=vgd_get(vgd,'VIPM - level ip1 list (m)',i1d)
   deallocate(i1d);allocate(i1d(1))
   stat_alloc_wrong_size=vgd_get(vgd,'VIPM - level ip1 list (m)',i1d)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)
   
   print*,'----------------------------------';flush(6)
   print*,'Testing get_allocate_r1d VCDM';flush(6)
   deallocate(r1d)
   stat_not_alloc=vgd_get(vgd,'VCDM - vertical coordinate (m)',r1d)
   stat_alloc_correct_size=vgd_get(vgd,'VCDM - vertical coordinate (m)',r1d)
   deallocate(r1d);allocate(r1d(1))
   stat_alloc_wrong_size=vgd_get(vgd,'VCDM - vertical coordinate (m)',r1d)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)

   print*,'----------------------------------';flush(6)
   print*,'Testing get_allocate_r81d VCDM';flush(6)
   deallocate(r81d)
   stat_not_alloc=vgd_get(vgd,'CA_M - vertical A coefficient (m)',r81d)
   stat_alloc_correct_size=vgd_get(vgd,'CA_M - vertical A coefficient (m)',r81d)
   deallocate(r81d);allocate(r81d(1))
   stat_alloc_wrong_size=vgd_get(vgd,'CA_M - vertical A coefficient (m)',r81d)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)

   print*,'----------------------------------';flush(6)
   print*,'Testing get_allocate_r83d VTBL';flush(6)
   deallocate(table_8)
   stat_not_alloc=vgd_get(vgd,'VTBL - vgrid_descriptor table',table_8)
   stat_alloc_correct_size=vgd_get(vgd,'VTBL - vgrid_descriptor table',table_8)
   deallocate(table_8);allocate(table_8(1,1,1))
   stat_alloc_wrong_size=vgd_get(vgd,'VTBL - vgrid_descriptor table',table_8)
   call check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)   

   deallocate(ip1_list,lev,lev_1d,i1d,r1d,r81d,table_8)

   status=VGD_OK

end function test_functions

subroutine check_status(stat,reshape_L)
   use vGrid_Descriptors, only: VGD_ERROR,VGD_OK
   implicit none
   integer :: stat
   logical :: reshape_L
   if(reshape_L)then
      if(stat==VGD_ERROR)then
         print*,'The above should not produce an error and it does not'
      endif
   else
      if(stat==VGD_ERROR)then
         print*,'The above error is normal'
      else
         print*,'The above should produce an error and does not'
         error stop 1
      endif
   endif   
   flush(6)
   return
end subroutine check_status
subroutine check_status2(stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size,reshape_L)
   use vGrid_Descriptors, only: VGD_ERROR,VGD_OK
   implicit none
   integer :: stat_not_alloc,stat_alloc_correct_size,stat_alloc_wrong_size
   logical :: reshape_L

   if(stat_not_alloc==VGD_ERROR)then
      print*,'The above should not produce an error and it does (not alloc)'
      error stop 1
   endif
   if(stat_alloc_correct_size==VGD_ERROR)then
      print*,'The above should not produce an error and does (alloc correct size)'
      error stop 1
   endif
   if(reshape_L)then
      if(stat_alloc_wrong_size==VGD_ERROR)then
         print*,'The above should not produce an error and it does (alloc wrong size)'
         error stop 1
      endif
   else
      if(stat_alloc_wrong_size==VGD_ERROR)then
         print*,'The above error is normal on (alloc wrong size)'
      else
         print*,'The above should produce an error and does not (alloc wrong size)'
         error stop 1
      endif
   endif   
   flush(6)
   return
end subroutine check_status2
