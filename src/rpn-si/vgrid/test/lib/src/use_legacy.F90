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
program use_legacy
   use Unit_Testing, only: ut_report
   use vGrid_Descriptors, only : vgd_putopt, VGD_ERROR
   implicit none
   integer :: istat,test_it
   logical :: OK=.true.

   istat = vgd_putopt("ALLOW_SIGMA",.true.)

   istat=test_it("data/dm_1001_from_model_run")
   if(istat.eq.VGD_ERROR) error stop 1

   istat=test_it("data/dm_1002_from_model_run_multi_pt")
   if(istat.eq.VGD_ERROR) error stop 1

   istat=test_it("data/dm_1002_from_model_run_ip1s_new_style")
   if(istat.eq.VGD_ERROR) error stop 1
   
   istat=test_it("data/dm_1003_from_pgsm_lam_east_ops")
   if(istat.eq.VGD_ERROR) error stop 1

   istat=test_it("data/dm_5001_from_model_run")
   if(istat.eq.VGD_ERROR) error stop 1

   call ut_report(OK,'Grid_Descriptors::vgd_vintage (1001 sigma)')   
end program use_legacy
!=========================================================================
!=========================================================================
integer function test_it(F_file) result(istat)
   
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,vgd_print,vgd_levels,VGD_OK,VGD_ERROR

   implicit none
   character(len=*), intent(in) :: F_file

   ! Local Variables
   type(vgrid_descriptor) :: vgd
   integer, parameter :: nmax=1000
   integer, dimension(nmax) :: liste
   integer :: fnom,fstouv,fstfrm,fstinl,ni,nj,nk,infon,fstluk,k,fclos
   integer :: stat,lu=10
   real, dimension(:,:,:), pointer :: levels,px
   real, parameter :: epsilon=6.0e-6
   real :: ff
      
   nullify(levels,px)

   istat=VGD_ERROR

   print*,'========================================================'
   print*,'========================================================'
   print*,'========================================================'
   print*,'========================================================'
   print*,'========================================================'
   print*,'========================================================'
   print*,'Testing file ',trim(F_file)

   stat=fnom(lu,F_file,"RND+R/O",0)
   if(stat.lt.0)then
      print*,'ERROR with fnom'
      return
   endif
   stat=fstouv(lu,"RND")
   if(stat.lt.0)then
      print*,'Error with fstouv'
      return
   endif
   
   stat=vgd_new(vgd,lu)
   if(stat.eq.VGD_ERROR)return
  
   stat=vgd_print(vgd,-1)  
   if(stat.eq.VGD_ERROR)return

   stat=fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','PX',liste,infon,nmax)
  
   allocate(levels(ni,nj,infon),px(ni,nj,infon))

   stat=vgd_levels(vgd,lu,liste(1:infon),levels)
   if(stat.eq.VGD_ERROR)return
   
   do k=1,infon
      stat=fstluk(px(1,1,k),liste(k),ni,nj,nk)
      if(stat.lt.0)return
   enddo
   
   do k=1,infon
      ff=abs(px(1,1,k)-levels(1,1,k)/100.)/px(1,1,k)
      if(ff.gt.epsilon)then        
         print*,'OUPS, these should be equal within epsilon',px(1,1,k),levels(1,1,k)/100.
         print*,'But found difference of',ff,' for k =',k
         return
      endif
   enddo
   
   stat=fstfrm(lu)
   stat=fclos(lu)
   
   lu=lu+1

   istat=VGD_OK

end function test_it

