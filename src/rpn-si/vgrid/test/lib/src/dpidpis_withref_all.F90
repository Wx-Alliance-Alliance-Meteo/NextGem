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
  use Unit_Testing, only: ut_report

  use vGrid_Descriptors, only: VGD_OK

  implicit none

  integer :: stat,lu=0,fnom,fstouv,fstfrm,test_dpidpis,i
  
  integer, parameter :: nfile=4

  character(len=4), dimension(nfile) :: vcode_S=(/"5002","5004","5005","5100"/)

  logical :: ok=.true.

  do i=1,nfile

     print*,'Testing ',"data/dm_"//vcode_S(i)//"_from_model_run"

     lu=10+i
     stat=fnom(lu,"data/dm_"//vcode_S(i)//"_from_model_run","RND",0)
     if(stat.lt.0)then
        print*,'ERROR with fnom on',"data/dm_"//vcode_S(i)//"_from_model_run"
        error stop 1
     endif
     stat=fstouv(lu,'RND')
     if(stat.lt.0)then
        print*,'No record in RPN file'
        error stop 1
     endif
     stat=test_dpidpis(lu)
     if(stat.ne.VGD_OK)ok=.false.

  enddo

  call ut_report(ok,message='Grid_Descriptors::vgd_dpidpis level calculation status')

  stat=fstfrm(lu)

end program tests


!=============================================================================================

integer function test_dpidpis(F_lu) result(stat)

   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_dpidpis,vgd_free,VGD_ERROR,VGD_OK
   

   implicit none
   
   integer :: F_lu
   
   ! Local variables
   integer :: ni,nj,nk,fstlir,fstinf,k,kind
   integer, parameter :: i0=20,j0=10
   integer, dimension(:), pointer :: ip1_list
   real, dimension(:,:,:), pointer :: dpidpis_cube
   real(kind=8), dimension(:,:,:), pointer :: dpidpis_cube_8
   real   :: epsilon  =1.e-6
   real, dimension(:,:), pointer :: p0,px
   real(kind=8), dimension(:,:), pointer :: p0_8
   real   :: w1,pres
   type(vgrid_descriptor) :: d
   logical :: ok
   real(kind=8), dimension(:), pointer :: coef_b
  
   nullify(ip1_list,dpidpis_cube,dpidpis_cube_8,p0,px,p0_8,coef_b)
  
   stat=VGD_ERROR

   ! Get dpidpis
   stat = vgd_new(d,unit=F_lu,format="fst")
   if(stat /= VGD_OK)return

   stat = vgd_get(d,key='CB_T - vertical B coefficient (t)',value=coef_b)
   if(stat /= VGD_OK)return
   
   stat = vgd_get(d,key='VIPT - level ip1 list (t)'        ,value=ip1_list)
   if(stat /= VGD_OK)return
   
   stat = fstinf(F_lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',"UU")
   if(stat <0)return
   
   allocate(p0(ni,nj),px(ni,nj),p0_8(ni,nj))
   
   stat = fstlir(p0,F_lu,ni,nj,nk,-1,'',-1,-1,-1,'','P0')
   if(stat < 0)return

   p0  =p0*100. !mb to Pa
   p0_8=p0
   
   !===============
   ! Real interface
   
   stat = vgd_dpidpis(d,sfc_field=p0,ip1_list=ip1_list,dpidpis=dpidpis_cube)
   if(stat.ne.VGD_OK)then
      print*,'ERROR: problem with vgd_dpidpis real'
      return
   endif
   print*,'size(dpidpis_cube)',size(dpidpis_cube)
   
   OK=.true.
   do k=1,size(coef_b)
      call convip(ip1_list(k),pres,kind,-1,"",.false.)
      if(kind.ne.5)cycle
      stat = fstlir(px,F_lu,ni,nj,nk,-1,'',ip1_list(k),-1,-1,'','PX')
      if(stat < 0)then
         print*,'ERROR problem with fstlir on PX for ip1 ->',ip1_list(k)
         return
      endif
      !
      ! If the difference in the value of the derivtive (dpidpis) is smaller than epsilon*p
      ! we consider this ok
      ! 
      w1=coef_b(k)*px(i0,j0)*100./p0(i0,j0)
      ! Cannot test zero so take a very small real value
      if( abs(w1) < 1.E-37 )then
         if ( abs(dpidpis_cube(i0,j0,k)) > 1.E-37 )then
            print*,'vgd_dpidpis real do not validate, expect', w1,' got ',dpidpis_cube(i0,j0,k)
            return
         endif
      else
         print*,'TOTO',abs(dpidpis_cube(i0,j0,k)- w1), px(i0,j0)*100.*epsilon
         if(abs(dpidpis_cube(i0,j0,k)- w1) > px(i0,j0)*100.*epsilon)then
            print*,'vgd_dpidpis real do not validate, expect', w1,' got ',dpidpis_cube(i0,j0,k),abs(dpidpis_cube(i0,j0,k)- w1)/w1, &
                 epsilon
         endif
      endif
   enddo
   
   !=================
   ! Real*8 interface
   
   stat = vgd_dpidpis(d,sfc_field=p0_8,ip1_list=ip1_list,dpidpis=dpidpis_cube_8)
   if(stat.ne.VGD_OK)then
      print*,'ERROR: problem with vgd_dpidpis real(kind=8)'
      return
   endif
   
   print*,'size(dpidpis_cube_8)',size(dpidpis_cube_8)
   
   OK=.true.
   do k=1,size(coef_b)
      call convip(ip1_list(k),pres,kind,-1,"",.false.)
      if(kind.ne.5)cycle
      stat = fstlir(px,F_lu,ni,nj,nk,-1,'',ip1_list(k),-1,-1,'','PX')
      if(stat < 0)then
         print*,'ERROR problem with fstlir on PX for ip1 ->',ip1_list(k)
         return
      endif
      w1=coef_b(k)*px(i0,j0)*100./p0_8(i0,j0)
      ! Note since px is a real the precision cannot be hier than real
      if( abs(w1) < 1.E-37)then
         if ( abs(dpidpis_cube_8(i0,j0,k)) > 1.d0*1.0E-37)then
            print*,'vgd_dpidpis real(kind=8) do not validate, expect', w1,' got ',dpidpis_cube_8(i0,j0,k)
            return
         endif
      else
         if(abs(dpidpis_cube_8(i0,j0,k)- w1) > px(i0,j0)*100.d0*epsilon)then
            print*,'vgd_dpidpis real(kind=8) do not validate, expect', w1,' got ',dpidpis_cube_8(i0,j0,k)
            return
         endif
      endif
   enddo

   deallocate(ip1_list,dpidpis_cube,dpidpis_cube_8,p0,px,p0_8,coef_b)
   
   stat=vgd_free(d)

   stat=VGD_OK

end function test_dpidpis
