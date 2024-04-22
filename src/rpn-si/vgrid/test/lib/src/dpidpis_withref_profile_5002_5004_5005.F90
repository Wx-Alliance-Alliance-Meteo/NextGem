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

   integer, parameter :: nfile=3
   integer :: stat,lu=0,fnom,fstouv,fstfrm,lutxt=69,ip1,ip2,test_it,i
   character(len=4), dimension(nfile) :: vcode_S=(/"5002","5004","5005"/)   
   logical :: ok=.true.
      
   do i=1,nfile
      
      print*,'Testing ',"data/dm_"//vcode_S(i)//"_from_model_run"           

      lu=10+i
      stat=fnom(lu,"data/dm_"//vcode_S(i)//"_from_model_run","RND",0)
      if(stat.lt.0)then
         print*,'ERROR with fnom',"data/dm_"//vcode_S(i)//"_from_model_run"
         error stop 1
      endif
      stat=fstouv(lu,'RND')
      if(stat.lt.0)then
         print*,'No record in RPN file'
         error stop 1
      endif
      open(unit=lutxt,file='data/dm_'//vcode_S(i)//'_ips.txt',status='OLD')
      read(lutxt,*) ip1,ip2
      close(lutxt)
      
      stat=test_it(lu,ip1,ip2)
      if(stat.ne.VGD_OK)ok=.false.
      
      stat=fstfrm(lu)
      
   enddo

   call ut_report(ok,message='Grid_Descriptors::vgd_dpidpis level calculation status')
      
end program tests


!===============================================================

integer function test_it(F_lu,F_ip1,F_ip2) result(stat)
   
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_dpidpis,vgd_free,VGD_ERROR,VGD_OK
   
   
   implicit none
   
   integer :: F_lu, F_ip1, F_ip2
   
   ! Local variables
   integer :: ni,nj,nk,fstlir,fstinf,k,kind
   integer, parameter :: i0=20,j0=10
   integer, dimension(:), pointer :: ip1_list
   real, dimension(:), pointer :: dpidpis_profil
   real(kind=8), dimension(:), pointer :: dpidpis_profil_8
   real :: epsilon=1.e-5
   real, dimension(:,:), pointer :: p0,px
   real :: local_pres,w1,pres
   real(kind=8) :: local_pres_8
   type(vgrid_descriptor) :: d
   real(kind=8), dimension(:), pointer :: coef_b   

   nullify(ip1_list,dpidpis_profil,dpidpis_profil_8,p0,px,coef_b)

   stat=VGD_ERROR

   ! Get dpidpis
   stat = vgd_new(d,unit=F_lu,format="fst",ip1=F_ip1,ip2=F_ip2)
   if(stat /= VGD_OK)return

   stat = vgd_get(d,key='CB_T - vertical B coefficient (t)',value=coef_b)
   if(stat /= VGD_OK)return

   stat = vgd_get(d,key='VIPT - level ip1 list (t)'        ,value=ip1_list)
   if(stat /= VGD_OK)return

   stat = fstinf(F_lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',"UU")
   if(stat <0)return

   allocate(p0(ni,nj),px(ni,nj))
   
   stat = fstlir(p0,F_lu,ni,nj,nk,-1,'',-1,-1,-1,'','P0')
   if(stat < 0)return
   p0 = p0*100. !mb to Pa
   
   local_pres  =p0(i0,j0)
   local_pres_8=p0(i0,j0)
   
   !===============
   ! Real interface
   
   stat = vgd_dpidpis(d,sfc_field=local_pres,ip1_list=ip1_list,dpidpis=dpidpis_profil)
   if(stat.ne.VGD_OK)then
      print*,'ERROR: problem with vgd_dpidpis'
      return
   endif
   
   print*,'size(dpidpis_profil)',size(dpidpis_profil)
   
   do k=1,size(coef_b)
      call convip(ip1_list(k),pres,kind,-1,"",.false.)
      if(kind.ne.5)cycle
      stat = fstlir(px,F_lu,ni,nj,nk,-1,'',ip1_list(k),-1,-1,'','PX')
      w1=coef_b(k)*px(i0,j0)*100./p0(i0,j0)
      print*,w1,dpidpis_profil(k)
      if(abs(dpidpis_profil(k)- w1)/w1>epsilon)then
         print*,'ERROR; dpidpis_profil do not validate'
         return
      endif
   enddo
   
   !=================
   ! Real*8 interface
   
   stat = vgd_dpidpis(d,sfc_field=local_pres_8,ip1_list=ip1_list,dpidpis=dpidpis_profil_8)
   if(stat.ne.VGD_OK)then
      print*,'ERROR: problem with vgd_dpidpis real(kind=8)'
      return
   endif
   
   print*,'size(dpidpis_profil_8)',size(dpidpis_profil_8)
   
   do k=1,size(coef_b)
      call convip(ip1_list(k),pres,kind,-1,"",.false.)
      if(kind.ne.5)cycle
      stat = fstlir(px,F_lu,ni,nj,nk,-1,'',ip1_list(k),-1,-1,'','PX')
      w1=coef_b(k)*px(i0,j0)*100./p0(i0,j0)
      print*,w1,dpidpis_profil_8(k)
      if(abs(dpidpis_profil_8(k)- w1)/w1>epsilon)then
         print*,'ERROR; dpidpis_profil_8 do not validate'
         return
      endif
   enddo
   
   deallocate(ip1_list,dpidpis_profil,dpidpis_profil_8,p0,px,coef_b)
   stat=vgd_free(d)

   stat=VGD_OK
       
end function test_it
