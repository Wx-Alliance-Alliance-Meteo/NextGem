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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels,vgd_get,vgd_print,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer, pointer, dimension(:) :: ip1s
  integer :: stat
  real, dimension(57) :: hyb= &
     (/0.0134575, 0.0203980, 0.0333528, 0.0472815, 0.0605295, 0.0720790, &
       0.0815451, 0.0889716, 0.0946203, 0.0990605, 0.1033873, 0.1081924, &
       0.1135445, 0.1195212, 0.1262188, 0.1337473, 0.1422414, 0.1518590, &
       0.1627942, 0.1752782, 0.1895965, 0.2058610, 0.2229843, 0.2409671, &
       0.2598105, 0.2795097, 0.3000605, 0.3214531, 0.3436766, 0.3667171, &
       0.3905587, 0.4151826, 0.4405679, 0.4666930, 0.4935319, 0.5210579, &
       0.5492443, 0.5780612, 0.6074771, 0.6374610, 0.6679783, 0.6989974, &
       0.7299818, 0.7591944, 0.7866292, 0.8123021, 0.8362498, 0.8585219, &
       0.8791828, 0.8983018, 0.9159565, 0.9322280, 0.9471967, 0.9609448, &
       0.9735557, 0.9851275, 0.9950425/)
  real, dimension(2,1) :: p0, p0l
  real :: rcoef1=0.,rcoef2=4.,rcoef3=0.,rcoef4=200.
  
  real, pointer, dimension(:) :: levels
  real, pointer, dimension(:,:,:) :: levels_3d

  real(kind=8) :: ptop=805d0,pref=100000d0
  logical :: OK=.true.
  integer :: test_5100
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file
  
  nullify(ip1s,levels,levels_3d)

  p0(:,1) =(/100000.,50000./)
  p0l(:,1)=(/100000.,90000./)

  ! Call with missing rcoef3 and rcoef4
  print*,'Following error on rcoef3 and rcoef4 is expected'
  if( vgd_new(vgd,kind=5,version=100,hyb=hyb,rcoef1=rcoef1,rcoef2=rcoef2,pref_8=pref,dhm=10.0,dht=2.0,ptop_out_8=ptop,avg_L=.true.)&
       /= VGD_ERROR)OK=.false.

  ! Construct a new set of vertical coordinate descriptors 5100
  stat = vgd_new(vgd,kind=5,version=100,hyb=hyb,rcoef1=rcoef1,rcoef2=rcoef2,rcoef3=rcoef3,rcoef4=rcoef4,pref_8=pref,dhm=10.0,&
       dht=2.0,ptop_out_8=ptop,avg_L=.true.)
  stat = vgd_print(vgd)
  file='data/data_constructor_gen_5100_avg.txt'
  stat = test_5100(vgd,file,write_control_L,stat)
  if(stat.eq.VGD_ERROR)OK=.false.

  stat = vgd_new(vgd,kind=5,version=100,hyb=hyb,rcoef1=rcoef1,rcoef2=rcoef2,rcoef3=rcoef3,rcoef4=rcoef4,pref_8=pref,dhm=10.0,&
       dht=2.0,ptop_out_8=ptop,avg_L=.false.)
  stat = vgd_print(vgd)
  file='data/data_constructor_gen_5100.txt'
  stat = test_5100(vgd,file,write_control_L,stat)
  if(stat.eq.VGD_ERROR)OK=.false.  

  if(.false.)then
     stat = vgd_get(vgd,'VIPM',ip1s)
     print*,'ip1s=',ip1s
     stat = vgd_levels(vgd,ip1s,levels,sfc_field=100000.,in_log=.false.,sfc_field_ls=100000.)
     print*,'levels=',levels
     
     stat = vgd_get(vgd,'VIPT',ip1s)
     print*,'ip1s=',ip1s
     stat = vgd_levels(vgd,ip1s,levels,sfc_field=100000.,in_log=.false.,sfc_field_ls=100000.)
     print*,'levels=',levels
     stat = vgd_levels(vgd,ip1s,levels_3d,sfc_field=p0,in_log=.false.,sfc_field_ls=p0l)
     print*,'levels_3d(1,1,1:size(levels_3d,dim=3))=',levels_3d(1,1,1:size(levels_3d,dim=3))
     print*,'levels_3d(2,1,1:size(levels_3d,dim=3))=',levels_3d(2,1,1:size(levels_3d,dim=3))
  endif

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer (5100) value')
end program constructor
!==============================================================================
!===============================================================================
!===============================================================================
integer function test_5100(F_d,F_file,F_write_control_L,F_stat) result(istat)
   !
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,VGD_ERROR,VGD_OK
   

   implicit none
   !
   type (vgrid_descriptor) :: F_d
   logical :: F_write_control_L
   character (len=256) :: F_file
   integer :: F_stat
   !
   ! Local variable
   !
   real, dimension(:), pointer :: vcdm,vcdt,work
   real(kind=8), dimension(:), pointer :: b_m_8,c_m_8,a_m_8,b_t_8,c_t_8,a_t_8,work_8
   integer, dimension(:), pointer :: vipm,vipt,work_i
   integer :: nl_m,nl_t,k,nk,kind,vers,ip1,my_ip1

   istat=VGD_OK

   flush(6)

   if(F_stat.eq.VGD_ERROR)then
      print*,'In test_5100: test not performed, since previous command failed'
      istat=VGD_ERROR
      return
   endif

   nullify(vcdm,vcdt,work,b_m_8,c_m_8,a_m_8,b_t_8,c_t_8,a_t_8,work_8,vipm,vipt,work_i)

   if( vgd_get(F_d,key='KIND - vertical coordinate ip1 kind' ,value=kind)   == VGD_ERROR) return
   if( vgd_get(F_d,key='VERS - vertical coordinate version'  ,value=vers)   == VGD_ERROR) return
   if( vgd_get(F_d,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CA_T - vertical A coefficient (t)'   ,value=a_t_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CB_M - vertical B coefficient (m)'   ,value=b_m_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CC_M - vertical C coefficient (m)'   ,value=c_m_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CB_T - vertical B coefficient (t)'   ,value=b_t_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CC_T - vertical C coefficient (t)'   ,value=c_t_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='VIPM - level ip1 list (m)'           ,value=vipm)   == VGD_ERROR) return
   if( vgd_get(F_d,key='VIPT - level ip1 list (t)'           ,value=vipt)   == VGD_ERROR) return
   if( vgd_get(F_d,key='VCDM - vertical coordinate (m)'      ,value=vcdm)   == VGD_ERROR) return
   if( vgd_get(F_d,key='VCDT - vertical coordinate (t)'      ,value=vcdt)   == VGD_ERROR) return
   if( vgd_get(F_d,key='NL_M - Number of vertical levels (m)',value=nl_m)   == VGD_ERROR) return
   if( vgd_get(F_d,key='NL_T - Number of vertical levels (t)',value=nl_t)   == VGD_ERROR) return
   if( vgd_get(F_d,key='IP_1 - record ip1'                   ,value=ip1)    == VGD_ERROR) return
   print*,'TESTING ',kind*1000+vers

   if(F_write_control_L)then
     open(unit=10,file=F_file)
     write(10,*)a_m_8
     write(10,*)a_t_8
     write(10,*)b_m_8
     write(10,*)b_t_8
     write(10,*)c_m_8
     write(10,*)c_t_8
     write(10,*)vipm
     write(10,*)vipt
     write(10,*)vcdm
     write(10,*)vcdt
     write(10,*)ip1
     close(10)
  endif

  if(nl_m.ne.size(b_m_8))then
     print*,'wrong size with NL_M, got',nl_m,' should be ',size(b_m_8)
     istat=VGD_ERROR
  endif
  if(nl_t.ne.size(b_t_8))then
     print*,'wrong size with NL_T, got',nl_t,' should be ',size(b_t_8)
     istat=VGD_ERROR
  endif

  nk=nl_m

  ! Check A
  open(unit=10,file=F_file,ACTION='READ')
  allocate(work_8(size(a_m_8)))
  print*,'Reading A M'
  read(10,*)work_8

  do k=1,nk
     if(abs(work_8(k)-a_m_8(k))/a_m_8(k)>100.*epsilon(a_m_8(k)))then
        istat=VGD_ERROR
        print*,'Probleme avec A M, pas dans les limites tollerees, k=',k
        print*,work_8(k),'vs'
        print*,a_m_8(k)
     endif
  enddo

  print*,'Reading A T'
  read(10,*)work_8
  do k=1,nk
     if(abs(work_8(k)-a_t_8(k))/a_t_8(k)>100.*epsilon(1.))then
        istat=VGD_ERROR
        print*,'Probleme avec A T, pas dans les limites tollerees, k=',k
        print*,work_8(k),'vs'
        print*,a_t_8(k)
     endif
  enddo
 
  ! Check B
  print*,'Reading B M'
  read(10,*)work_8
  do k=1,nk
     if(b_m_8(k).eq.0.)then
        if(work_8(k).ne.0.)then
           istat=VGD_ERROR
           print*,'Probleme avec B M, pas egal a zero ',work_8(k),' pour k=',k
        endif
     else
        if(abs(work_8(k)-b_m_8(k))/b_m_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B M, pas dans les limites tollerees, k=',k
           print*,work_8(k),'vs'
           print*,b_m_8(k)
        endif
     endif
  enddo
  print*,'Reading B T'
  read(10,*)work_8
  do k=1,nk
     if(b_t_8(k).eq.0.)then
        if(work_8(k).ne.0.)then
           istat=VGD_ERROR
           print*,'Probleme avec B T, pas egal a zero ',work_8(k),' pour k=',k
        endif
     else
        if(abs(work_8(k)-b_t_8(k))/b_t_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B T, pas dans les limites tollerees, k=',k
           print*,work_8(k),'vs'
           print*,b_t_8(k)
        endif
     endif
  enddo
  
  ! Check C
  print*,'Reading C M'
  read(10,*)work_8
  do k=1,nk
     if(c_m_8(k).eq.0.)then
        if(work_8(k).ne.0.)then
           istat=VGD_ERROR
           print*,'Probleme avec C M, pas egal a zero ',work_8(k),' pour k=',k
        endif
     else
        if(abs(work_8(k)-c_m_8(k))/c_m_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec C M, pas dans les limites tollerees, k=',k
           print*,work_8(k),'vs'
           print*,c_m_8(k)
        endif
     endif
  enddo
  print*,'Reading C T'
  read(10,*)work_8
  do k=1,nk
     if(c_t_8(k).eq.0.)then
        if(work_8(k).ne.0.)then
           istat=VGD_ERROR
           print*,'Probleme avec C t, pas egal a zero ',work_8(k),' pour k=',k
        endif
     else 
        if(abs(work_8(k)-c_t_8(k))/c_t_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec C T, pas dans les limites tollerees, k=',k
           print*,work_8(k),'vs'
           print*,c_t_8(k)
        endif
     endif
  enddo

  ! Check IPs
  allocate(work_i(size(vipm)))
  print*,'Reading IP M'
  read(10,*)work_i
  do k=1,nk
     if(work_i(k).ne.vipm(k))then
        istat=VGD_ERROR
        print*,'Probleme avec IP pour k=',k,' :'
        print*,work_i(k),'vs'
        print*,vipm(k)
     endif
  enddo
  print*,'Reading IP T'
  read(10,*)work_i
  do k=1,nk
     if(work_i(k).ne.vipt(k))then
        istat=VGD_ERROR
        print*,'Probleme avec IP pour k=',k,' :'
        print*,work_i(k),'vs'
        print*,vipt(k)
     endif
  enddo

  ! Check vcdm
  allocate(work(size(vcdm)))
  print*,'Reading VCDM'
  read(10,*)work
  do k=1,nk
     if(vcdm(k).eq.0.)then
        if(abs(work(k)-vcdm(k))>100.*epsilon(vcdm(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees, k=',k
           print*,work(k),'vs'
           print*,vcdm(k)
        endif
     else
        if(abs(work(k)-vcdm(k))/vcdm(k)>100.*epsilon(vcdm(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees, k=',k
           print*,work(k),'vs'
           print*,vcdm(k)
        endif
     endif
  enddo

  print*,'Reading VCDT'
  read(10,*)work
  do k=1,nk
     if(vcdt(k).eq.0.)then
        if(abs(work(k)-vcdt(k))>100.*epsilon(vcdt(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees, k=',k
           print*,work(k),'vs'
           print*,vcdt(k)
        endif
     else
        if(abs(work(k)-vcdt(k))/vcdt(k)>100.*epsilon(vcdt(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdt, pas dans les limites tollerees, k=',k
           print*,work(k),'vs'
           print*,vcdt(k)
        endif
     endif
  enddo  

  print*,'Reading VCDT'
  read(10,*)my_ip1
  if(my_ip1.ne.ip1)then
      istat=VGD_ERROR
     print*,'Probleme avec ip1, expected ',my_ip1,' got',ip1
  endif

  close(10)

  deallocate(vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8,vipm,vipt,work_i)

end function test_5100
