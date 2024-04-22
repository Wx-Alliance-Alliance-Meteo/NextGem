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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  !
  implicit none
  !
  type(vgrid_descriptor) :: d,d2
  integer :: stat,i,test_2001
  integer, dimension(:), pointer :: ip1s
  integer, dimension(:), pointer :: vipm,vipt,work_i
  real(kind=8), dimension(4) :: pres=(/700.,850.,925.,1000./),b=(/0.,0.,0.,0./)
  real(kind=8), dimension(:,:,:), pointer :: tbl
  real(kind=8), dimension(:), pointer :: a_m_8,a_t_8,b_m_8,b_t_8,work_8
  real, dimension(:), pointer :: hyb
  real :: my_real
  logical :: OK
  character(len=1) :: dum_S
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file
  !  
  nullify(ip1s,vipm,vipt,work_i,tbl,a_m_8,a_t_8,b_m_8,b_t_8,work_8,hyb)
  !
  ! Construct a new set of 3D coordinate descriptors
  pres = pres*100. !convert mb to Pa
  allocate(ip1s(size(pres)))
  do i=1,size(ip1s)
     my_real=pres(i)/100.
     call convip(ip1s(i),my_real,2,2,dum_S,.false.)
  enddo
  !
  stat = vgd_new(d,kind=2,version=1,nk=size(pres),ip1_m=ip1s,a_m_8=pres,b_m_8=b)
  stat = vgd_get(d,key='VTBL - vertical coordinate table',value=tbl)
  OK=abs(tbl(2,2,1)-pres(1)) < epsilon(pres)
  !
  file='data/data_constructor_build_2001.txt'
  stat = test_2001(d,file,write_control_L)
  if(stat.eq.VGD_ERROR)OK=.false.
  !
  allocate(hyb(size(pres)))
  hyb=pres/100.
  !
  stat = vgd_new(d2,kind=2,version=1,hyb=hyb)
  !
  if(stat == VGD_ERROR)then
     print*,'Error avec vgd_new'
     error stop 1
  endif
  !
  file='data/data_constructor_build_2001.txt'
  stat = test_2001(d,file,.false.)
  if(stat.eq.VGD_ERROR)OK=.false.  
  !
  call ut_report(OK,'Grid_Descriptors::vgd_new vertical build initializer (2001) value')
  !
end program constructor
!==============================================================================
!===============================================================================
!===============================================================================
integer function test_2001(F_d,F_file,F_write_control_L) result(istat)
   !
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,VGD_ERROR,VGD_OK
   

   implicit none
   !
   type (vgrid_descriptor) :: F_d
   logical :: F_write_control_L
   character (len=256) :: F_file
   !
   ! Local variable
   !
   real, dimension(:), pointer :: vcdm,vcdt,work
   real(kind=8), dimension(:), pointer :: b_m_8,a_m_8,b_t_8,a_t_8,work_8
   integer, dimension(:), pointer :: vipm,vipt,work_i
   integer :: nl_m,nl_t,k,nk,stat
   !
   istat=VGD_ERROR
   !
   nullify(vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8,vipm,vipt,work_i)
   stat = vgd_get(F_d,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8)
   stat = vgd_get(F_d,key='CA_T - vertical A coefficient (t)'   ,value=a_t_8)
   stat = vgd_get(F_d,key='CB_M - vertical B coefficient (m)'   ,value=b_m_8)
   stat = vgd_get(F_d,key='CB_T - vertical B coefficient (t)'   ,value=b_t_8)
   stat = vgd_get(F_d,key='VIPM - level ip1 list (m)'           ,value=vipm)
   stat = vgd_get(F_d,key='VIPT - level ip1 list (t)'           ,value=vipt)
   stat = vgd_get(F_d,key='VCDM - vertical coordinate (m)'      ,value=vcdm)
   stat = vgd_get(F_d,key='VCDT - vertical coordinate (t)'      ,value=vcdt)
   stat = vgd_get(F_d,key='NL_M - Number of vertical levels (m)',value=nl_m)
   stat = vgd_get(F_d,key='NL_T - Number of vertical levels (t)',value=nl_t)
   !
   if(F_write_control_L)then
      open(unit=10,file=F_file)
      write(10,*)a_m_8
      write(10,*)vipm
      write(10,*)vcdm
      close(10)
   endif
   !
   if(nl_m.ne.size(b_m_8))then
      print*,'Problem with NL_M should be',size(b_m_8),' got',nl_m
      error stop 1
   endif
   if(nl_t.ne.size(b_m_8))then
      print*,'Problem with NL_T should be',size(b_m_8),' got',nl_t
      error stop 1
   endif
   !
   nk=size(a_m_8)
   !   
   ! Check A
   open(unit=10,file=F_file,ACTION='READ')
   allocate(work_8(size(a_m_8)))
   print*,'Reading B'
   read(10,*)work_8
   do k=1,nk
      if(abs(work_8(k)-a_m_8(k))/a_m_8(k)>10.*epsilon(a_m_8(k)))then
         print*,'Probleme avec A, pas dans les limites tollerees'
         print*,work_8(k),'vs',a_m_8(k)
         error stop 1
      endif
   enddo
   if(any(a_m_8.ne.a_t_8))then
      print*,'Probleme avec A T, dot egaler A M :'
      print*,a_t_8
      print*,a_m_8
      error stop 1
   endif
   !
   ! Check B
   if(any(b_m_8.ne.0.d0))then
      print*,'Probleme avec B M, doit egaler 0.0d :',b_m_8
      error stop 1
   endif
   if(any(b_t_8.ne.0.d0))then
      print*,'Probleme avec B T, doit egaler 0.0d :',b_t_8
      error stop 1
   endif
   !
   ! Check IP
   allocate(work_i(size(vipm)))
   print*,'Reading IP'
   read(10,*)work_i
   do k=1,nk
      if(work_i(k).ne.vipm(k))then
         print*,'Probleme avec IP :'
         print*,work_i(k),'vs',vipm(k)
         error stop 1
      endif
   enddo
   if(any(vipm.ne.vipt))then
      print*,'Probleme avec ip1 T , doit egaler ip1 M :'
      print*,vipt
      print*,vipm
      error stop 1
   endif
   !
   ! Check vcdm
   allocate(work(size(vcdm)))
   print*,'Reading VCDM'
   read(10,*)work
   do k=1,nk
      if(abs(work(k)-vcdm(k))/vcdm(k)>10.*epsilon(vcdm(k)))then
         print*,'Probleme avec vcdm, pas dans les limites tollerees'
         print*,work(k),'vs',vcdm(k)
         error stop 1
      endif
   enddo
   if(any(vcdm.ne.vcdt))then
      print*,'Probleme avec VCDT, doit egaler VCDM'
      print*,vcdt
      print*,vcdm
      error stop 1
   endif
   !
   close(10)
   !
   istat=VGD_OK
   !
end function test_2001
