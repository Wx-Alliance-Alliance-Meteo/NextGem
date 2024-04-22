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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_free,vgd_get,VGD_ERROR
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer, parameter :: nk=27
  integer :: stat,k,nl_m=-1,nl_t=-1
  integer, dimension(:), pointer :: vipm,vipt,work_i
  real :: my_epsilon=1.e-6
  real, dimension(nk) :: sigma,sigma2
  real, dimension(:), pointer :: vcdm,vcdt,work
  real(kind=8), dimension(:), pointer :: b_m_8,a_m_8,b_t_8,a_t_8,work_8
  logical :: OK=.true.
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file

  nullify(vipm,vipt,work_i,vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8)

  sigma=(/0.011000, 0.027000, 0.051000, 0.075000, 0.101000, 0.127000,&
       0.155000, 0.185000, 0.219000, 0.258000, 0.302000, 0.351000,&
       0.405000, 0.460000, 0.516000, 0.574000, 0.631000, 0.688000,&
       0.744000, 0.796000, 0.842000, 0.884000, 0.922000, 0.955000,&
       0.980000, 0.993000, 1.000000/)

  stat = vgd_new(vgd,kind=1,version=1,hyb=sigma)

  stat = vgd_get(vgd,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8)
  stat = vgd_get(vgd,key='CA_T - vertical A coefficient (t)'   ,value=a_t_8)
  stat = vgd_get(vgd,key='CB_M - vertical B coefficient (m)'   ,value=b_m_8)
  stat = vgd_get(vgd,key='CB_T - vertical B coefficient (t)'   ,value=b_t_8)
  stat = vgd_get(vgd,key='VIPM - level ip1 list (m)'           ,value=vipm)
  stat = vgd_get(vgd,key='VIPT - level ip1 list (t)'           ,value=vipt)
  stat = vgd_get(vgd,key='VCDM - vertical coordinate (m)'      ,value=vcdm)
  stat = vgd_get(vgd,key='VCDT - vertical coordinate (t)'      ,value=vcdt)
  stat = vgd_get(vgd,key='NL_M - Number of vertical levels (m)',value=nl_m)
  stat = vgd_get(vgd,key='NL_T - Number of vertical levels (t)',value=nl_t)

  file='data/data_constructor_gen_1001.txt'
  if(write_control_L)then
     open(unit=10,file=file)
     write(10,*)b_m_8
     write(10,*)vipm
     write(10,*)vcdm
     close(10)
  endif

  if(nl_m.ne.size(b_m_8))then
     OK=.false.
     print*,'Problem with NL_M 1001 should be',size(b_m_8),' got',nl_m
  endif
  if(nl_t.ne.size(b_m_8))then
     OK=.false.
     print*,'Problem with NL_T 1001 should be',size(b_m_8),' got',nl_t
  endif

  ! Check A
  if(any(a_m_8.ne.0.d0))then
     OK=.false.
     print*,'Probleme avec A M, doit egaler 0.0d :',a_m_8
  endif
  if(any(a_t_8.ne.0.d0))then
     OK=.false.
     print*,'Probleme avec A T, doit egaler 0.0d :',a_t_8
  endif

  ! Check B
  open(unit=10,file=file,ACTION='READ')
  allocate(work_8(size(b_m_8)))
  print*,'Reading B'
  read(10,*)work_8
  do k=1,nk
     if(abs(work_8(k)-b_m_8(k))/b_m_8(k)>my_epsilon)then
        OK=.false.
        print*,'Probleme avec B, pas dans les limites tollerees'
        print*,work_8(k),'vs',b_m_8(k)
     endif
  enddo

  if(any(b_m_8.ne.b_t_8))then
     OK=.false.
     print*,'Probleme avec B T, dot egaler B M :'
     print*,b_t_8
     print*,b_m_8
  endif

  ! Check IP
  allocate(work_i(size(vipm)))
  print*,'Reading IP'
  read(10,*)work_i
  do k=1,nk
     if(work_i(k).ne.vipm(k))then
        OK=.false.
        print*,'Probleme avec IP :'
        print*,work_i(k),'vs',vipm(k)
     endif
  enddo
  if(any(vipm.ne.vipt))then
     OK=.false.
     print*,'Probleme avec ip1 A , doit egaler ip1 M :'
     print*,vipt
     print*,vipm
  endif

  ! Check vcdm
  allocate(work(size(vcdm)))
  print*,'Reading VCDM'
  read(10,*)work
  do k=1,nk
     if(abs(work(k)-vcdm(k))/vcdm(k)>my_epsilon)then
        OK=.false.
        print*,'Probleme avec vcdm, pas dans les limites tollerees'
        print*,work(k),'vs',vcdm(k)
     endif
  enddo  
  if(any(vcdm.ne.vcdt))then
     OK=.false.
     print*,'Probleme avec VCDT, doit egaler VCDM'
     print*,vcdt
     print*,vcdm
  endif
  
  print*,'The following error is normal'
  print*,'Testing last sig level not at 1.0'
  stat = vgd_new(vgd,kind=1,version=1,hyb=sigma(1:size(sigma)-1))
  if(stat.ne.VGD_ERROR)OK=.false.

  print*,'The following error is normal'
  print*,'Testing sig levels not in order'
  sigma2=sigma
  sigma2(size(sigma)-3)=sigma2(size(sigma)-2)
  stat = vgd_new(vgd,kind=1,version=1,hyb=sigma2)
  if(stat.ne.VGD_ERROR)OK=.false.
  deallocate(a_m_8,a_t_8,b_m_8,b_t_8,vipm,vipt,vcdm)
  stat = vgd_free(vgd)

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer (1001) value')

end program constructor
