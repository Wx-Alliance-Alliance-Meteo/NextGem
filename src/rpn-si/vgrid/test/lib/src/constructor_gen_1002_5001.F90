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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_putopt, vgd_free,VGD_ERROR
  use Unit_Testing, only: ut_report
  
   
  implicit none

  type(vgrid_descriptor) :: d
  integer :: stat
  real :: my_epsilon=1.e-6
  !From /home/binops/afsi/sio/datafiles/constants/modeles/GEMDM/glb_033L80_g1_20090922/gem_settings.nml
  real, dimension(28) :: hyb2,hyb=(/&
       0.0125000,0.0233625,0.0391625,0.0628625,0.0865625,&
       0.1122375,0.1379125,0.1655625,0.1951875,0.2287625,&
       0.2672750,0.3107250,0.3591125,0.4124375,0.4667500,&
       0.5220500,0.5793250,0.6356125,0.6919000,0.7472000,&
       0.7985500,0.8439750,0.8854500,0.9229750,0.9555625,&
       0.9802499,0.9930875,1.0000000&
       /)
  real, dimension(28) :: hyb_N=(/&
       0.000,   0.011,    0.027,    0.051,    0.075,&
       0.101,   0.127,    0.155,    0.185,    0.219,&
       0.258,   0.302,    0.351,    0.405,    0.460,&
       0.516,   0.574,    0.631,    0.688,    0.744,&
       0.796,   0.842,    0.884,    0.922,    0.955,&
       0.980,   0.993,    1.000&
       /)
  real :: rcoef1=1.6
  real(kind=8) :: ptop=-1,ptop_N=1000.d0,pref=80000.d0,w1
  real(kind=8), dimension(:), pointer :: b_m_8,a_m_8
  logical :: OK=.true.
  integer :: test_1002_5001
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file

  nullify(b_m_8,a_m_8)

  stat = vgd_putopt('ALLOW_RESHAPE',.true.)

  !
  ! Testing 1002
  !    NOTE : le changement du calcul des A et des B pour 1002 fait
  !           que les donnees de control n<ont que 4 chiffre d'identiques car
  !           les intrants de ces calculs sont a 32 bits.
  !
  print*,''
  print*,'===== 1002 Test on Error ptop<0, error message is normal ===='
  stat = vgd_new(d,kind=1,version=2,hyb=hyb_N,ptop_8=-1.d0)
  if(stat.ne.VGD_ERROR)OK=.false.  
  !
  print*,''  
  print*,'===== 1002 Test on complete level set ===='
  stat = vgd_new(d,kind=1,version=2,hyb=hyb_N,ptop_8=ptop_N)
  file='data/data_constructor_gen_1002.txt'
  stat = test_1002_5001(d,file,write_control_L)
  if(stat.eq.VGD_ERROR)OK=.false.
  print*,''  
  print*,'===== 1002 Test on incomplete level set ===='
  stat = vgd_new(d,kind=1,version=2,hyb=hyb_N(10:28),ptop_8=ptop_N)
  stat = vgd_get(d,key='CB_M - vertical B coefficient',value=b_m_8)
  stat = vgd_get(d,key='CA_M - vertical A coefficient',value=a_m_8)
  w1=0.2189999967813492
  if(.not.(abs(b_m_8(1)-w1)/w1<my_epsilon))then
     OK=.false.
     print*,'Probleme avec B 1002 (incomplete level set)'
     print*,b_m_8(1),' should equal ',w1
  endif
  w1=780.9996604919434
  if(.not.(abs(a_m_8(1)-w1)/w1<my_epsilon))then
     OK=.false.
     print*,'Probleme avec A 1002 (incomplete level set)'
     print*,a_m_8(1),' should equal ',w1
  endif  
  !
  ! Testing 5001
  !
  print*,''
  print*,'===== Test on Error ptop>hyb(1), error message is normal ===='
  ptop=100000.

  stat = vgd_new(d,kind=5,version=1,hyb=hyb,rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat.ne.VGD_ERROR)OK=.false.
  print*,''
  print*,'===== Test on Error last hyb not equal 1.0, error message is normal ===='  
  stat = vgd_new(d,kind=5,version=1,hyb=hyb(1:size(hyb)-1),rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat.ne.VGD_ERROR)OK=.false.
  print*,''
  print*,'===== hyb not monotone, error message is normal ===='  
  hyb2=hyb
  hyb2(size(hyb)-3)=hyb(size(hyb)-2) 
  stat = vgd_new(d,kind=5,version=1,hyb=hyb2,rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat.ne.VGD_ERROR)OK=.false.

  print*,''
  print*,'===== Test on Error ptop<=0, error message is normal ===='
  ptop=-1.
  stat = vgd_new(d,kind=5,version=1,hyb=hyb,rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat.ne.VGD_ERROR)OK=.false.
  !
  print*,''
  print*,'===== Test on Error ptop lower than first level error message is normal ===='
  ptop=hyb(1)*pref+10.
  stat = vgd_new(d,kind=5,version=1,hyb=hyb,rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat.ne.VGD_ERROR)OK=.false.
  !
  print*,''
  print*,'===== Test complete level set ===='
  ptop=hyb(1)*pref
  w1=ptop
  stat = vgd_new(d,kind=5,version=1,hyb=hyb,rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  if(stat .eq. VGD_ERROR)OK=.false.
  if(.not.(abs(ptop-w1)/w1<my_epsilon))then
     OK=.false.
     print*,'Probleme avec ptop 5001'
     print*,ptop,' should equal ',w1
  endif
  file='data/data_constructor_gen_5001.txt'
  stat = test_1002_5001(d,file,write_control_L)
  if(stat.eq.VGD_ERROR)OK=.false.
  print*,''
  print*,'===== Test incomplete level set ===='
  stat = vgd_new(d,kind=5,version=1,hyb=hyb(10:28),rcoef1=rcoef1,ptop_8=ptop,pref_8=pref)
  stat = vgd_get(d,key='CB_M - vertical B coefficient',value=b_m_8)
  stat = vgd_get(d,key='CA_M - vertical A coefficient',value=a_m_8)
  w1=8.8046513497829437E-002
  if(.not.(abs(b_m_8(1)-w1)/w1<my_epsilon))then
     OK=.false.
     print*,'Probleme avec B 5001'
     print*,b_m_8(1),' should equal ',w1
  endif
  w1=11257.31885433197
  if(.not.(abs(a_m_8(1)-w1)/w1<my_epsilon))then
     OK=.false.
     print*,'Probleme avec A 5001'
     print*,a_m_8(1),' should equal ',w1
  endif

  stat = vgd_free(d)
  deallocate(b_m_8,a_m_8)

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer (5001) value')

end program constructor
!==============================================================================
!===============================================================================
!===============================================================================
integer function test_1002_5001(F_d,F_file,F_write_control_L) result(istat)
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
   integer :: nl_m,nl_t,k,nk,stat,kind,vers
   real :: my_real = 0

   istat=VGD_OK

   nullify(vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8,vipm,vipt,work_i)

   stat = vgd_get(F_d,key='KIND - vertical coordinate ip1 kind' ,value=kind)
   stat = vgd_get(F_d,key='VERS - vertical coordinate version'  ,value=vers)   
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
   
   print*,'TESTING ',kind*1000+vers

   if(F_write_control_L)then
     open(unit=10,file=F_file)
     write(10,*)a_m_8
     write(10,*)b_m_8
     write(10,*)vipm
     write(10,*)vcdm
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
  print*,'Reading A'
  read(10,*)work_8
  do k=1,nk
     
     if(a_m_8(k).eq.0.)then
        if(abs(work_8(k)-a_m_8(k))>1000.*epsilon(my_real))then
           istat=VGD_ERROR
           print*,'Probleme avec A, pas dans les limites tollerees'
           print*,work_8(k),'vs',a_m_8(k)
        endif
     else
        if(abs(work_8(k)-a_m_8(k))/a_m_8(k)>1000.*epsilon(my_real))then
           istat=VGD_ERROR
           print*,'Probleme avec A, pas dans les limites tollerees'
           print*,work_8(k),'vs',a_m_8(k)
        endif
     endif
  enddo
  if(any(a_m_8.ne.a_t_8))then
     istat=VGD_ERROR
     print*,'Probleme avec A T, dot egaler A M :'
     print*,a_t_8
     print*,a_m_8
  endif
 
  ! Check B
  print*,'Reading B'
  read(10,*)work_8
  do k=1,nk
     if(b_m_8(k).eq.0.)then
         if(abs(work_8(k)-b_m_8(k))>10.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B, pas dans les limites tollerees'
           print*,work_8(k),'vs',b_m_8(k)
        endif
     else
        if(abs(work_8(k)-b_m_8(k))/b_m_8(k)>10.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B, pas dans les limites tollerees'
           print*,work_8(k),'vs',b_m_8(k)
        endif
     endif
  enddo
  if(any(b_m_8.ne.b_t_8))then
     istat=VGD_ERROR
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
        istat=VGD_ERROR
        print*,'Probleme avec IP:'
        print*,work_i(k),'vs',vipm(k)
     endif
  enddo
  if(any(vipm.ne.vipt))then
     istat=VGD_ERROR
     print*,'Probleme avec ip1 A , doit egaler ip1 M :'
     print*,vipt
     print*,vipm
  endif

  ! Check vcdm
  allocate(work(size(vcdm)))
  print*,'Reading VCDM'
  read(10,*)work
  do k=1,nk
     if(vcdm(k).eq.0.)then
        if(abs(work(k)-vcdm(k))>10.*epsilon(vcdm(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs',vcdm(k)
        endif
     else
        if(abs(work(k)-vcdm(k))/vcdm(k)>10.*epsilon(vcdm(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs',vcdm(k)
        endif
     endif
  enddo  
  if(any(vcdm.ne.vcdt))then
     istat=VGD_ERROR
     print*,'Probleme avec VCDT, doit egaler VCDM'
     print*,vcdt
     print*,vcdm
  endif

  deallocate(vcdm,vcdt,work,b_m_8,a_m_8,b_t_8,a_t_8,work_8,vipm,vipt,work_i)

end function test_1002_5001
