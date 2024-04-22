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

  ! Revision : Andre Plante test on B instead of A since A not sensitive to rcoefs

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_print,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer :: test_hgts, stat
  real, dimension(80) :: hgts= &
(/ 6.62078E+04, 6.29395E+04, 5.95406E+04, 5.60645E+04, 5.26367E+04, 4.93408E+04, 4.63244E+04, 4.35656E+04, 4.10827E+04, 3.88781E+04,&
 3.68699E+04, 3.50538E+04, 3.33957E+04, 3.18990E+04, 3.05104E+04, 2.92230E+04, 2.80314E+04, 2.69268E+04, 2.59083E+04, 2.49636E+04,&
 2.41007E+04, 2.33054E+04, 2.25525E+04, 2.18511E+04, 2.12170E+04, 2.06171E+04, 2.00673E+04, 1.95593E+04, 1.90950E+04, 1.86683E+04,&
 1.82816E+04, 1.79227E+04, 1.76021E+04, 1.73039E+04, 1.70024E+04, 1.66971E+04, 1.63648E+04, 1.60072E+04, 1.56247E+04, 1.52095E+04,&
 1.47589E+04, 1.42767E+04, 1.37586E+04, 1.31993E+04, 1.26335E+04, 1.20691E+04, 1.15349E+04, 1.10244E+04, 1.05204E+04, 1.00221E+04,&
 9.51942E+03, 9.02303E+03, 8.53161E+03, 8.05913E+03, 7.59641E+03, 7.13710E+03, 6.68393E+03, 6.23286E+03, 5.78976E+03, 5.35941E+03,&
 4.93475E+03, 4.51786E+03, 4.10324E+03, 3.69539E+03, 3.29876E+03, 2.91503E+03, 2.56689E+03, 2.24003E+03, 1.95540E+03, 1.68822E+03,&
 1.43746E+03, 1.22341E+03, 1.02339E+03, 8.36754E+02, 6.73131E+02, 5.21595E+02, 3.81767E+02, 2.53260E+02, 1.45511E+02, 4.83247E+01/)


  real :: rcoef1=0.,rcoef2=4.,rcoef3=0.,rcoef4=200.
  
  logical :: OK=.true.
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file

  ! Construct a new set of vertical coordinate descriptors 21101 Gal-Chen
  if( vgd_new(d,kind=21,version=1,hyb=hgts,rcoef1=rcoef1,rcoef2=rcoef2,rcoef3=rcoef3,rcoef4=rcoef4,dhm=10.0,dht=1.5) == &
       VGD_ERROR )OK=.false.
  if( vgd_print(d,convip_L=.true.) == VGD_ERROR )then
     print*,'ERROR'
     error stop 1
  endif
  if( vgd_print(21101) == VGD_ERROR )then
     print*,'ERROR'
     error stop 1
  endif

  file='data/data_constructor_gen_21001_SLEVE.txt'
  stat = test_hgts(d,file,write_control_L)
  if(stat.eq.VGD_ERROR)OK=.false.

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer SLEVE like coordinate')
end program constructor
!==============================================================================
!===============================================================================
!===============================================================================
integer function test_hgts(F_d,F_file,F_write_control_L) result(istat)
   !
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_new,vgd_free,operator(==),VGD_ERROR,VGD_OK
   
   implicit none
   !
   type (vgrid_descriptor) :: F_d
   logical :: F_write_control_L
   character (len=256) :: F_file
   !
   ! Local variable
   !
   real, dimension(:), pointer :: vcdm,vcdt,work
   real(kind=8), dimension(:), pointer :: b_m_8,a_m_8,c_m_8,a_t_8,b_t_8,c_t_8,work_8
   integer, dimension(:), pointer :: vipm,vipt,work_i
   integer :: nl_m,nl_t,k,nk,kind,vers,ip1,my_ip1
   real(kind=8), dimension(:,:,:), pointer :: table
   type (vgrid_descriptor) :: vgrid_rebuilt

   nullify(vcdm,vcdt,work,a_m_8,b_m_8,c_m_8,a_t_8,b_t_8,c_t_8,work_8,vipm,vipt,work_i,table)

   istat=VGD_ERROR

   flush(6)

   if(vgd_get(F_d,'VTBL',table) == VGD_ERROR)return
   if(vgd_new(vgrid_rebuilt,table) ==  VGD_ERROR)return
   if (vgrid_rebuilt == F_d) then
      ! Ok do noting
   else
      print*,'ERROR: rebuilding table'
      return
   endif
   if(vgd_free(vgrid_rebuilt) ==  VGD_ERROR)return
   deallocate(table)
   
   if( vgd_get(F_d,key='KIND - vertical coordinate ip1 kind' ,value=kind)   == VGD_ERROR) return
   if( vgd_get(F_d,key='VERS - vertical coordinate version'  ,value=vers)   == VGD_ERROR) return
   if( vgd_get(F_d,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CA_T - vertical A coefficient (t)'   ,value=a_t_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CB_M - vertical B coefficient (m)'   ,value=b_m_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CB_T - vertical B coefficient (t)'   ,value=b_t_8)  == VGD_ERROR) return
   if( vgd_get(F_d,key='CC_M - vertical C coefficient (m)'   ,value=c_m_8)  == VGD_ERROR) return
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

  istat=VGD_OK

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
     if(a_m_8(k).eq.0.)then
        if(abs(work_8(k)-a_m_8(k))>100.*epsilon(a_m_8(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,work_8(k),'vs'
           print*,a_m_8(k)
        endif
     else
        if(abs(work_8(k)-a_m_8(k))/a_m_8(k)>100.*epsilon(a_m_8(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec A M, pas dans les limites tollerees'
           print*,work_8(k),'vs'
           print*,a_m_8(k)
        endif
     endif
  enddo

  print*,'Reading A T'
  read(10,*)work_8
  do k=1,nk
     if(a_t_8(k).eq.0.)then
        if(abs(work_8(k)-a_t_8(k))>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec A T, pas dans les limites tollerees'
           print*,work_8(k),'vs'
           print*,a_t_8(k)
        endif
     else
        if(abs(work_8(k)-a_t_8(k))/a_t_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec A T, pas dans les limites tollerees'
           print*,work_8(k),'vs'
           print*,a_t_8(k)
        endif
     endif
  enddo

  ! Check B
  print*,'Reading B M'
  read(10,*)work_8
  do k=1,nk
     if(b_m_8(k).eq.0.)then
        if(work_8(k).ne.0.)then
           istat=VGD_ERROR
           print*,'Probleme avec B M, pas egal a zero ',work_8(k)
        endif
     else
        if(abs(work_8(k)-b_m_8(k))/b_m_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B M, pas dans les limites tollerees'
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
           print*,'Probleme avec B T, pas egal a zero ',work_8(k)
        endif
     else
        if(abs(work_8(k)-b_t_8(k))/b_t_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec B T, pas dans les limites tollerees'
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
           print*,'Probleme avec C M, pas egal a zero ',work_8(k)
        endif
     else
        if(abs(work_8(k)-c_m_8(k))/c_m_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec C M, pas dans les limites tollerees'
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
           print*,'Probleme avec C T, pas egal a zero ',work_8(k)
        endif
     else
        if(abs(work_8(k)-c_t_8(k))/c_t_8(k)>100.*epsilon(1.))then
           istat=VGD_ERROR
           print*,'Probleme avec C T, pas dans les limites tollerees'
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
        print*,'Probleme avec IP M:'
        print*,work_i(k),'vs'
        print*,vipm(k)
     endif
  enddo
  print*,'Reading IP T'
  read(10,*)work_i
  do k=1,nk
     if(work_i(k).ne.vipt(k))then
        istat=VGD_ERROR
        print*,'Probleme avec IP T:'
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
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs'
           print*,vcdm(k)
        endif
     else
        if(abs(work(k)-vcdm(k))/vcdm(k)>100.*epsilon(vcdm(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
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
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs'
           print*,vcdt(k)
        endif
     else
        if(abs(work(k)-vcdt(k))/vcdt(k)>100.*epsilon(vcdt(k)))then
           istat=VGD_ERROR
           print*,'Probleme avec vcdt, pas dans les limites tollerees'
           print*,work(k),'vs'
           print*,vcdt(k)
        endif
     endif
  enddo

  print*,'Reading ip1'
  read(10,*)my_ip1
  if(my_ip1.ne.ip1)then
      istat=VGD_ERROR
     print*,'Probleme avec ip1, expected ',ip1,' got',my_ip1     
  endif

  close(10)

  deallocate(vcdm,vcdt,work,a_m_8,b_m_8,a_t_8,b_t_8,work_8,vipm,vipt,work_i)

end function test_hgts
