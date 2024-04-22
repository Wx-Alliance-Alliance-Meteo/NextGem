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

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write
  use Unit_Testing, only: ut_report
  

  implicit none

  type(vgrid_descriptor) :: d  
  integer, parameter :: luo=10,nk=3
  integer :: stat,ip1,k
  integer :: fnom,fstouv,fstfrm
  real :: pp
  real(kind=8), dimension(:), pointer :: a_m_8,b_m_8
  real(kind=8), dimension(:), pointer :: a_t_8, b_t_8
  integer, dimension(:), pointer :: ip1_m,ip1_t
  character(len=10) :: blk_S

  nullify(a_m_8,b_m_8,a_t_8,b_t_8,ip1_m,ip1_t)
  
  stat=fnom(luo,"data_out/dm_5002","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(luo,'RND')

  allocate(ip1_m(nk+1),a_m_8(nk+1),b_m_8(nk+1),ip1_t(nk+2),a_t_8(nk+2),b_t_8(nk+2))

  ! This was obtaines from the listing of a model run
!0:  prgenab Ver_a_8%m(k),Ver_b_8%m(k),k=1,G_nk+1
!0:             1    9.210340386877345        9.6040478224904788E-003
!0:             2    10.59663474799723        0.4828198537713364     
!0:             3    11.28978192855718        0.9280402876755142     
!0:             4    11.51292546497023         1.000000000000000     

!0:  prgenab Ver_a_8%t(k),Ver_b_8%t(k),k=1,G_nk+2
!0:             1    7.873515992234325        4.8020239112452394E-003
!0:             2    9.903487567437288        0.2462119507969134     
!0:             3    10.94320833827721        0.7054300707234253     
!0:             4    11.40135369676370        0.9640201438377572     
!0:             5    11.51292546497023         1.000000000000000     

  a_m_8=(/9.210340386877345,10.59663474799723,11.28978192855718,11.51292546497023/)
  b_m_8=(/9.6040478224904788E-003,0.4828198537713364,0.9280402876755142,1.000000000000000/)
  
  a_t_8=(/7.873515992234325,9.903487567437288,10.94320833827721,11.40135369676370,11.51292546497023/)
  b_t_8=(/4.8020239112452394E-003,0.2462119507969134,0.7054300707234253,0.9640201438377572,1.000000000000000/)

  call convip(ip1,pp,5,0,blk_S,.false.)

  do k=1,nk
     pp=exp(a_m_8(k))/100000.
     call convip(ip1_m(k),pp,5,1,blk_S,.false.)
     !print*,pp,ip1_m(k)
  enddo
  ip1_m(nk+1)=93423264
  ip1_t(nk+2)=93423264
  ! The following is not the exacte code pour good for the exemple
  do k=1,nk+1
     pp=exp(a_t_8(k))/100000.
     call convip(ip1_t(k),pp,5,1,blk_S,.false.)
     !print*,pp,ip1_t(k)
  enddo
  
  stat=vgd_new(&
       d,                  &
       kind    = 5,        &
       version = 2,        &
       nk      = nk,       &
       ip1    = 123,       &
       ip2    = 456,       &
       ptop_8  = 690.d0,   &
       pref_8  = 80000.d0, &
       rcoef1  = 1.0,      &
       rcoef2  = 15.0,     &
       a_m_8   = a_m_8,    &
       b_m_8   = b_m_8,    &
       a_t_8   = a_t_8,    &
       b_t_8   = b_t_8,    &
       ip1_m   = ip1_m,    &
       ip1_t   = ip1_t)

  ! Construct a new set of 3D coordinate descriptors
  stat = vgd_write(d,luo,format='fst')

  call ut_report(stat,'Grid_Descriptors, vgd_new, vgd_write')

  stat=fstfrm(luo)

end program constructor
