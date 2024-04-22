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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_dpidpis,VGD_OK
  use Unit_Testing, only: ut_report
  

  implicit none

  integer :: stat,lu=0,fnom,fstouv,fstfrm,k
  integer, parameter :: i0=1,j0=1
  integer, dimension(:), pointer :: ip1_list
  real, dimension(:,:,:), pointer :: dpidpis_cube
  real(kind=8), dimension(:,:,:), pointer :: dpidpis_cube_8
  real :: eps=1.e-6
  type(vgrid_descriptor) :: d
  logical :: ok
  real(kind=8), dimension(:), pointer :: coef_b
  
  nullify(ip1_list,dpidpis_cube,dpidpis_cube_8,coef_b)

  stat=fnom(lu,"data/dm_5001_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Get dpidpis
  stat = vgd_new(d,unit=lu,format="fst")
  if(stat/=VGD_OK)then
     print*,'ERROR with vgd_new'
     error stop 1
  endif  

  stat = vgd_get(d,key='CB_M - vertical B coefficient (m)',value=coef_b)
  stat = vgd_get(d,key='VIPM - level ip1 list (m)'        ,value=ip1_list)

  !===============
  ! Real interface
  
  stat = vgd_dpidpis(d,ip1_list=ip1_list,dpidpis=dpidpis_cube)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_dpidpis'
     stat=fstfrm(lu)
     error stop 1
  endif

  print*,'size(dpidpis_cube)',size(dpidpis_cube)

  OK=.true.
  do k=1,size(coef_b)
     print*,coef_b(k),dpidpis_cube(i0,j0,k)
     if(coef_b(k).le.epsilon(coef_b(k)))then
        if(abs(coef_b(k)-dpidpis_cube(i0,j0,k))>epsilon(coef_b(k)))then
           print*,'OUPS'
           OK=.false.
           exit
        endif
     else
        if(abs(dpidpis_cube(i0,j0,k)- coef_b(k))/coef_b(k)>eps)then
           print*,'OUPS'
           OK=.false.
           exit
        endif
     endif
  enddo

  !=================
  ! Real*8 interface
  
  stat = vgd_dpidpis(d,ip1_list=ip1_list,dpidpis=dpidpis_cube_8)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_dpidpis real(kind=8)'
     stat=fstfrm(lu)
     error stop 1
  endif

  print*,'size(dpidpis_cube_8)',size(dpidpis_cube_8)

  OK=.true.
  do k=1,size(coef_b)
     print*,coef_b(k),dpidpis_cube_8(i0,j0,k)
     if(coef_b(k).le.epsilon(coef_b(k)))then
        if(abs(coef_b(k)-dpidpis_cube_8(i0,j0,k))>epsilon(coef_b(k)))then
           print*,'dpidpis for real(kind=8) do not validate'
           OK=.false.
           exit
        endif
     else
        if(abs(dpidpis_cube_8(i0,j0,k)- coef_b(k))/coef_b(k)>eps)then
           print*,'dpidpis for real(kind=8) do not validate'
           OK=.false.
           exit
        endif
     endif
  enddo

  call ut_report(ok,message='Grid_Descriptors::vgd_dpidpis level calculation status')

  stat=fstfrm(lu)

end program tests
