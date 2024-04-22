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
  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels,VGD_OK
  use Unit_Testing, only: ut_report

  implicit none

  integer :: stat,lu=0,ni,nj,fnom,fstouv,fstfrm,fstlir,fstinf,k,lutxt=69,ip1,ip2,nk
  integer, parameter :: i0=20,j0=10
  integer, dimension(:), pointer :: ip1_list
  real, dimension(:), pointer :: pres_profil
  real :: epsilon=0.01
  real, dimension(:,:), pointer :: p0,px
  real, dimension(:,:,:), pointer :: lev
  real :: local_pres
  type(vgrid_descriptor) :: d
  logical :: ok

  nullify(ip1_list,pres_profil,p0,px,lev)

  stat=fnom(lu,"data/gg_5001_from_model_run_plus_toc","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.lt.0)then
     print*,'No record in RPN file'
     error stop 1
  endif
  open(unit=lutxt,file='data/gg_5001_ips.txt',status='OLD')
  read(lutxt,*) ip1,ip2
  close(lutxt)

  ! Get physical levelling information
  stat = vgd_new(d,unit=lu,format="fst",ip1=ip1,ip2=ip2)
  
  stat = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',"UU")

  allocate(p0(ni,nj),ip1_list(5),px(ni,nj))

  ip1_list=(/95332840,95345840,95356840,95366840,93423264/)

  stat = fstlir(p0,lu,ni,nj,nk,-1,'',-1,-1,-1,'','P0')
  p0 = p0*100. !mb to Pa
  
  local_pres=p0(i0,j0)  

  stat = vgd_levels(d,sfc_field=local_pres,ip1_list=ip1_list,levels=pres_profil)
  if(stat.ne.VGD_OK)then
     print*,'ERROR: problem with vgd_levels'
     stat=fstfrm(lu)
     error stop 1
  endif

  pres_profil=pres_profil*.01
  print*,'size(pres_profil)',size(pres_profil)

  OK=.true.
  do k=1,nk
     stat = fstlir(px,lu,ni,nj,nk,-1,'',ip1_list(k),-1,-1,'','PX')
     print*,px(i0,j0),pres_profil(k)
     if(abs(pres_profil(k)-px(i0,j0))>epsilon)then
        print*,'OUPS'
        OK=.false.
     endif
  enddo


  call ut_report(ok,message='Grid_Descriptors::vgd_levels level calculation status')

  stat=fstfrm(lu)

end program tests
