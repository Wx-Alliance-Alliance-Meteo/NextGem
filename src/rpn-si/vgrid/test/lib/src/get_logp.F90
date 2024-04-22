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
program get_logp
  use Vgrid_Descriptors, only: Vgrid_descriptor,Vgd_free,Vgd_new,Vgd_get,Vgd_print,VGD_OK,vgd_putopt
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10
  integer :: stat,i
  integer :: fnom,fstouv,fstfrm,fclos
  logical :: ok,logp_L

  integer, parameter :: nfiles=15
  character(len=200), dimension(nfiles) :: files=(/&
       "data/dm_1001_from_model_run           ",&
       "data/dm_1002_from_model_run           ",&
       "data/dm_2001_from_editfst             ",&
       "data/dm_4001_from_model_run           ",&
       "data/dm_5001_from_model_run           ",&
       "data/dm_5002_from_model_run           ",&
       "data/dm_5003_from_model_run           ",&
       "data/dm_5004_from_model_run           ",&
       "data/dm_5005_from_model_run           ",&
       "data/dm_5100_from_model_run           ",&
       "data/dm_5999_from_model_run           ",&
       "data/dm_21001_from_model_run_SLEVE    ",&
       "data/dm_21001_from_model_run_NON_SLEVE",&
       "data/dm_21002_from_model_run_SLEVE    ",&
       "data/dm_21002_from_model_run_NON_SLEVE"&
       /)

  logical, dimension(nfiles) :: is_in_log =(/&
       .false., &
       .false., &
       .false., &
       .false., &
       .false., &
       .true., &
       .true., &
       .true., &
       .true., &
       .true., &
       .false., &
       .false., &
       .false., &
       .false., &
       .false. &
       /)

  stat = vgd_putopt("ALLOW_SIGMA",.true.)

  ok=.true.

  do i=1,nfiles
     print*,'======================================================'
     print*,'file = ',trim(files(i))
     stat=fnom(lu+i,files(i),"RND",0)
     if(stat.lt.0)then
        print*,'ERROR with fnom on file ',trim(files(i))
        error stop 1
     endif
     stat=fstouv(lu+i,'RND')
     if(stat.le.0)then
        print*,'No record in RPN file ',trim(files(i))
        error stop 1
     endif     
     stat = vgd_new(d,unit=lu+i,format="fst",ip1=-1,ip2=-1)
     if(stat.ne.VGD_OK)then
        print*,'ERROR: problem with vgd_new on file ',trim(files(i))
        stat=fstfrm(lu+i)
        error stop 1
     endif
     stat = vgd_get(d,key='LOGP',value=logp_L)
     if(stat /= VGD_OK )then
        print*,'Error with vgd_get on LOGP on file ',trim(files(i))
        ok=.false.
     endif
     if(logp_L)then
        if(.not.is_in_log(i))then
           print*,'ERROR file ',trim(files(i))," is in log but vgd_get on 'LOGP' returned .false."
           ok=.false.
        endif
     else
        if(is_in_log(i))then
           print*,'ERROR file ',trim(files(i))," is not in log but vgd_get on 'LOGP' returned .true."
           ok=.false.
        endif
     endif
     stat=fstfrm(lu+i)  
     stat=fclos(lu+i)

  end do

  call ut_report(ok,'Grid_Descriptors, vgd_get, LOGP 5002')  

end program get_logp
