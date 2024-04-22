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
program get_levels_from_keys

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels,VGD_ERROR
  use Unit_Testing, only: ut_report

  implicit none

  type(vgrid_descriptor) :: d
  integer, parameter :: lu=10, nmax=1000, lutxt=69
  integer, dimension(nmax) :: liste
  integer :: stat,infon,ier,i,j,k
  integer :: fnom,fstouv,fstfrm,fstinf,fstluk,fstinl,fstprm
  real, dimension(:,:,:), pointer :: pres
  real, dimension(:,:), pointer :: ff
  ! Variable for fstprm
  integer ::dateo, datev, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
       ig2, ig3, ig4, ip1, ip2, ip3, key, lng, nbits,&
       ni,  nj, nk, npas, swa, ubc
  character(len=12) :: etiket
  character(len=4)  :: nomvar
  character(len=2)  :: typvar
  character(len=1)  :: grtyp

  nullify(pres,ff)

  stat=fnom(lu,"data/dm_5002_from_model_run","RND",0)
  if(stat.lt.0)then
     print*,'ERROR with fnom'
     error stop 1
  endif
  stat=fstouv(lu,'RND')
  if(stat.le.0)then
     print*,'No record in RPN file'
     error stop 1
  endif

  ! Get vertical grid descriptor
  stat = vgd_new(d,unit=lu,format="fst")
  if(stat < 0 )then
     print*,'Problem getting vertical grid descriptor'
     error stop 1
  endif

  ier = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',liste,infon,nmax)

  if(infon == 0 )then
     print*,'pas de record de TT'
     error stop 1
  endif

  stat = vgd_levels(d,unit=lu,fstkeys=liste(1:infon),levels=pres)
  if(stat < 0 )then
     print*,'Problem with vgd_levels'
     error stop 1
  endif
  
  pres=pres/100.

  allocate(ff(ni,nj))

  do k=1,infon

     ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
     ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
     dltf,ubc,extra1,extra2,extra3)
     if(ier.lt.0)then
        print*,'Problem with fstprm on TT, ip1=',ip1
        error stop 1
     endif
     call incdatr(datev,dateo,deet*npas/3600.d0)
     key=fstinf(lu,ni,nj,nk,datev,' ',ip1,ip2,-1,typvar,'PX')     
     ier = fstluk(ff,key,ni,nj,nk)
      if(ier.lt.0)then
        print*,'Problem with fstinf on PX, ip1=',ip1
        error stop 1
     endif

     do j=1,nj
        do i=1,ni
           if(abs((ff(i,j)-pres(i,j,k))/ff(i,j))>1.0e-6)then
              stat=VGD_ERROR
              print*,'Difference in pressure is too large at'
              print*,'i,j,k,ff(i,j),pres(i,j,k)',i,j,k,ff(i,j),pres(i,j,k)
              exit
           endif
        enddo
     enddo

  enddo

  call ut_report(stat,'Grid_Descriptors, vgd_new')

  stat=fstfrm(lu)

end program get_levels_from_keys
