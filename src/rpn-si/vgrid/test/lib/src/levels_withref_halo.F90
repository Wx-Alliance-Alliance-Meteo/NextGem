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
program level_withref_5001

  use vGrid_Descriptors, only: vgd_putopt,VGD_OK,VGD_ERROR
  use Unit_Testing, only: ut_report
  
  implicit none
  integer :: stat,ier,chek_levels_withref

  stat=VGD_OK

  ier=vgd_putopt('ALLOW_SIGMA',.true.)
  ier=chek_levels_withref('data/dm_1001_from_model_run','')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=vgd_putopt('ALLOW_SIGMA',.false.)

  ier=chek_levels_withref('data/dm_1002_from_model_run','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_1002_from_model_run_above_0.5','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/2001_from_model_run','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_5001_from_model_run','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_5002_from_model_run','data/dm_5002_ips.txt')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  
  ier=chek_levels_withref('data/dm_5003_from_model_run','data/dm_5003_ips.txt')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  
  ier=chek_levels_withref('data/dm_5002_from_model_run_ig4_ip1_link','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_21001_from_model_run_NON_SLEVE','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_21001_from_model_run_SLEVE','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_21002_from_model_run_NON_SLEVE','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=chek_levels_withref('data/dm_21002_from_model_run_SLEVE','')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  call ut_report(stat,'Grid_Descriptors, vgd_new')

end program level_withref_5001
!====================================================================
!====================================================================
!====================================================================

integer function chek_levels_withref(F_fst,F_ips) result(status)

   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_putopt,vgd_levels,vgd_get,vgd_free,VGD_ERROR,VGD_OK, VGD_NO_REF_NOMVAR
   

   implicit none  

   character(len=*) :: F_fst,F_ips

   ! Local variables
   integer :: minx, maxx, miny, maxy, G_nk
   integer, save :: lu=10   
   integer :: fnom,fstouv,fstfrm,lutxt=69,kind
   type(vgrid_descriptor) :: vgd
   integer, parameter :: nmax=1000
   integer, dimension(nmax) :: liste
   integer :: ier,fstinl,fstprm,fstinf,fstluk,infon,i,j,k
   real, dimension(:,:,:), pointer :: pres
   real(kind=8), dimension(:,:,:), pointer :: pres_8
   real, dimension(:,:), pointer :: p0, p0ls,px
   real(kind=8), dimension(:,:), pointer :: p0_8, p0ls_8
   real :: epsilon=5.0e-4,pppp, factor
   integer, dimension(:), pointer :: ip1s
   ! Variable for fstprm, sorry...
   integer ::dateo, datev, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
        ig2, ig3, ig4, ip1, ip2, ip3, key, lng, nbits,&
        ni,  nj, nk, npas, swa, ubc
   character(len=12) :: etiket
   character(len=4)  :: nomvar, rfld, nomvar_metric
   character(len=2)  :: typvar
   character(len=1)  :: grtyp, dummy_S
   logical :: sfc_field_ls_L

   status=VGD_ERROR   

   nullify(pres, pres_8, p0, p0ls, px, p0_8, p0ls_8, ip1s)

   ier=vgd_putopt('ALLOW_RESHAPE',.true.)

   lu=lu+1

   print*,'Tests '//trim(F_fst)

   ier=fnom(lu,F_fst,"RND",0)
   if(ier.lt.0)then
      print*,'(Test) ERROR with fnom on ',trim(F_fst)
      return
   endif
   ier=fstouv(lu,'RND')
   if(ier.le.0)then
      print*,'(Test) No record in RPN file ',trim(F_fst),ier
      return
   endif
   if(trim(F_ips).eq.'')then
      ip1=-1; ip2=-1     
   else
      open(unit=lutxt,file=F_ips,status='OLD')
      read(lutxt,*) ip1,ip2
      close(lutxt)
   endif

   ! Get vertical grid descriptor
   ier = vgd_new(vgd,unit=lu,format="fst",ip1=ip1,ip2=ip2)
   if(ier == VGD_ERROR )then
      print*,'(Test) Problem getting vertical grid descriptor'
      return
   endif

   ier = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',liste,infon,nmax)
   if(infon == 0 )then
      print*,'(Test) pas de record de TT'      
      return
   endif
   allocate(ip1s(infon))
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      ip1s(k)=ip1
   end do
   
   if(.true.)then
      minx=-1; maxx=ni-2
      miny=-1; maxy=nj-2      
   else
      minx=1; maxx=ni
      miny=1; maxy=nj
   endif
   G_nk=infon
      
   
   sfc_field_ls_L=.false.
   call convip_plus(ip1,pppp,kind,-1,dummy_S,.false.)
   if(kind .eq. 2)then
      ! pressure levels do not have reference surface fields
      allocate(p0(minx:maxx,miny:maxy),px(minx:maxx,miny:maxy),p0ls(minx:maxx,miny:maxy),pres(minx:maxx,miny:maxy,G_nk))
      factor = 100.
      p0=0.
      p0ls=0.
   else
      ier = vgd_get(vgd,'RFLD',rfld)
      if(ier == VGD_ERROR)then
         print*,'(Test) Problem with vgd_get "RFLD" for file ',F_fst
         return
      endif
      key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfld)
      allocate(p0(minx:maxx,miny:maxy),px(minx:maxx,miny:maxy),p0ls(minx:maxx,miny:maxy),pres(minx:maxx,miny:maxy,G_nk))
      ier = fstluk(p0,key,ni,nj,nk)
      if(ier.lt.0)then
         print*,'(Test) Problem with fstluk on ',rfld
         return
      endif
      if(trim(rfld) == "P0")then
         nomvar_metric = "PX"
         factor = 100.
         p0=p0*100.
      else
         nomvar_metric = "GZ"
         factor = 10.
      endif
      ier = vgd_get(vgd,'RFLS',rfld,quiet=.true.)
      if(rfld /= VGD_NO_REF_NOMVAR)sfc_field_ls_L=.true.
      if(sfc_field_ls_L)then
         key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfld)
         ier = fstluk(p0ls,key,ni,nj,nk)
         if(ier.lt.0)then
            print*,'(Test) Problem with fstluk on ',rfld
            return
         endif
      else
         p0ls = 0.
      endif
   endif

   !print*,'p0(minx,miny),pres(minx,miny,G_nk)',p0(minx,miny),pres(minx,miny,G_nk)
   !print*,'shape(p0)',shape(p0)
   !print*,'shape(pres)',shape(pres)
   !print*,'lbound(p0,1),ubound(p0,1)',lbound(p0,1),ubound(p0,1)
   !print*,'lbound(pres,1),ubound(pres,1)',lbound(pres,1),ubound(pres,1)
      
   allocate(p0_8(minx:maxx,miny:maxy),p0ls_8(minx:maxx,miny:maxy),pres_8(minx:maxx,miny:maxy,G_nk))
   p0_8=p0
   p0ls_8=p0ls
   ! Test 32 bits interface
   if(sfc_field_ls_L)then
      ier = vgd_levels(vgd,ip1s,pres,p0,sfc_field_ls=p0ls)
   else
      ier = vgd_levels(vgd,ip1s,pres,p0)
   endif
   if(ier == VGD_ERROR )then
      print*,'(Test) Problem with vgd_levels 32 bits'
      return
   endif

   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      if(ier.lt.0)then
         print*,'(Test) Problem with fstprm on TT, ip1=',ip1
         return
      endif
      if(kind .eq. 2)then
         call convip_plus(ip1,pppp,kind,-1,dummy_S,.false.)
         px = pppp
      else
         ! Surface pressure must be equal to P0 for surface levels
         call convip_plus(ip1,pppp,kind,-1,dummy_S,.false.)
         if(abs(pppp-1.).lt.epsilon)then
            do j=miny,maxy
               do i=minx,maxx
                  if(pres(i,j,k).ne.p0(i,j))then
                     print*,'(Test) Surface pressure must be exacly equal to p0'
                     print*,'i,j,k,pres(i,j,k),p0(i,j)',i,j,k,pres(i,j,k),p0(i,j)
                     return
                  endif
               enddo
            enddo
         endif
         call incdatr(datev,dateo,deet*npas/3600.d0)
         key=fstinf(lu,ni,nj,nk,datev,' ',ip1,ip2,-1,typvar,nomvar_metric)
         if(key < 0 )then
            print*,'Cannot find nomvar ',nomvar_metric,' for ip1',ip1,' in file ',F_fst
            return
         endif
         ier = fstluk(px,key,ni,nj,nk)
      endif
      px = px * factor
      if(ier.lt.0)then
         print*,'(Test) Problem with fstinf on PX, ip1=',ip1
         return
      endif
      do j=miny,maxy
         do i=minx,maxx
            if(abs((px(i,j)-pres(i,j,k))/px(i,j))>epsilon)then
               print*,'(Test) 32 bits: Difference in pressure is too large at'
               print*,'i,j,k,px(i,j),pres(i,j,k)',i,j,k,px(i,j),pres(i,j,k)
               return
            endif
         enddo
      enddo
   enddo

   ! Test 64 bits interface
   if(sfc_field_ls_L)then
      ier = vgd_levels(vgd,ip1s,pres_8,p0_8,sfc_field_ls=p0ls_8)
   else
      ier = vgd_levels(vgd,ip1s,pres_8,p0_8)
   endif
   if(ier == VGD_ERROR )then
      print*,'(Test) Problem with vgd_levels 64 bits'
      return
   endif 
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      if(ier.lt.0)then
         print*,'Problem with fstprm on TT, ip1=',ip1
         return
      endif
      if(kind .eq. 2)then
         call convip_plus(ip1,pppp,kind,-1,dummy_S,.false.)
         px = pppp
      else
         ! Surface pressure must be equal to P0 for surface levels
         call convip_plus(ip1,pppp,kind,-1,dummy_S,.false.)
         if(abs(pppp-1.).lt.epsilon)then
            do j=miny,maxy
               do i=minx,maxx
                  if(pres_8(i,j,k).ne.p0_8(i,j))then
                     print*,'Surface pressure must be exacly equal to p0_8'
                     print*,'File:',F_fst
                     print*,'i,j,k,pres_8(i,j,k),p0_8(i,j)',i,j,k,pres_8(i,j,k),p0_8(i,j)
                     return
                  endif
               enddo
            enddo
         endif
         call incdatr(datev,dateo,deet*npas/3600.d0)
         key=fstinf(lu,ni,nj,nk,datev,' ',ip1,ip2,-1,typvar,nomvar_metric)
         if(key < 0 )then
            print*,'Cannot find nomvar ',nomvar_metric,' for ip1',ip1,' in file ',F_fst
            return
         endif
         ier = fstluk(px,key,ni,nj,nk)
      endif
      px = px * factor
      if(ier.lt.0)then
         print*,'(Test) Problem with fstinf on PX, ip1=',ip1
         return
      endif
      do j=miny,maxy
         do i=minx,maxx
            if(abs((px(i,j)-pres_8(i,j,k))/px(i,j))>epsilon)then
               print*,'(Test) 64 bits: Difference in pressure is too large at'
               print*,'i,j,k,px(i,j),pres_8(i,j,k)',i,j,k,px(i,j),pres_8(i,j,k)
               return
            endif
         enddo
      enddo
   enddo
   deallocate(ip1s,px,p0,p0_8,pres,pres_8)
   
   ier=vgd_free(vgd)

   ier=fstfrm(lu)
   status=VGD_OK   
   
end function chek_levels_withref
