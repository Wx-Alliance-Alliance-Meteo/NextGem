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
program levels_withref_profile_all

  use vGrid_Descriptors, only: vgd_putopt,VGD_OK,VGD_ERROR
  use Unit_Testing, only: ut_report
  
  implicit none
  integer :: stat,ier,check_levels_withref

  stat=VGD_OK

  ier = vgd_putopt("ALLOW_SIGMA",.true.)
  
  ier=check_levels_withref('data/dm_1001_from_model_run','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier =  vgd_putopt("ALLOW_SIGMA",.false.)

  ier=check_levels_withref('data/dm_1002_from_model_run','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_4001_from_model_run','','GZ')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_5001_from_model_run','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_5002_from_model_run','data/dm_5002_ips.txt','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5002_from_model_run','data/dm_5002_ips.txt','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  
  ier=check_levels_withref('data/dm_5003_from_model_run','data/dm_5003_ips.txt','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5003_from_model_run','data/dm_5003_ips.txt','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  
  ier=check_levels_withref('data/dm_5004_from_model_run','data/dm_5004_ips.txt','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5004_from_model_run','data/dm_5004_ips.txt','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_5005_from_model_run','data/dm_5005_ips.txt','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5005_from_model_run','data/dm_5005_ips.txt','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_5100_from_model_run','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5100_from_model_run','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_5999_from_model_run','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_5999_from_model_run','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_21001_from_model_run_SLEVE','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_21001_from_model_run_SLEVE','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_21001_from_model_run_NON_SLEVE','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_21001_from_model_run_NON_SLEVE','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_21002_from_model_run_SLEVE','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_21002_from_model_run_SLEVE','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  ier=check_levels_withref('data/dm_21002_from_model_run_NON_SLEVE','','TT')
  if(ier==VGD_ERROR)stat=VGD_ERROR
  ier=check_levels_withref('data/dm_21002_from_model_run_NON_SLEVE','','UU')
  if(ier==VGD_ERROR)stat=VGD_ERROR

  call ut_report(stat,'Grid_Descriptors, vgd_new')

end program levels_withref_profile_all
!====================================================================
!====================================================================
!====================================================================

integer function check_levels_withref(F_fst,F_ips,F_var) result(status)

   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_levels,vgd_get,vgd_print,VGD_LEN_RFLD,VGD_NO_REF_NOMVAR,VGD_ERROR,VGD_OK
   
  
   implicit none  

   character(len=*) :: F_fst,F_ips,F_var

   ! Local variables
   integer, save :: lu=10   
   integer :: fnom,fstouv,fstfrm,lutxt=69,kind
   type(vgrid_descriptor) :: vgd
   integer, parameter :: nmax=1000
   integer, dimension(nmax) :: liste
   integer :: ier,fstinl,fstprm,fstinf,fstluk,infon,k,i,j,i0,j0,vcode
   real, dimension(:), pointer :: pres
   real(kind=8), dimension(:), pointer :: pres_8
   real, dimension(:,:), pointer :: p0,p0ls,px
   real :: epsilon=5.0e-6,pppp,p0_point,p0ls_point,fact
   real(kind=8) :: p0_point_8,p0ls_point_8,fact_8
   integer, dimension(:), pointer :: ip1s
   character(len=VGD_LEN_RFLD) :: rfld_S, rfls_S
   ! Variable for fstprm, sorry...
   integer ::dateo, datev, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
        ig2, ig3, ig4, ip1, ip2, ip3, key, lng, nbits,&
        ni,  nj, nk, npas, swa, ubc
   character(len=12) :: etiket
   character(len=4)  :: nomvar,nomvar_metric
   character(len=2)  :: typvar
   character(len=1)  :: grtyp, dummy_S
   logical :: two_refs_L, ok
   
   nullify(pres,pres_8,p0,p0ls,px,ip1s)

   status=VGD_ERROR   

   print*,'================================================================'
   print*,'Testing, level: ',trim(F_var),' on file ',trim(F_fst)

   lu=lu+1
   
   ier=fnom(lu,F_fst,"RND+R/O",0)
   if(ier.lt.0)then
      print*,'ERROR in test with fnom on ',trim(F_fst)
      return
   endif
   ier=fstouv(lu,'RND')
   if(ier.le.0)then
      print*,'ERROR in test: No record in RPN file ',trim(F_fst),ier
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
   ier = vgd_print(vgd)
   if(ier == VGD_ERROR )then
      print*,'ERROR in test: Problem getting vertical grid descriptor'
      print*,'FILE: ',trim(F_fst)
      return
   endif
   
   ier = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',F_var,liste,infon,nmax)
   if(ier>0)then
      print*,'ERROR in test: problem with fstinl on ',F_var
      return
   endif
   if(infon == 0 )then
      print*,'ERROR in test: pas de record de ',F_var
      print*,'FILE: ',trim(F_fst)
      return
   endif
   allocate(ip1s(infon))
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      ip1s(k)=ip1
   end do
   
   ier = vgd_get(vgd, "RFLD", rfld_S, .true.);
   if( rfld_S == VGD_NO_REF_NOMVAR )then
      ier = vgd_get(vgd, "VCOD", vcode);
      if(vcode == 4001)then         
         ! Allocate P0, used only to get horizontale size         
         key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','GZ')
         allocate(p0(ni,nj),px(ni,nj))
         p0=0.; p0_point=0.; i0=1; j0=1
      else
         print*,'ERROR in test: in test reference field, Vcode not supported: ',trim(rfld_S)
         return
      endif
   else
      key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfld_S)
      allocate(p0(ni,nj),px(ni,nj))
      ier = fstluk(p0,key,ni,nj,nk)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstluk on ',rfld_S
         print*,'FILE: ',trim(F_fst)
         return
      endif
      if(trim(rfld_S) == "P0")then
         ! Find lowest pressure point
         p0_point=10000.
         do j=1,nj
            do i=1,ni
               if(p0(i,j)<p0_point)then
                  i0=i; j0=j; p0_point=p0(i,j)
               endif
            end do
         end do
         p0_point=p0(i0,j0)*100.
      elseif(trim(rfld_S) == "ME")then
         ! Find heighest height point
         p0_point=0
         do j=1,nj
            do i=1,ni
               if(p0(i,j)>p0_point)then
                  i0=i; j0=j; p0_point=p0(i,j)
               endif
            end do
         end do
         p0_point=p0(i0,j0)
      else
         print*,'ERROR in test: reference field not supported ',trim(rfld_S)
         return
      endif
      p0_point_8=p0_point
   endif
   
   two_refs_L = .false.
   ier = vgd_get(vgd, "RFLS", rfls_S, .true.);
   if( rfls_S /= VGD_NO_REF_NOMVAR )then      
      two_refs_L = .true.
      key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfls_S)
      allocate(p0ls(ni,nj))
      ier = fstluk(p0ls,key,ni,nj,nk)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstluk on ',rfls_S
         print*,'FILE: ',trim(F_fst)
         return
      endif
      if(trim(rfls_S) == "P0LS")then
         p0ls_point=p0ls(i0,j0)*100.
      else
         p0ls_point=p0ls(i0,j0)
      endif
      p0ls_point_8=p0ls_point
   endif
   ! Test 32 bits interface
   if(two_refs_L)then
      ier = vgd_levels(vgd,ip1s,pres,p0_point,sfc_field_ls=p0ls_point)
   else
      ier = vgd_levels(vgd,ip1s,pres,p0_point)
   endif
   if(ier == VGD_ERROR )then
      print*,'ERROR in test: Problem with vgd_levels 32 bits'
      print*,'FILE: ',trim(F_fst)
      return
   endif
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstprm on ',F_var,', ip1=',ip1
         print*,'FILE: ',trim(F_fst)
         return
      endif
      ! Surface pressure must be equal to P0 for surface levels
      call convip(ip1,pppp,kind,-1,dummy_S,.false.)
      !print*,'ip1, level = ',ip1,pppp
      if(abs(pppp-1.).lt.epsilon.or. ( abs(pppp).lt.epsilon .and. trim(rfld_S) == "ME") )then
         if(pres(k).ne.p0_point)then
            print*,'ERROR in test: 32 bits level at surface must be exacly equal to ',rfld_S
            print*,'k,pres(k),p0_point',k,pres(k),p0_point
            print*,'FILE: ',trim(F_fst)
            return
         endif
      endif
      call incdatr(datev,dateo,deet*npas/3600.d0)
      if(trim(rfld_S) == "P0")then
         nomvar_metric="PX"
         fact=.01
         fact_8=.01
      elseif(trim(rfld_S) == "ME")then
         nomvar_metric="GZ"
         fact=.1
         fact_8=.1d0
      elseif(trim(rfld_S) == trim(VGD_NO_REF_NOMVAR))then
         ier = vgd_get(vgd, "VCOD", vcode);
         if(vcode == 4001)then
            nomvar_metric="GZ"
            fact=.1
            fact_8=.1d0
         else
            print*,'ERROR in test: in test reference field, Vcode not supported: ',trim(rfld_S)
            return
         endif
      else
         print*,'ERROR in test: in test reference field not supported: ',trim(rfld_S)
         return
      endif
      key=fstinf(lu,ni,nj,nk,datev,' ',ip1,ip2,-1,typvar,nomvar_metric)     
      ier = fstluk(px,key,ni,nj,nk)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstinf on',trim(nomvar_metric), ' , ip1=',ip1
         print*,'FILE: ',trim(F_fst)
         return
      endif
      ok=.true.
      if( abs(px(i0,j0)) < epsilon )then
         ! Value is probably zero
         if( abs(px(i0,j0)-pres(k)*fact) > epsilon)ok=.false.
      else
         if(abs((px(i0,j0)-pres(k)*fact)/abs(px(i0,j0)))>epsilon)ok=.false.
      endif
      if(.not.OK)then
         if(trim(rfld_S) == "P0")then
            print*,'ERROR in test: 32 bits: Difference in pressure is too large at'
            print*,'k,px(i0,j0),pres(k)*fact.',k,px(i0,j0),pres(k)*fact
         elseif(trim(rfld_S) == "ME")then
            print*,'ERROR in test: 32 bits: Difference in heights is too large at'
            print*,'k,gz(i0,j0),heights(k)*fact',k,px(i0,j0),pres(k)*fact
         else
            print*,'ERROR in test rfld_S ',trim(rfld_S),' no supported'
         endif         
         print*,'ERROR in test: on file: ',F_fst
         return
      endif
   enddo

   ! Test 64 bits interface
   if(two_refs_L)then
      ier = vgd_levels(vgd,ip1s,pres_8,p0_point_8,sfc_field_ls=p0ls_point_8)
   else
      ier = vgd_levels(vgd,ip1s,pres_8,p0_point_8)
   endif
   if(ier == VGD_ERROR )then
      print*,'ERROR in test: Problem with vgd_levels 64 bits'
      print*,'FILE: ',trim(F_fst)
      return
   endif 
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstprm on ',F_var,' ip1=',ip1
         print*,'FILE: ',trim(F_fst)
         return
      endif
      ! Surface pressure must be equal to P0 for surface levels
      call convip(ip1,pppp,kind,-1,dummy_S,.false.)
      if(abs(pppp-1.).lt.epsilon.or. ( abs(pppp).lt.epsilon .and. trim(rfld_S) == "ME") )then
         if(pres_8(k).ne.p0_point_8)then
            print*,'ERROR in test: 64 bits level at surface must be exacly equal to ',trim(rfld_S)
            print*,'File:',F_fst
            print*,'k,pres_8(k),p0_point_8',k,pres_8(k),p0_point_8
            return
         endif
      endif
      call incdatr(datev,dateo,deet*npas/3600.d0)
      key=fstinf(lu,ni,nj,nk,datev,' ',ip1,ip2,-1,typvar,nomvar_metric)     
      ier = fstluk(px,key,ni,nj,nk)
      if(ier.lt.0)then
         print*,'ERROR in test: Problem with fstinf on ',trim(nomvar_metric),', ip1=',ip1
         print*,'FILE: ',trim(F_fst)
         return
      endif
      ok=.true.
      if( abs(px(i0,j0)) < epsilon )then
         ! Value is probably zero
         if( abs(px(i0,j0)-pres_8(k)*fact) > epsilon)ok=.false.
      else
         if(abs((px(i0,j0)-pres_8(k)*fact_8)/abs(px(i0,j0)))>epsilon)ok=.false.
      endif
      if(.not.ok)then
         if(trim(rfld_S) == "P0")then
            print*,'ERROR in test: 46 bits: Difference in pressure is too large at'
            print*,'k,px(i0,j0),pres_8(k)*fact_8',k,px(i0,j0),pres_8(k)*fact_8
         endif
         if(trim(rfld_S) == "ME")then
            print*,'ERROR in test: 64 bits: Difference in heights is too large at'
            print*,'k,gz(i0,j0),heights_8(k)*fact_8',k,px(i0,j0),pres_8(k)*fact_8
         endif
         print*,'ERROR in test: on FILE: ',trim(F_fst)
         return
      endif
   enddo
   deallocate(ip1s,px,p0,pres,pres_8)

   print*,'   TEST OK!'

   ier=fstfrm(lu)

   status=VGD_OK   
   
end function check_levels_withref

