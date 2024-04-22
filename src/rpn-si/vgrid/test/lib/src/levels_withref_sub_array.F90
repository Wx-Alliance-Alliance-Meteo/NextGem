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
module mod_levels_withref_sub_array
contains 
   integer function test_levels(levels,l_ni,l_nj,G_nk,liste,infon,lu) result(status)
      use vgrid_descriptors, only: VGD_ERROR, VGD_OK
      implicit none
      real, dimension(:,:,:), pointer :: levels
      integer :: infon
      integer :: l_ni,l_nj,G_nk,lu, liste(infon), datev
      ! Local variables
      integer :: i, j, k, fstprm, fstluk, fstinf, ier
      real, dimension(:,:), pointer :: wk
      real :: epsilon=5.0e-4
      ! Variable for fstprm, sorry...
      integer ::dateo, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
           ig2, ig3, ig4, ip1, ip2, ip3, key, lng, nbits,&
           nii,  njj, nkk, npas, swa, ubc
      character(len=12) :: etiket
      character(len=4)  :: nomvar
      character(len=2)  :: typvar
      character(len=1)  :: grtyp
      
      status = VGD_ERROR

      nullify(wk)

      if(infon > G_nk)then
         print*,'nk size error in test'
         return
      endif

      do k=1,infon
         ier = fstprm(liste(k),dateo,deet,npas,nii,njj,nkk,nbits,datyp, &
              ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
              dltf,ubc,extra1,extra2,extra3)
         call incdatr(datev,dateo,deet*npas/3600.d0)
         key = fstinf(lu,nii,njj,nkk,datev,' ',ip1,ip2,-1,typvar,'PX')
         if(key < 0)then
            print*,'Cannot find PX for datev,ip1,ip2,typvar',datev,ip1,ip2,typvar
            error stop 1
         endif
         if(k == 1) allocate(wk(nii,njj))
         ier = fstluk(wk,key,nii,njj,nkk)
         wk = wk*100.
         do j=1,l_nj
            do i=1,l_ni
               if(abs((wk(i,j)-levels(i,j,k))/wk(i,j))>epsilon)then
                  print*,'(Test) Difference in pressure is too large at'
                  print*,'i,j,k,px(i,j),levels(i,j,k)',i,j,k,wk(i,j),levels(i,j,k)
                  error stop 1
               endif
            enddo
         enddo
      end do      
      deallocate(wk)
      status = VGD_OK      
   end function test_levels
   
end module mod_levels_withref_sub_array
!========================================================
program levels_withref_sub_array
   use vgrid_descriptors, only: vgrid_descriptor, vgd_new, vgd_get, &
        vgd_levels, VGD_ERROR
   use mod_levels_withref_sub_array, only : test_levels
   use Unit_Testing, only: ut_report
   implicit none
   integer, parameter :: lu=10
   integer, parameter :: nmax=1000
   integer, dimension(nmax) :: liste
   integer :: fnom, fstouv, fstinf, fstluk, fstinl, fstprm
   integer :: ier, infon, k, l_ni, l_nj, G_nk
   integer, dimension(:), pointer :: ip1s
   real, dimension(:,:), pointer :: p0, p0ls, wk
   real, dimension(:,:), pointer :: p0_p, p0ls_p   
   real, dimension(:,:,:), pointer :: levels, levels_p
   character(len=128) :: fst_S
   character(len=4) :: rfld_S, rfld_ls_S
   type(vgrid_descriptor) :: vgd
   ! Variable for fstprm, sorry...
   integer ::dateo, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
        ig2, ig3, ig4, ip1, ip2, ip3, key, lng, nbits,&
        ni,  nj, nk, npas, swa, ubc
   character(len=12) :: etiket
   character(len=4)  :: nomvar
   character(len=2)  :: typvar
   character(len=1)  :: grtyp

   nullify(p0, p0_p, p0ls, p0ls_p, wk)

   fst_S="data/dm_5100_from_model_run"

   if(fnom(lu,fst_S,"RND",0) < 0)then
      print*,'(Test) ERROR with fnom on ',trim(fst_S)
      error stop 1
   endif
   if(fstouv(lu,'RND') < 0)then
      print*,'(Test) No record in RPN file ',trim(fst_S)
      error stop 1
   endif
   
   ier = fstinl(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ','TT',liste,infon,nmax)
   l_ni=ni; l_nj=nj; G_nk=infon
   if(infon == 0 )then
      print*,'(Test) pas de record de TT'      
      error stop 1
   endif
   allocate(ip1s(infon))
   do k=1,infon
      ier = fstprm(liste(k),dateo,deet,npas,ni,nj,nk,nbits,datyp, &
           ip1,ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng, &
           dltf,ubc,extra1,extra2,extra3)
      ip1s(k)=ip1
   end do
   allocate(wk(l_ni,l_nj),p0(0:l_ni+1,0:l_nj+1),p0ls(0:l_ni+1,0:l_nj+1), &
        levels(0:l_ni+1,0:l_nj+1,1:G_nk))
   !allocate(wk(l_ni,l_nj),p0(l_ni,l_nj),p0ls(l_ni,l_nj), &
   !     levels(l_ni,l_nj,G_nk))

   if(vgd_new(vgd,unit=lu,format="fst") == VGD_ERROR )then
      print*,'(Test) Problem getting vertical grid descriptor'
      error stop 1
   endif        
   ier = vgd_get(vgd,'RFLD',rfld_S)
   ier = vgd_get(vgd,'RFLS',rfld_ls_S)   
   print*,"rfld_S='",rfld_S,"'"
   print*,"rfld_ls_S='",rfld_ls_S,"'"
   
   key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfld_S)
   if(key < 0)then
      print*,'(Test) Problem getting ',rfld_S
      error stop 1
   endif
   if( fstluk(wk,key,ni,nj,nk) < 0)then
      print*,'(Test) Problem with fstluk on ',rfld_S
      error stop 1
   endif
   p0 = -99.
   p0(1:l_ni,1:l_nj) = wk(1:l_ni,1:l_nj)*100.
   key = fstinf(lu,ni,nj,nk,-1,' ',-1,-1,-1,' ',rfld_ls_S)
   if(key < 0)then
      print*,'(Test) Problem getting ',rfld_ls_S
      error stop 1
   endif
   if( fstluk(wk,key,ni,nj,nk) < 0)then
      print*,'(Test) Problem with fstluk on ',rfld_ls_S
      error stop 1
   endif
   p0ls = -99.
   p0ls(1:l_ni,1:l_nj) = wk(1:l_ni,1:l_nj)*100.
   
   p0_p => p0(1:l_ni,1:l_nj)
   p0ls_p => p0ls(1:l_ni,1:l_nj)
   levels_p => levels(1:l_ni,1:l_nj,1:G_nk)
   if( vgd_levels(vgd,ip1s,levels_p,p0_p,sfc_field_ls=p0ls_p) == VGD_ERROR)then
      error stop 1
   endif
   if( test_levels(levels_p,l_ni,l_nj,G_nk,liste,infon,lu) == VGD_ERROR)then
      print*,'ERROR in test with vgd_levels sub array on levels, p0 and p0ls'
      error stop 1
   endif
   !Test profile
   p0_p => p0(1:1,1:1)
   p0ls_p => p0ls(1:1,1:1)
   levels_p => levels(1:1,1:1,1:G_nk)
   if( vgd_levels(vgd,ip1s,levels_p,p0_p,sfc_field_ls=p0ls_p) == VGD_ERROR)then
      error stop 1
   endif
   if( test_levels(levels_p,1,1,G_nk,liste,infon,lu) == VGD_ERROR)then
      print*,'ERROR in test with vgd_levels profile sub array on levels, p0 and p0ls'
      error stop 1
   endif

   p0_p => p0(1:l_ni,1:l_nj)
   p0ls_p => p0ls(1:l_ni,1:l_nj)
   levels_p => levels(1:l_ni,1:l_nj,1:G_nk)
   ! Test mix sub array and array
   deallocate(levels)
   allocate(levels(l_ni,l_nj,G_nk))
   if( vgd_levels(vgd,ip1s,levels,p0_p,sfc_field_ls=p0ls_p) == VGD_ERROR)then
      error stop 1
   endif
   if( test_levels(levels,l_ni,l_nj,G_nk,liste,infon,lu) == VGD_ERROR)then
      print*,'ERROR in test with vgd_levels sub array on p0 and p0ls'
      error stop 1
   endif   
   wk(1:l_ni,1:l_nj)=p0(1:l_ni,1:l_nj)
   if( vgd_levels(vgd,ip1s,levels,wk,sfc_field_ls=p0ls_p) == VGD_ERROR)then
      error stop 1
   endif   
   if( test_levels(levels,l_ni,l_nj,G_nk,liste,infon,lu) == VGD_ERROR)then
      print*,'ERROR in test with vgd_levels sub array on p0ls'
      error stop 1
   endif
   wk(1:l_ni,1:l_nj)=p0ls(1:l_ni,1:l_nj)
   if( vgd_levels(vgd,ip1s,levels,p0_p,sfc_field_ls=wk) == VGD_ERROR)then
      error stop 1
   endif   
   if( test_levels(levels,l_ni,l_nj,G_nk,liste,infon,lu) == VGD_ERROR)then
      print*,'ERROR in test with vgd_levels sub array on p0'
      error stop 1
   endif

   call ut_report(.true.,'Grid_Descriptors, sub array')
end program levels_withref_sub_array
!===================================================================
