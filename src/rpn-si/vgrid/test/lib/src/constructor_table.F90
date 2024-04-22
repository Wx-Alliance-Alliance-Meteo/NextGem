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
module mod_constructor_table
  
  type FSTD_ext
     integer :: ig1,ig2,ig3,ig4,dateo,deet,npas,datyp,nbits,ni,nj,nk
     integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3,datev
     character(len=1) :: grtyp
     character(len=2) :: typvar
     character(len=4) :: nomvar
     character(len=12) :: etiket
  end type FSTD_ext
contains
   integer function my_fstprm(fstkey,record) result(status)
      ! Use fstprm function to get information about the record
      integer, intent(in) :: fstkey               !Key from FST file record
      type(FSTD_ext) :: record                    !Record information
      integer :: error
      integer, external :: fstprm,fstinf
      status = -1
      error=fstprm(fstkey,record%dateo,record%deet,record%npas, &
           record%ni,record%nj,record%nk,record%nbits,record%datyp,record%ip1,record%ip2, &
           record%ip3,record%typvar,record%nomvar,record%etiket,record%grtyp, &
           record%ig1,record%ig2,record%ig3,record%ig4,record%swa, &
           record%lng,record%dltf,record%ubc,record%extra1,record%extra2, &
           record%extra3)
      if (error < 0) then
         write(6,*) 'cannot fstprm for fstkey ',fstkey
         return
      end if
      !nhours=record%deet*record%npas/3600.d0
      !call incdatr(record%datev,record%dateo,nhours)
      status = 0
   end function my_fstprm
end module mod_constructor_table
!================================================================================================
program constructor
   use vGrid_Descriptors, only:VGD_OK
   use Unit_Testing, only: ut_report
   implicit none
   integer, parameter :: nversion=11
   integer stat,i,test_it
   logical :: OK=.true.
   character(len=10), dimension(nversion) :: vcode_S=(/"4001 ","5002 ","5003 ","5004 ","5005 ","5100 ","5999 ","21001", &
        "21001","21002","21002"/)
   character(len=60), dimension(nversion) :: suffix_S=(/"          ","          ","          " ,"          ","          ", &
        "          ","          " ,"_SLEVE    ","_NON_SLEVE","_SLEVE    ","_NON_SLEVE"/)

   do i=1,nversion
      
      print*,"============="//trim(vcode_S(i))//"=============="
      stat=test_it("data/dm_"//trim(vcode_S(i))//"_from_model_run"//trim(suffix_S(i)))
      if(stat.ne.VGD_OK)then
         OK=.false.
         print*,'ERROR with ',trim(vcode_S(i))
      endif
   enddo

   call ut_report(OK,'Grid_Descriptors, rebuild from table')    

end program constructor
!================================================================================================
integer function test_it(F_file_from_model_run) result(status)
 
   use mod_constructor_table, only:FSTD_ext,my_fstprm
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_get,vgd_free,operator(==),vgd_write,VGD_OK,VGD_ERROR
   use Unit_Testing, only: ut_report
   
   implicit none
   
   character(len=*) :: F_file_from_model_run

   ! Local variables

   type(vgrid_descriptor) :: vgrid,vgrid_rebuilt
   integer, parameter :: nmax = 1000
   integer, dimension(nmax) :: liste
   integer :: lu=0,lu2=0
   integer :: ier,stat,infon,i,ni,nj,nk,key
   integer :: fnom,fstouv,fstfrm,fstinl,fsteff,fstinf
   real(kind=8), dimension(:,:,:), pointer :: table
   type(FSTD_ext) :: prm_orig,prm_rebuilt
   
   status=VGD_ERROR
   
   nullify(table)

   ier=fnom(lu,trim(F_file_from_model_run),"RND",0)
   if(ier.lt.0)then
      print*,'ERROR with fnom'
      return
   endif
   ier=fstouv(lu,'RND')
   if(ier.le.0)then
      print*,'No record in RPN file ',trim(F_file_from_model_run)
      return
   endif   
   
   ! Construct a new set of 3D coordinate descriptors
   stat=0  
   if(VGD_OK.ne.vgd_new(vgrid,unit=lu,format="fst"))return
   if(VGD_OK.ne.vgd_get(vgrid,'VTBL',table))return
   if(VGD_OK.ne.vgd_new(vgrid_rebuilt,table))return
   if (vgrid_rebuilt == vgrid) then
      ! Ok do noting
   else
      print*,'ERROR: rebuilding table'
      return
   endif

   ! Write table in fst and check its fst parameters matches
   call system('rm -f trash')
   ier=fnom(lu2,"trash","RND",0)
   if(ier.lt.0)then
      print*,'ERROR with fnom 2'
      return
   endif
   ier=fstouv(lu2,'RND')
   ! Remove all possible records
   ier = fstinl(lu2,ni,nj,nk,-1,' ',-1,-1,-1,' ',' ',liste,infon,nmax)
   do i=1,infon
      ier = fsteff(liste(i))
   enddo

   if (vgd_write(vgrid_rebuilt,lu2,'fst') /= 0) then
      print*, 'ERROR with vgd_write'
      return
   endif

   ! Get parameters of original !!
   key = fstinf(lu ,ni,nj,nk,-1,' ',-1,-1,-1,' ','!!')
   ier = my_fstprm(key,prm_orig)
   if(ier.lt.0)then
      print*,'ERROR with my_fstprm 1'
      return
   endif
   
   ! !! is alone in its file, read it
   key = fstinf(lu2,ni,nj,nk,-1,' ',-1,-1,-1,' ',' ')
   ier = my_fstprm(key,prm_rebuilt)
   if(ier.lt.0)then
      print*,'ERROR with my_fstprm 2'
      return
   endif

   if(trim(prm_rebuilt%nomvar).ne.trim(prm_orig%nomvar))then
      print*,'nomvar do not match, expected ',trim(prm_orig%nomvar),' got ',trim(prm_rebuilt%nomvar)
      return
   endif
   !if(trim(prm_rebuilt%etiket).ne.trim(prm_orig%etiket))then
   !   print*,'Etiket do not match, expected ',trim(prm_orig%etiket),' got ',trim(prm_rebuilt%etiket)
   !   return
   !endif
   if(prm_rebuilt%datyp.ne.prm_orig%datyp)then
      print*,'datyp do not match, expected ',prm_orig%datyp,' got ',prm_rebuilt%datyp
      return
   endif
   
   deallocate(table)
   ier = vgd_free(vgrid)
   ier = vgd_free(vgrid_rebuilt)

   ier=fstfrm(lu)
   call fclos(lu)
   ier=fstfrm(lu2)
   call fclos(lu2)
   
   status=VGD_OK
   
end function test_it
