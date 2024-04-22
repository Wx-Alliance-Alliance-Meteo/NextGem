#include "vgrid_build_info.h"

program convert_toctoc_5002
   use app
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,vgd_print,VGD_OK
   use, intrinsic :: iso_fortran_env
   implicit none
   integer, parameter :: ncle=3,lui=10,nmax=1000
   integer, dimension(nmax) :: liste
   integer, dimension(2*nmax) :: key_to_erase
   integer :: npos,istat,fnom,fstouv,fstfrm,infon,fstprm,fstinl,fstluk,fstinf,fstecr,fsteff,ind_to_erase
   integer :: i,key,luo
   character(len=256), dimension(ncle) :: cle,val,def
   logical :: info_L
   real, allocatable, dimension(:,:) :: tocsf,p0
   real(kind=REAL64), allocatable, dimension(:,:) :: toctoc
   real :: dummy
   real(kind=REAL64) :: nhours
   logical :: samefile_L
   character(len=12), parameter :: version='1.1.0'

   ! For fstprm
   integer :: ig1,ig2,ig3,ig4,dateo,deet,npas,datyp,nbits
   integer :: ni,nj,nk,ni2,nj2,nk2,nkmod
   integer :: ip1,ip2,ip3,swa,lng,dltf,ubc,extra1,extra2,extra3,datev
   real(kind=REAL64), dimension(:), pointer :: a_m_8, b_m_8, a_t_8, b_t_8
   character(len=1) :: grtyp
   character(len=2) :: typvar
   character(len=4) :: nomvar
   character(len=12) :: etiket
   type(vgrid_descriptor) :: gdv
   
   cle=(/'s.      ','d.      ','samefile'/)
   val=(/'undef   ','undef   ','NO      '/)
   def=(/'undef   ','undef   ','YES     '/)

   nullify(a_m_8, b_m_8, a_t_8, b_t_8)

   ind_to_erase=0

   !==========================================================================
   app_ptr=app_init(0,'r.convert_toc_toc_5002',version,'Display verticalcoordinate information',BUILD_TIMESTAMP)
   call app_start()
   !==========================================================================

   ! Get keys
   npos=-111
   call ccard(cle,def,val,ncle,npos) 
   
   info_L=.false.
   if(val(1).eq.'undef')then
      write(app_msg,*) 'usage :'//EOL//&
         '      r.convert_toctoc_5002 -s source_file -d destination_file'//EOL//&
         '   or r.convert_toctoc_5002 -s source_file -samefile'
      call app_log(APP_VERBATIM,app_msg)
      istat=app_end(-1)
      stop
   endif
   
   if(val(2).eq.'undef'.and.val(3).eq.'NO')then
      write(app_msg,*) 'usage :'//EOL//&
         '      r.convert_toctoc_5002 -s source_file -d destination_file'
      call app_log(APP_VERBATIM,app_msg)
      istat=app_end(-1)
      stop
   endif
   
   samefile_L=val(3).eq.'YES'

   !==========================================================================
   ! Open files
   if(samefile_L)then
      istat=fnom(lui,val(1),'RND',0)
   else
      istat=fnom(lui,val(1),'RND+R/O',0)
   endif
   if(istat.lt.0)then
      istat=app_end(-1)
      error stop 1
   endif
   istat=fstouv(lui,'RND')
   if(istat.le.0)then
      istat=fstfrm(lui)
      istat=app_end(-1)
      error stop 1
   endif

   !==========================================================================
   ! Look for !!SF if not present then no convertion is needed
   istat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','!!SF',liste,infon,nmax)
   if(istat.lt.0)then
      istat=app_end(-1)
      error stop 1
   endif
   if(infon.eq.0)then
      call app_log(APP_WARNING,'Noting to convert')
      istat=fstfrm(lui)
      istat=app_end(-1)
      stop
   endif
   !
   ! Open output file
   if(samefile_L)then
      luo=lui
   else
      luo=lui+1
      istat=fnom(luo,val(2),'RND',0)
      if(istat.lt.0)then
         istat=app_end(-1)
         error stop 1
      endif
      istat=fstouv(luo,'RND')
      if(istat.lt.0)then
         istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
      endif 
   endif
   ! Loop on !!SF found
   do i=1,infon
      istat=fstprm(liste(i),dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,&
           ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
           ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)
      if(allocated(tocsf))deallocate(tocsf)
      allocate(tocsf(ni,nj))
      istat=fstluk(tocsf,liste(i),ni,nj,nk)
      if(istat.lt.0)then
         istat=fstfrm(lui)
         if(samefile_L)istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
      endif
      ! Look for matching P0, if not there write it
      nhours=deet*npas/3600.d0     
      call incdatr(datev,dateo,nhours)
      key=fstinf(lui,ni2,nj2,nk2,datev,' ',ip1,ip2,ip3,' ','P0')
      if(key.le.0)then
         ! Compute P0 form !!SF
         write(app_msg,*) 'Computing P0 from !!SF for hour ',nhours
         call app_log(APP_INFO,app_msg)
         if(allocated(p0))deallocate(p0)
         allocate(p0(ni,nj))
         p0=exp(tocsf)*1000.
         istat=fstecr(p0,dummy,-32,luo,dateo,deet,npas,ni,nj,nk,&
              ip1,ip2,ip3,typvar,'P0',etiket,grtyp,ig1,ig2,ig3,ig4,datyp,.true.)
      endif
      if(samefile_L)then
         ind_to_erase=ind_to_erase+1
         key_to_erase(ind_to_erase)=liste(i)
      endif
   enddo
   !
   !==========================================================================
   ! Find all !! and convert them
   istat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','!!',liste,infon,nmax)
   if(istat.lt.0)then
      istat=app_end(-1)
      error stop 1
   endif
   if(infon.lt.1)then
      write(app_msg,*) 'ERROR: cannot find !! in file ',val(1)
      call app_log(APP_ERROR,app_msg)
      istat=app_end(-1)
      error stop 1
   endif
   ! Loop on !! found
   do i=1,infon
      istat=fstprm(liste(i),dateo,deet,npas,ni,nj,nk,nbits,datyp,ip1,&
           ip2,ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3, &
           ig4,swa,lng,dltf,ubc,extra1,extra2,extra3)      
      ! Decode !! and reencode it to the format 5002      
      if(allocated(toctoc))deallocate(toctoc)
      allocate(toctoc(ni,nj))
      istat=fstluk(toctoc,liste(i),ni,nj,nk)
      if(istat.lt.0)then
         istat=fstfrm(lui)
         if(.not.samefile_L)istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
      endif
      nkmod=(nj-3)/2
      write(app_msg,*) 'njmod',nkmod
      call app_log(APP_DEBUG,app_msg)

      ! Looks like last momentum level is badly described in !! !!sf format
      toctoc(1,nkmod+1)=toctoc(1,nj)
      toctoc(2,nkmod+1)=toctoc(2,nj)
      toctoc(3,nkmod+1)=toctoc(3,nj)

      allocate(a_m_8(nkmod+1),b_m_8(nkmod+1),a_t_8(nkmod+2),b_t_8(nkmod+2))
      a_m_8 = toctoc(2,1:nkmod+1)
      b_m_8 = toctoc(3,1:nkmod+1)
      a_t_8 = toctoc(2,nkmod+2:nj)
      b_t_8 = toctoc(3,nkmod+2:nj)
      istat=vgd_new(&
           gdv,&            
           kind     = 5,&
           version  = 2,&
           nk       = nkmod,&
           ip1      = ip1,&
           ip2      = ip2,&
           ptop_8   = dble(ig2),&
           pref_8   = 100000.d0,&
           rcoef1   = real(ig3)/1000.,&    
           rcoef2   = real(ig4)/1000.,&
           a_m_8    = a_m_8,&
           b_m_8    = b_m_8,&   
           a_t_8    = a_t_8,&
           b_t_8    = b_t_8,&
           ip1_m    = nint(toctoc(1,1:nkmod+1)),&
           ip1_t    = nint(toctoc(1,nkmod+2:nj)))
      if(istat.ne.VGD_OK)then
         istat=fstfrm(lui)
         if(.not.samefile_L)istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
      endif
      istat = vgd_print(gdv)
      call app_log(APP_INFO,'Writing converted !!')
      istat = vgd_write(gdv,luo,'fst')
      if(istat.ne.VGD_OK)then
         istat=fstfrm(lui)
         if(.not.samefile_L)istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
         endif
      if(samefile_L)then
         ind_to_erase=ind_to_erase+1
         key_to_erase(ind_to_erase)=liste(i)
      endif
   enddo
   
   do i=1,ind_to_erase
      istat=fsteff(key_to_erase(i))
      if(istat.lt.0)then
         call app_log(APP_WARNING,'cannot erase record')
      endif
   enddo

   !==========================================================================
   ! close file(s)
   istat=fstfrm(lui)
   if(.not.samefile_L)istat=fstfrm(luo)
   istat=app_end(-1)
end program convert_toctoc_5002
