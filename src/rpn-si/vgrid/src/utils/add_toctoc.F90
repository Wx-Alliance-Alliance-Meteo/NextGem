#include "vgrid_build_info.h"

program add_toctoc
   use app
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_write,vgd_print,vgd_putopt,VGD_ERROR
   implicit none
   integer, parameter :: ncle=5,lui=10,nmax=1000
   integer, dimension(nmax) :: liste
   integer :: fnom,fstouv,infon,fstinl,fstfrm,ind_to_erase
   integer :: npos,istat,luo,ni,nj,nk,kind
   character(len=256), dimension(ncle) :: cle,val,def
   logical :: samefile_L, allow_sigma_L
   character(len=12), parameter :: version='2.2.0'
   type(vgrid_descriptor) :: vgd

   cle=(/'s.         ','d.         ','kind       ','samefile   ','allow_sigma'/)
   val=(/'undef      ','undef      ','undef      ','NO         ','NO         '/)
   def=(/'undef      ','undef      ','undef      ','YES        ','YES        '/)

   ind_to_erase=0

   !==========================================================================
   app_ptr=app_init(0,'r.add_toctoc',version,'Add vertical descriptor (!!)',BUILD_TIMESTAMP)
   call app_start()
   !==========================================================================
   
   ! Get keys
   npos=-111
   call ccard(cle,def,val,ncle,npos) 
   
   if(val(1).eq.'undef')then
      call help(version)
      istat=app_end(-1)
      stop
   endif
   
   if(val(2).eq.'undef'.and.val(3).eq.'NO')then
      call help(version)
      istat=app_end(-1)
      stop
   endif

   kind=-1
   if(val(3).ne.'undef')then
      read(val(3),*)kind
      write(app_msg,*) 'Only adding level descriptor for kind',kind
      call app_log(APP_INFO,app_msg)
   endif  

   samefile_L=val(4).eq.'YES'

   allow_sigma_L = trim(val(5)) == "YES"
   call app_log(APP_INFO,'Construction of sigma coordinate is allowed')

   !==========================================================================
   ! Open files
   if(samefile_L)then
      istat=fnom(lui,val(1),'RND',0)
   else
      istat=fnom(lui,val(1),'RND+R/O',0)
   endif
   if(istat.lt.0)then
      write(app_msg,*) 'Problem with fnom on file ',trim(val(1))
      call app_log(APP_ERROR,app_msg)
      istat=app_end(-1)
      error stop 1
   endif
   istat=fstouv(lui,'RND')
   if(istat.le.0)then
      write(app_msg,*) 'No record in RPN file ',trim(val(1))
      call app_log(APP_ERROR,app_msg)
      istat=fstfrm(lui)
      istat=app_end(-1)
      error stop 1
   endif

   !==========================================================================
   ! Look for !!SF if present then tell user to use convert_toctoc_5002
   istat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','!!SF',liste,infon,nmax)
   if(infon.gt.0)then
      call app_log(APP_ERROR,'record !!SF is present, use convert_toctoc_5002 instead')
      istat=app_end(-1)
      error stop 1
   endif
   !
   istat=fstinl(lui,ni,nj,nk,-1,' ',-1,-1,-1,' ','!!',liste,infon,nmax)
   if(infon.gt.0)then
       call app_log(APP_INFO,'Record !! already there, noting to do')
       goto 999
   endif
   !
   ! Open output file
   if(samefile_L)then
      luo=lui
   else
      luo=lui+1
      istat=fnom(luo,val(2),'RND',0)
      if(istat.lt.0)then
         write(app_msg,*) 'Problem with fnom on file ',trim(val(2))
         call app_log(APP_ERROR,app_msg)
         istat=app_end(-1)
         error stop 1
      endif
      istat=fstouv(luo,'RND')
      if(istat.lt.0)then
         write(app_msg,*) 'Problem with fstouv on ',trim(val(2))
         call app_log(APP_ERROR,app_msg)
         istat=fstfrm(luo)
         istat=app_end(-1)
         error stop 1
      endif 
   endif

   istat = vgd_putopt("ALLOW_SIGMA",allow_sigma_L)
   if(istat.eq.VGD_ERROR)then
      istat=app_end(-1)
      error stop 1
   endif

   istat=vgd_new(vgd,lui,kind=kind)
   if(istat.eq.VGD_ERROR)then
      istat=app_end(-1)
      error stop 1
   endif

   istat=vgd_print(vgd)
   if(istat.eq.VGD_ERROR)then
      istat=app_end(-1)
      error stop 1
   endif

   istat=vgd_write(vgd,luo)
   if(istat.eq.VGD_ERROR)then
      istat=app_end(-1)
      error stop 1
   endif

   !==========================================================================
   ! close file(s)
   if(.not.samefile_L)istat=fstfrm(luo)
999 continue
   istat=fstfrm(lui)
   istat=app_end(-1)
end program add_toctoc

subroutine help (version)
   use app
   implicit none
   character(len=*) :: version
   integer :: istat
   write(app_msg,*) 'Usage  :'//EOL//&
      '   r.add_toctoc -s source_file -d destination_file'//EOL//&
      ' Options:'//EOL//&
      '   -kind         narrow selection of ip1 of kind kind in reconstruction'//EOL//&
      '                 This is handy, and some times necessairy when there'//EOL//&
      '                 are mixed kind in the input file'//EOL//&
      '   -samefile     put constructed !! in input file (-s input_file),'//EOL//&
      '                 do not supply argument -d with this option'//EOL//&
      '   -allow_sigma  allow construction of sigma levels, default to false'
   call app_log(APP_VERBATIM,app_msg)

   end subroutine help
