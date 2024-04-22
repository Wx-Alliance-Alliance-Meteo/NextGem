#include "vgrid_build_info.h"

program vcode
   !
   use app
   use vGrid_Descriptors, only: vgd_print,VGD_ERROR
   !
   implicit none
   !
   integer, parameter :: lu=10,ncle=1
   integer :: npos,stat,code
   character(len=12), parameter :: version='1.1.0'
   character(len=256), dimension(ncle) :: cle,val,def
  
   !==========================================================================
   app_ptr=app_init(0,'r.vcode',version,'Display vertical coordinate information',BUILD_TIMESTAMP)
   call app_start()
   !==========================================================================
   
   ! Get keys
   cle=(/'-'/)
   val=(/' '/)
   def=(/' '/)
   npos=1
   call ccard(cle,def,val,ncle,npos) 
   
   if(npos.eq.0)val='-1'
   read(val(1),*)code
   stat=vgd_print(code)
   
   stat=app_end(-1)

   
end program vcode
