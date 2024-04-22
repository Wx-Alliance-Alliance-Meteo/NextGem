program gem
   use app
   use iso_fortran_env
   implicit none

#include <gem_build_info.h>
#include <gemdyn_version.inc>
#include <rpnphy_version.inc>
#include <modelutils_version.inc>

   integer(kind=int32) ierror
   integer :: colors(3), COMMs(3), wnum, wme, cnum, cme, un_out
   character(len=256) :: component_S

   app_ptr=app_init(0,PROJECT_NAME_STRING,VERSION,PROJECT_DESCRIPTION_STRING,BUILD_TIMESTAMP)
   call app_libregister(APP_LIBVGRID,HAVE_VGRID)
   call app_libregister(APP_LIBTDPACK,HAVE_TDPACK)
   call app_libregister(APP_LIBGEMDYN,GEMDYN_VERSION_S)
   call app_libregister(APP_LIBRPNPHY,RPNPHY_VERSION_S)
   call app_libregister(APP_LIBMDLUTIL,MODELUTILS_VERSION_S)

!MPI_THREAD_SINGLE: Only one thread will execute. 
!MPI_THREAD_FUNNELED: The process may be multi-threaded, but only the main thread will make MPI calls (all MPI calls are funneled to the main thread). 
!MPI_THREAD_SERIALIZED: The process may be multi-threaded, and multiple threads may make MPI calls, but only one at a time: MPI calls are not made concurrently from two distinct threads (all MPI calls are serialized). 
!MPI_THREAD_MULTIPLE: Multiple threads may call MPI, with no restrictions.
   component_S= 'GEMDM' ; colors= (/1,3,4/)
   call MiMd_init (component_S,colors,3,COMMs,wnum, wme, cnum, cme)
   un_out= -1
   if ( cme == 0 ) un_out= 6
   call gemtime ( 6, ' ', .false. )
   if (un_out>0) call gemtime ( un_out, 'STARTING GEMDM DOMAINS', .false. )

   call app_start()
 
   ! Initialize: Domain, MPI, processor topology and ptopo.cdk
   call init_component(COMMs,3)

   ! Establish: model configuration, domain decomposition and model geometry
   call set_world_view()
  
   ! Initialize the ensemble prevision system
   call itf_ens_init()
  
   ! Initialize the physics parameterization package
   call itf_phy_init()
  
   ! Initialize tracers
   call tracers()
  
   ! Setup main memory
   call main_gmm_storage()
   call set_dyn_opr()
  
   ! Run GEM
   call gem_ctrl()
  
   ! Terminate
   call stop_world_view()
   if (un_out>0) call gemtime ( un_out, 'ENDING GEMDM DOMAINS', .true. )
   call MPI_FINALIZE(ierror)  

   app_status=app_end(-1)
end program gem

