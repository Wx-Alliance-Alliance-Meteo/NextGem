#include "vgrid_build_info.h"

module mod_vgrid_sample
  use app
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: gen_it, sample_vgrid
  type vcode_eta
     real(kind=REAL64) :: ptop_8 = -1.d0
     character(len=12) :: etiket = 'ETA'
  end type vcode_eta
  type vcode_5001
     real :: rcoef = 1.6
     character(len=12) :: etiket = 'HYB_REGULAR'
  end type vcode_5001
  type vcode_5002
     real :: rcoef1 = 0., rcoef2 = 5.
     character(len=12) :: etiket = 'HYB_STG_CP'
  end type vcode_5002
  type vcode_5005
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., dhm=10., dht=1.5
     character(len=12) :: etiket = 'HYB_CP'
  end type vcode_5005
  type vcode_5100
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., rcoef3 = 0., rcoef4 = 100., &
          dhm=10., dht=1.5
     character(len=12) :: etiket = 'HYB_CP_SLEVE'
  end type vcode_5100
  type vcode_21001_NON_SLEVE
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., &
          dhm=10., dht=1.5
     character(len=12) :: etiket = 'HYB_GC'
  end type vcode_21001_NON_SLEVE
  type vcode_21001_SLEVE
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., rcoef3 = 0., rcoef4 = 100., &
          dhm=10., dht=1.5
     character(len=12) :: etiket = 'HYB_GC_SLEVE'
  end type vcode_21001_SLEVE
  type vcode_21002_NON_SLEVE
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., &
          dhm=10., dht=1.5, dhw=1.
     character(len=12) :: etiket = 'HYB_LZ'
  end type vcode_21002_NON_SLEVE
  type vcode_21002_SLEVE
     real :: hyb_flat = -1.
     real :: rcoef1 = 0., rcoef2 = 5., rcoef3 = 0., rcoef4 = 100., &
          dhm=10., dht=1.5, dhw=1.
     character(len=12) :: etiket = 'HYB_LZ_SLEVE'
  end type vcode_21002_SLEVE
  
  type(vcode_eta),  public :: vc_eta
  type(vcode_5001), public :: vc_5001
  type(vcode_5002), public :: vc_5002
  type(vcode_5005), public :: vc_5005
  type(vcode_5100), public :: vc_5100
  type(vcode_21001_NON_SLEVE), public :: vc_21001_NON_SLEVE
  type(vcode_21001_SLEVE), public :: vc_21001_SLEVE
  type(vcode_21002_NON_SLEVE), public :: vc_21002_NON_SLEVE
  type(vcode_21002_SLEVE), public :: vc_21002_SLEVE
  integer :: nk
  real, dimension(1000) :: levs, levs_eta
  
contains

  integer function gen_it(F_dir, F_name) result(status)

    ! Andre Plante Sept 2018

    use vgrid_descriptors, only: vgrid_descriptor, vgd_new, vgd_putopt, &
         vgd_print, vgd_write, vgd_free, VGD_ERROR, VGD_OK, &
         vgd_stda76_hgts_from_pres_list, VGD_STDA76_SFC_P, vgd_put, &
         vgd_stda76, vgd_get

    use, intrinsic :: iso_fortran_env
    implicit none
    character(len=*) :: F_dir, F_name

    ! Local variables
    integer, save :: lu = 10
    real, dimension(nk) :: eta
    real, dimension(nk-1) :: hyb_m1, zeta, hgts, hgts_rev, pres
    real(kind=REAL64) :: ptop_8, ptop_out_8, pref_8
    real, dimension(:), pointer :: temp
    integer :: stat, fnom, fstouv, fstfrm, k
    integer, dimension(:), pointer :: ip1s
    type(vgrid_descriptor) :: vgd
    status = VGD_ERROR

    nullify(temp,ip1s)

    lu = lu + 1
    write(app_msg,*) 'lu =',lu
    call app_log(APP_DEBUG,app_msg)

    stat = fnom(lu, trim(F_dir)//'/'//trim(F_name), "RND", 0)
    if(stat < 0)then
      stat=app_end(-1)
       error stop 1
    endif
    stat = fstouv(lu, 'RND')
    if(stat < 0)then
      stat=app_end(-1)
       error stop 1
    endif

    pres = levs(1:nk-1) * VGD_STDA76_SFC_P
    ! Compute stda 1976 heights from pres pressure profile
    if(vgd_stda76_hgts_from_pres_list(hgts, pres, nk-1) == &
         VGD_ERROR)then
            stat=app_end(-1)
            error stop 1
    endif
    do k=1,nk-1
       hgts_rev(k) = hgts(nk-k)
    end do

    select case(trim(F_name))
    case('pressure')
       if( vgd_new(vgd, kind=2, version=1, hyb=levs(1:nk) * VGD_STDA76_SFC_P/100.) == VGD_ERROR)then
         write(app_msg,*) 'building ',trim(F_name)
         call app_log(APP_ERROR,app_msg)
         stat=app_end(-1)
         error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", "PRESSURE") == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif
    case('sigma')
       if( vgd_new(vgd, kind=1 , version=1, hyb=levs(1:nk)) == VGD_ERROR)then
         write(app_msg,*) 'building ',trim(F_name)
         call app_log(APP_ERROR,app_msg)
         stat=app_end(-1)
         error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", "SIGMA") == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif
    case('eta')
       eta = (levs(1:nk) - levs(1)) / (levs(nk) - levs(1))
       ptop_8 = levs(1) * VGD_STDA76_SFC_P * 1.d0
       if(vc_eta%ptop_8 < 0 )then
          vc_eta%ptop_8 = levs(1) * VGD_STDA76_SFC_P * 1.d0
       else
          if( abs(vc_eta%ptop_8 - ptop_8) > 1.e-5)then
            write(app_msg,*) 'value of eta ptop inconsistant ',trim(F_name),', expected ',ptop_8,' got ', vc_eta%ptop_8,&
             '. If you give levs in namelist do not set vc_eta%ptop_8 value'
            call app_log(APP_ERROR,app_msg)
            stat=app_end(-1)
            error stop 1
          endif
       endif
       if( vgd_new(vgd, kind=1 , version=2, ptop_8=vc_eta%ptop_8, hyb=eta) == &
            VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_eta%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif
    case('4001')
       if( vgd_new(vgd, kind=4 , version=1, hyb=hgts_rev) == &
            VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", "M ABV SFC") == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif
    case('5001')
       ! Note : here we want hyb level to be equal to sigam level at
       !        p0 = VGD_STDA76_SFC_P. If we set pref = VGD_STDA76_SFC_P, then hyb = levs.
       ptop_8=levs(1)*100000.
       pref_8 = VGD_STDA76_SFC_P * 1.d0
       if( vgd_new(vgd, kind=5 , version=1, ptop_8=ptop_8, pref_8=pref_8, &
            rcoef1=vc_5001%rcoef, hyb=levs(1:nk)) == &
            VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_5001%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

      case('5002')
       ! Note : here we want hyb level to be equal to sigma level at
       !        p0 = VGD_STDA76_SFC_P. If we set pref = VGD_STDA76_SFC_P,
       !        then zeta = ln(pres)
       !        zeta = zetas + ln(eta) => eta = levs
       hyb_m1 = levs(1:nk-1)
       zeta = log(hyb_m1*VGD_STDA76_SFC_P)
       ptop_8 = exp(zeta(1) - .5*(zeta(2) - zeta(1)))    
       pref_8 = VGD_STDA76_SFC_P * 1.d0
       if( vgd_new(vgd, kind=5 , version=2, ptop_8=ptop_8, pref_8=pref_8, &
            rcoef1=vc_5002%rcoef1, rcoef2=vc_5002%rcoef2, hyb=hyb_m1) == &
            VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_5002%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('5005')
       ! Note : here we want hyb level to be equal to sigam level at
       !        p0 = VGD_STDA76_SFC_P. If we set pref = VGD_STDA76_SFC_P,
       !        then zeta = ln(pres)
       !        zeta = zetas + ln(eta) => eta = levs
       hyb_m1 = levs(1:nk-1)
       zeta = log(hyb_m1*VGD_STDA76_SFC_P)
       pref_8 = VGD_STDA76_SFC_P
       if( vgd_new(vgd, kind=5 , version=5, ptop_out_8=ptop_out_8, &
            pref_8=pref_8, rcoef1=vc_5005%rcoef1, rcoef2=vc_5005%rcoef2, &
            hyb=hyb_m1, dhm=vc_5005%dhm, dht=vc_5005%dht, hyb_flat=vc_5005%hyb_flat) &
            == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_5005%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('5100')
       ! Note : here we want hyb level to be equal to sigam level at
       !        p0 = VGD_STDA76_SFC_P. If we set pref = VGD_STDA76_SFC_P,
       !        then zeta = ln(pres)
       !        zeta = zetas + ln(eta) => eta = levs
       hyb_m1 = levs(1:nk-1)
       zeta = log(hyb_m1*VGD_STDA76_SFC_P)
       pref_8 = VGD_STDA76_SFC_P
       if( vgd_new(vgd, kind=5 , version=100, ptop_out_8=ptop_out_8, &
            pref_8=pref_8, rcoef1=vc_5100%rcoef1, rcoef2=vc_5100%rcoef2, &
            rcoef3=vc_5100%rcoef3, rcoef4=vc_5100%rcoef4, &
            hyb=hyb_m1, dhm=vc_5100%dhm, dht=vc_5100%dht,hyb_flat=vc_5100%hyb_flat)  == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_5100%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('21001_NON_SLEVE')
       if( vgd_new(vgd, kind=21 , version=1, rcoef1=vc_21001_NON_SLEVE%rcoef1,&
            rcoef2=vc_21001_NON_SLEVE%rcoef2, rcoef3=-1., &
            rcoef4=-1., hyb=hgts, dhm=vc_21001_NON_SLEVE%dhm,&
            dht=vc_21001_NON_SLEVE%dht,hyb_flat=vc_21001_NON_SLEVE%hyb_flat) == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
                error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_21001_NON_SLEVE%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('21001_SLEVE')
       if( vgd_new(vgd, kind=21 , version=1, rcoef1=vc_21001_SLEVE%rcoef1,&
            rcoef2=vc_21001_SLEVE%rcoef2, rcoef3=vc_21001_SLEVE%rcoef3, &
            rcoef4=vc_21001_SLEVE%rcoef4, hyb=hgts, dhm=vc_21001_SLEVE%dhm,&
            dht=vc_21001_SLEVE%dht,hyb_flat=vc_21001_SLEVE%hyb_flat) == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
                stat=app_end(-1)
                error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_21001_SLEVE%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('21002_NON_SLEVE')
       if( vgd_new(vgd, kind=21 , version=2, rcoef1=vc_21002_NON_SLEVE%rcoef1,&
            rcoef2=vc_21002_NON_SLEVE%rcoef2, rcoef3=-1., &
            rcoef4=-1., hyb=hgts, dhm=vc_21002_NON_SLEVE%dhm,&
            dht=vc_21002_NON_SLEVE%dht, dhw=vc_21002_NON_SLEVE%dhw,hyb_flat=vc_21002_NON_SLEVE%hyb_flat)&
            == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
       endif
       if( vgd_put(vgd, "ETIKET", vc_21002_NON_SLEVE%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case('21002_SLEVE')
       if( vgd_new(vgd, kind=21 , version=2, rcoef1=vc_21002_SLEVE%rcoef1,&
            rcoef2=vc_21002_SLEVE%rcoef2, rcoef3=vc_21002_SLEVE%rcoef3, &
            rcoef4=vc_21002_SLEVE%rcoef4, hyb=hgts, dhm=vc_21002_SLEVE%dhm,&
            dht=vc_21002_SLEVE%dht, dhw=vc_21002_SLEVE%dhw,hyb_flat=vc_21002_SLEVE%hyb_flat) &
            == VGD_ERROR)then
               write(app_msg,*) 'building ',trim(F_name)
               call app_log(APP_ERROR,app_msg)
       stat=app_end(-1)
               error stop 1
       endif
       !if( vgd_get(vgd,'VIPM',ip1s) == VGD_ERROR)then
       !   print*,'Erreur with vgd_get on VIPM'
       !   error stop 1
       !endif
       !if( vgd_stda76(vgd,ip1s,temp,'PRESSURE',sfc_pres=100000.) == VGD_ERROR)then
       !   print*,'Erreur with vgd_stda76'
       !   error stop 1
       !endif
       !print*,'ip1s',ip1s
       !print*,'temp=',temp
       if( vgd_put(vgd, "ETIKET", vc_21002_SLEVE%etiket) == VGD_ERROR) then
         stat=app_end(-1)
         error stop 1
       endif

    case DEFAULT
      write(app_msg,*) 'Please add ', trim(F_name)
      call app_log(APP_INFO,app_msg)
       return
    end select

    if( vgd_write(vgd, lu, 'fst') == VGD_ERROR)then
       write(app_msg,*) 'writing vgrid for ',trim(F_name)
       call app_log(APP_ERROR,app_msg)
       stat=app_end(-1)
        error stop 1
    endif
    stat = vgd_print(vgd)
    stat = vgd_free(vgd)

    stat = fstfrm(lu)

    status = VGD_OK

  end function gen_it

  subroutine sample_vgrid()
    
    ! Andre Plante Sept 2018
    
   ! use mod_vgrid_sample, only: gen_it
    use vgrid_descriptors, only: vgd_putopt, VGD_STDA76_SFC_P, VGD_ERROR
    
    implicit none
    
    integer, parameter :: ncle=1, lunml = 9
    integer :: stat, npos
    integer :: exdb, exfin, system, k
    integer, dimension(1) :: dummy
    character(len=12), parameter :: version='1.1.0'
    character(len=32), parameter :: nml_file = 'vgrid_sample.nml'
    character(len=612), dimension(ncle) :: cle, val, def 
    namelist /cfg/ levs, levs_eta, vc_eta, vc_5001, vc_5002, vc_5005, vc_5100, &
       vc_21001_NON_SLEVE, vc_21001_SLEVE, vc_21002_NON_SLEVE, vc_21002_SLEVE
    
    cle = (/'out_dir.'/)
    val = (/'undef   '/)
    def = (/'undef   '/)
    
    !==========================================================================
    app_ptr=app_init(0,'r.vgrid_sample',version,'',BUILD_TIMESTAMP)
    call app_start()
    !==========================================================================
    
    npos = -111
    call ccard(cle, def, val, ncle, npos)
    
    if(trim(val(1)) == 'undef')then
      call app_log(APP_VERBATIM,'Usage: r.vgrid_sample -out_dir output_directory')
      stat=app_end(-1)
      stop
    endif
    
    if( system('mkdir -p '//trim(val(1))) /= 0 )then
      write(app_msg,*) 'Problem with directory ',trim(val(1)),' Is it a file, do you have write permission?'
      call app_log(APP_ERROR,app_msg)
      stat=app_end(-1)
      error stop 1
    endif
    
    
    levs = -1
    levs_eta = -1
    open (lunml, file=nml_file, delim='APOSTROPHE', STATUS='OLD', &
         iostat=stat)
    nk = 0
    if(stat /= 0)then
      write(app_msg,*) 'cannot open namlist file ',trim(nml_file),'. Using default vertical grid parameters'
    else
       read(lunml, nml=cfg)
       dummy=minloc(levs)
       nk = dummy(1) - 1
       if( nk /= 0)then
         call app_log(APP_INFO,'Using user defined sigma levels to set all vcode levels')
       else
          dummy=minloc(levs_eta)
          nk = dummy(1) - 1
          if( nk /= 0 )then
             call app_log(APP_INFO,'Using user defined eta levels to set all vcode levels')
             if(vc_eta%ptop_8 < 0.d0)then
               write(app_msg,*) 'you must set vc_eta%ptop_8 [Pa] in namelist when setting levs_eta'
               call app_log(APP_ERROR,app_msg)
               stat=app_end(-1)
               error stop 1
             endif
             do k = 1, nk
                levs(k) = levs_eta(k) *  (1.d0 - vc_eta%ptop_8/VGD_STDA76_SFC_P) + vc_eta%ptop_8/VGD_STDA76_SFC_P
             end do
          else
            call app_log(APP_INFO,'Using default sigma levels to set all vcode levels')
          endif
       endif
    end if
    if( nk == 0)then
       nk = 12
       levs(1:nk) = (/.01, .025, .05, .09, 0.15, 0.23, 0.35, 0.5, 0.65, 0.8, 0.9, 1.0/)
    endif
    
    stat = vgd_putopt("ALLOW_SIGMA",.true.)
    
    stat = gen_it(val(1), 'pressure')
    stat = gen_it(val(1), 'sigma')
    stat = gen_it(val(1), 'eta') 
    stat = gen_it(val(1), '4001') 
    stat = gen_it(val(1), '5001')
    stat = gen_it(val(1), '5002') 
    stat = gen_it(val(1), '5005') 
    stat = gen_it(val(1), '5100')
    stat = gen_it(val(1), '21001_NON_SLEVE') 
    stat = gen_it(val(1), '21001_SLEVE') 
    stat = gen_it(val(1), '21002_NON_SLEVE')
    stat = gen_it(val(1), '21002_SLEVE') 
    
    stat=app_end(-1)
    
  end subroutine sample_vgrid
end module mod_vgrid_sample
  
program call_sample_vgrid
    
  use mod_vgrid_sample, only : sample_vgrid
  call sample_vgrid()
  
end program call_sample_vgrid


!=====================================
