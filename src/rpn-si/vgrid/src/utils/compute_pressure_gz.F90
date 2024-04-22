#include "vgrid_build_info.h"

module mod_comp_pres_gz
  use app
  use vGrid_Descriptors, only: vgrid_descriptor, VGD_ERROR, VGD_OK
  implicit none
  private
  public :: cpg_print_doc,cpg_process_arguments,cpg_get_vgd_levels,cpg_rpn,cpg_get_rec,&
       cpg_get_px_gz
  integer, public, parameter :: cpg_lui=10
  integer, public :: cpg_luo=11,cpg_kind,cpg_version
  character(len=12), public :: cpg_etiket_S
  character(len=20), public :: cpg_levels_S
  logical, public :: cpg_samefile_L,cpg_allow_sigma_L,cpg_zap_etiket_L,&
       cpg_is_pressure_L,cpg_compute_gz_L
  type(vgrid_descriptor), public :: cpg_vgd
  type cpg_rpn
     real, dimension(:,:), pointer :: data=>null()
     integer ::dateo, datev, datyp, deet, dltf, extra1, extra2, extra3, ig1,&
          ig2, ig3, ig4, ip1, ip2, ip3, iun, key, lng, nbits,&
          ni,  nj, nk, npak, npas, swa, ubc
     character(len=12) :: etiket
     character(len=4)  :: nomvar
     character(len=2)  :: typvar
     character(len=1)  :: grtyp, ctype
     logical :: rewrit
  end type cpg_rpn
contains
  subroutine cpg_print_doc()
    write(6,'("   ")')
    write(6,'("   Usage:")')
    write(6,'("      Mandatory arguments:")')
    write(6,'("         -compute (pressure,gz)")')
    write(6,'("         -s source_file")')
    write(6,'("         -d destination_file (not mandatory if option -samefile is")')
    write(6,'("            used see below)")')
    write(6,'("   Options:")')
    write(6,'("         -samefile")')
    write(6,'("            Write results in input file, with this option argument -d must not be used.")')
    write(6,'("         -levels (NOMVAR,MOMENTUM,THERMO,VERTICAL_VELOCITY,ALL_LEVELS)")')
    write(6,'("            Nomvar or keyword to define which levels to get the results on.")')
    write(6,'("            If nomvar is given, e.g. -levels TT, then the computation will be done for all")')
    write(6,'("            TT levels in the input file.")')
    write(6,'("            One of the following keywords may be used:")')
    write(6,'("               MOMENTUM          to get the results on the momentum levels")')
    write(6,'("               THERMO            to get the results on the thermodynamic levels")')
    write(6,'("               VERTICAL_VELOCITY to get the results on the vertical-velocity levels")')
    write(6,'("               ALL_LEVELS        to get the results on all levels (default)")')
    write(6,'("         -kind (kind)")')
    write(6,'("            Usefull if there are more then one ip1 kind or !! in the input file. The kind value")')
    write(6,'("            passed to the program via this option is used to determine the taret levels ip1 kind")')
    write(6,'("            in the input file.")')
    write(6,'("         -version (vertical grid descriptor version number)")')
    write(6,'("            This is similar to the -kind option. The value passed to this option will be used")')
    write(6,'("            to target a specific vertical grid descriptor in the input file. option -kind must")')
    write(6,'("            also be used along with this option")')
    write(6,'("         -allow_sigma")')
    write(6,'("            By default, the computation for the sigma coordinate is not allowed in order to")')
    write(6,'("            avoid recurent errors done by users. For example if an eta file is passed to this")')
    write(6,'("            program without the record PT then a sigma coordinate will be faultly detected")')
    write(6,'("            leading to error in resulting values.")')
    write(6,'("            If one is certain that a sigma coordinate is present in the input file, then")')
    write(6,'("            this option may be used to unloock the safety mesure")')
    write(6,'("         -etiket ETIKET")')
    write(6,'("            To change the output record etiket. By default the output etiket will be:")')
    write(6,'("              The one of P0, or GZ id option -levels is not used.")')
    write(6,'("              The one of record nomvar if option -levels (nomvar) is used.")')
    write(6,'("         -verbose")')
    write(6,'("            Make program output more informations")')
    write(6,'("   ")')
  end subroutine cpg_print_doc
  !===================================================================================
  integer function cpg_process_arguments() result(status)
    use vGrid_Descriptors, only: vgd_putopt
    implicit none
    integer, parameter :: ncle=10, nmax=1000
    integer :: stat,npos
    integer :: fnom,fstouv,fstfrm
    character(len=256), dimension(ncle) :: cle,val,def
    status=VGD_ERROR
    cle=(/'s.         ','d.         ','samefile   ','levels     ','verbose    ','kind       ','version    ',&
         'allow_sigma','etiket     ','compute    '/)
    val=(/'undef      ','undef      ','undef      ','ALL_LEVELS ','undef      ','undef      ','undef      ',&
         'NO         ','undef      ','undef      '/)
    def=(/'undef      ','undef      ','YES        ','ALL_LEVELS ','INFO       ','undef      ','undef      ',&
         'YES        ','undef      ','undef      '/)
    npos=-111
    call ccard(cle,def,val,ncle,npos)
    stat=app_loglevel(trim(val(5))) 

    if(trim(val(10)).eq.'undef')then
       write(app_msg,*) 'Argument -compute (pressure,gz) is mandatory'
       call app_log(APP_ERROR,app_msg)
       call cpg_print_doc
       return
    endif
    if(trim(val(10)) == 'PRESSURE')then
       cpg_compute_gz_L=.false.
    elseif(trim(val(10)) == 'GZ')then
       cpg_compute_gz_L=.true.
    else
       write(6,'("Invalid value given to argument -compute, valid values are pressure or gz)")')
       return
    endif       
    if(trim(val(1)).eq.'undef')then
       call cpg_print_doc
       return
    endif
    cpg_samefile_L = trim(val(3)) == 'YES'
    if(cpg_samefile_L)then
       write(app_msg,*) 'Results will be written in input file '//trim(val(1))
       call app_log(APP_INFO,app_msg)
      endif
    if(trim(val(2)).eq.'undef' .and. .not. cpg_samefile_L)then
       write(app_msg,*) 'rgument -d (destination_file) is mandatory'
       call app_log(APP_ERROR,app_msg)
       call cpg_print_doc
       return
    endif
    cpg_levels_S=trim(val(4))
    if(trim(cpg_levels_S) == "ALL_LEVELS")then
      call app_log(APP_INFO,'Results will be computed on all levels')
    else
       write(app_msg,*) 'Results will be computed on '//trim(cpg_levels_S)//' levels'
       call app_log(APP_INFO,app_msg)
      endif
    if(trim(val(6)).eq.'undef')then
       cpg_kind=-1
    else
       read(val(6),*)cpg_kind
       write(app_msg,*) 'Looking for vertical descriptor of kind ',cpg_kind
       call app_log(APP_INFO,app_msg)
      endif
    if(trim(val(7)).eq.'undef')then
       cpg_version=-1
    else
       read(val(7),*)cpg_version
    endif
    if(cpg_kind.eq.-1.and.cpg_version.ne.-1)then
       call app_log(APP_ERROR,'argument -kind must be used with option -version')
         return
    endif
    if(cpg_version.ne.-1)then
       write(app_msg,*) 'Looking for vertical descriptor of version ',cpg_version
       call app_log(APP_INFO,app_msg)
      endif
    if(.not.cpg_samefile_L)then
       stat=fnom(cpg_lui,val(1),"RND+R/O",0)
       stat=fnom(cpg_luo,val(2),'RND',0)
       if(stat.lt.0)then
          write(app_msg,*) 'Problem with fnom on file ',trim(val(2))
          call app_log(APP_ERROR,app_msg)
          return
       endif
       stat=fstouv(cpg_luo,'RND')
       if(stat.lt.0)then
          write(app_msg,*) 'problem with fstouv on ',trim(val(2))
          call app_log(APP_ERROR,app_msg)
               stat=fstfrm(cpg_luo)
          return
       endif
    else
       stat=fnom(cpg_lui,val(1),"RND",0)
       cpg_luo=cpg_lui
    endif
    if(stat.lt.0)then
       call app_log(APP_ERROR,'ERROR with fnom on cpg_lui')
         return
    endif
    stat=fstouv(cpg_lui,'RND')
    if(stat.le.0)then
       call app_log(APP_WARNING,'No record in RPN file')
         return
    endif
    if(trim(val(8)).eq.'NO')then
       cpg_allow_sigma_L=.false.
    else if (trim(val(8)).eq.'YES')then
       cpg_allow_sigma_L=.true.
    else
       write(app_msg,*) "option -allow_sigma does't take any value, got ",trim(val(8))
       call app_log(APP_ERROR,app_msg)
         return
    endif
    stat = vgd_putopt("ALLOW_SIGMA",cpg_allow_sigma_L)
    if(val(9) == 'undef      ')then
       cpg_zap_etiket_L=.false.
    else
       cpg_zap_etiket_L=.true.
       cpg_etiket_S=trim(val(9))
       write(app_msg,*) 'Etiket of output records will be zap to ',cpg_etiket_S
       call app_log(APP_ERROR,app_msg)
    endif
    status=VGD_OK
  end function cpg_process_arguments
  !===================================================================================
  integer function cpg_print_prm(record) result(status)
    implicit none
    type(cpg_rpn) :: record
    ! Local variables
    status = VGD_ERROR
    print*,'dateo                =',record%dateo
    print*,'deet                 =',record%deet
    print*,'npas                 =',record%npas
    print*,'datev                =',record%datev
    print*,'ni,nj,nk             =',record%ni,record%nj,record%nk
    print*,'nbits                =',record%nbits
    print*,'datyp                =',record%datyp
    print*,'ip1,ip2,ip3          =',record%ip1,record%ip2,record%ip3
    print*,'typvat,nomvar        =',record%typvar,' ',record%nomvar
    print*,'etiket               =',record%etiket
    print*,'grtyp,ig1,ig2,ig3,ig4=',record%grtyp,record%ig1,record%ig2,record%ig3,record%ig4
    status = VGD_OK
  end function cpg_print_prm
  !===================================================================================
  integer function cpg_cp_params(F_po,F_pi,grid_only_L) result(status)
    implicit none
    type(cpg_rpn), intent(inout) :: F_po
    type(cpg_rpn), intent(in)  :: F_pi
    logical, intent(in), optional :: grid_only_L
    !
    ! Local variables
    !
    logical :: my_grid_only_L
    status = VGD_ERROR
    my_grid_only_L=.false.      
    if(present(grid_only_L))my_grid_only_L=grid_only_L
    !
    F_po%ni=F_pi%ni; F_po%nj=F_pi%nj; F_po%nk=F_pi%nk
    F_po%grtyp=F_pi%grtyp; F_po%ig1=F_pi%ig1; F_po%ig2=F_pi%ig2; F_po%ig3=F_pi%ig3; F_po%ig4=F_pi%ig4       
    if(.not.my_grid_only_L)then
       F_po%dateo=F_pi%dateo; F_po%deet=F_pi%deet; F_po%npas=F_pi%npas
       F_po%nbits=F_pi%nbits; F_po%datyp=F_pi%datyp
       F_po%ip1=F_pi%ip1; F_po%ip2=F_pi%ip2 ;F_po%ip3=F_pi%ip3
       F_po%typvar=F_pi%typvar;F_po%nomvar=F_pi%nomvar;F_po%etiket=F_pi%etiket
       F_po%swa=F_pi%swa; F_po%lng=F_pi%lng; F_po%dltf=F_pi%dltf; F_po%ubc=F_pi%ubc 
       F_po%extra1=F_pi%extra1; F_po%extra2=F_pi%extra2; F_po%extra3=F_pi%extra3
       F_po%datev=F_pi%datev
    endif
    !
    status = VGD_OK
  end function cpg_cp_params
  !===================================================================================
  integer function cpg_get_rec(F_f, F_lui, match_prm) result(status)
    type(cpg_rpn) :: F_f
    integer, intent(in) :: F_lui
    type(cpg_rpn), optional :: match_prm
    ! Local variables      
    integer, dimension(1000) :: keyList
    integer :: fstinl, ni, nj, nk, count, ier, fstluk
    type(cpg_rpn) :: prm
    status = VGD_ERROR
    if( present(match_prm) )then
       ier = cpg_cp_params(prm,match_prm)
    else
       prm%datev=-1; prm%etiket=' '; prm%ip1=-1;  prm%ip2=-1; prm%ip3=-1;
    endif
    if( fstinl(F_lui, ni, nj, nk, prm%datev, prm%etiket, prm%ip1, prm%ip2, prm%ip3, ' ',F_f%nomvar,&
         keyList,count,size(keyList)) < 0 )then
            call app_log(APP_ERROR,'ERROR in cpg_get_rec with fstinl')
       return
    endif
    if( count .eq. 0 )then
       write(app_msg,*) 'cpg_get_rec, no ', F_f%nomvar,' in input file with the following research key: datev = ', &
          prm%datev,', etiket = ',prm%etiket,',ip1 2 3 = ', prm%ip1, prm%ip2, prm%ip3
       call app_log(APP_ERROR,app_msg)
       return
    endif
    if( count > 1 )then
      write(app_msg,*) 'cpg_get_rec: more than one nomvar "'//F_f%nomvar//'" in inout file with the following research key:'
      call app_log(APP_ERROR,app_msg)
      write(app_msg,*) '      datev = ',prm%datev,', etiket = ',prm%etiket,',ip1 2 3 = ', prm%ip1, prm%ip2, prm%ip3
      call app_log(APP_VERBATIM,app_msg)
      write(app_msg,*) 'Use option -datev to select only one ', F_f%nomvar,' (Not implemented Yet, ask developer)'
      call app_log(APP_VERBATIM,app_msg)
       return
    endif
    if( associated(F_f%data) .and. &
         ( size(F_f%data,dim=1) /= ni .or. size(F_f%data,dim=2) /= nj ) )then
       deallocate(F_f%data)
    endif
    if( .not. associated(F_f%data) )then
       allocate(F_f%data(ni,nj),stat=ier)
       if(ier /= 0)then
          write(app_msg,*) 'cpg_get_rec, cannot alloacte F_f%data of size ->',ni,' x',nj
          call app_log(APP_ERROR,app_msg)
              return
       endif
    endif
    if( cpg_fstprm(keyList(1),F_f) == VGD_ERROR )then
       call app_log(APP_ERROR,'cpg_get_rec, with cpg_fstprm')
        return
    endif
    if( fstluk(F_f%data,keyList(1),ni,nj,nk) < 0 )then         
       write(app_msg,*) 'cpg_get_rec, with fstluk on ', F_f%nomvar
       call app_log(APP_ERROR,app_msg)
         return
    endif
    status = VGD_OK
  end function cpg_get_rec
  !===================================================================================
  integer function cpg_fstecr(F_lu,rec,F_rewrit,data) result(status)
    implicit none
    integer :: F_lu
    type(cpg_rpn) :: rec
    logical :: F_rewrit
    real, target, dimension(rec%ni,rec%nj), optional :: data
    ! Internal variables
    integer :: ier,fstecr
    real :: dumr
    real, dimension(:,:), pointer :: ptr
    status = VGD_ERROR      
    if(present(data))then
       ptr=>data
    else
       ptr=>rec%data
    endif
    ier = fstecr(ptr,dumr,-rec%nbits,F_lu,rec%dateo,rec%deet,rec%npas,rec%ni,rec%nj,&
         rec%nk,rec%ip1,rec%ip2,rec%ip3,rec%typvar,rec%nomvar,rec%etiket,rec%grtyp,&
         rec%ig1,rec%ig2,rec%ig3,rec%ig4,rec%datyp,F_rewrit)
    if(ier.ne.0)then
       write(6,*) 'ERROR: in cpg_fstecr, ier=',ier
       return
    endif
    status = VGD_OK
  end function cpg_fstecr
  !===================================================================================
  integer function cpg_sort_key_by_levels(sorted_keys,unsorted_keys,&
       sorted_levels,remove_duplicate_L,&
       sort_as_kind,only_kind,sorted_ip1s) result(status)
    implicit none
    integer, optional,intent(inout), pointer, dimension(:) :: sorted_keys, sorted_ip1s
    integer, optional,intent(in) , dimension(:) :: unsorted_keys      
    logical, optional,intent(in) :: remove_duplicate_L
    real, optional, pointer, dimension(:) ::sorted_levels
    integer, optional,intent(in) :: sort_as_kind,only_kind
    ! Local variables
    integer :: nkns,nks,ier,kind,i,k,itempo,iptr,l_sort_as_kind,l_only_kind,count,ig1,ig2
    integer, pointer, dimension(:) :: keys,ip1s
    real, pointer, dimension(:) :: hyb
    real :: tempo
    character (len=1) :: dummy_S
    logical :: l_remove_duplicate_L,reverse_L
    type (cpg_rpn) :: prm
    external :: convip
    nullify(prm%data,keys,hyb,ip1s)
    status=VGD_ERROR
    if(.not.present(sorted_keys))then
       write(app_msg,*) 'cpg_sort_key_by_levels: argument sorted_keys must be passed'
       call app_log(APP_ERROR,app_msg)
         return
    endif
    if(.not.present(unsorted_keys))then
       write(app_msg,*) 'cpg_sort_key_by_levels: argument unsorted_keys must be passed'
       call app_log(APP_ERROR,app_msg)
        return
    endif
    l_remove_duplicate_L=.false.
    if(present(remove_duplicate_L))l_remove_duplicate_L=remove_duplicate_L
    l_sort_as_kind=-1
    if(present(sort_as_kind))l_sort_as_kind=sort_as_kind
    l_only_kind=-1
    if(present(only_kind))l_only_kind=only_kind
    nkns=size(unsorted_keys)
    if(nkns == 0)then
       call app_log(APP_ERROR,'cpg_sort_key_by_levels: unsorted_keys list is empty')
         return
    endif
    write(app_msg,*) 'There are ',nkns,' unsorted levels'
    call app_log(APP_INFO,app_msg)
    allocate(keys(nkns),hyb(nkns),ip1s(nkns),stat=ier)
    if (ier /= 0)then
      call app_log(APP_ERROR,'cpg_sort_key_by_levels: allocation problem 1 in cpg_sort_key_by_levels')
       return         
    endif
    keys=unsorted_keys
    ! Sort levels by hyb values
    get_hyb: do k=1,nkns
       ier=cpg_fstprm(unsorted_keys(k),prm)
       if (ier /= 0)then
          write(app_msg,*) 'cpg_sort_key_by_levels: with cpg_fstprm on unsorted_keys k=',k
          call app_log(APP_ERROR,app_msg)
               return         
       endif
       call convip(prm%ip1,hyb(k),kind,-1,dummy_S,.false.)
       ip1s(k)=prm%ip1
       write(app_msg,*) 'Unsorted levels, k=',k,hyb(k),ip1s(k)
       call app_log(APP_DEBUG,app_msg)
       if(l_only_kind.ne.-1 .and. kind.ne.l_only_kind)then
         write(app_msg,*) '("Level",i10,f15.7," is not of kind ",i2," removing from list")',ip1s(k),hyb(k),l_only_kind
         call app_log(APP_DEBUG,app_msg)
          if(kind.ne.l_only_kind)then
             hyb(k)=-1
             ip1s(k)=-1
          endif
       endif
    enddo get_hyb
    if(l_only_kind.ne.-1)then
       count=0
       do k=1,nkns
          if(hyb(k).eq.-1)then
             ! look for a replacement
             do i=k+1,nkns
                if(hyb(i).ne.-1)then
                   hyb(k)=hyb(i)
                   ip1s(k)=ip1s(i)
                   hyb(i)=-1
                   ip1s(i)=-1
                   keys(k)=keys(i)
                   count=count+1
                   exit
                endif
             end do
          else
             count=count+1
          end if
          if(count.eq.0.and.k.eq.1)then
             write(app_msg,*) 'cpg_sort_key_by_levels: no level of kind ',l_only_kind,' found'
             call app_log(APP_WARNING,app_msg)
                     if (associated(sorted_keys))deallocate(sorted_keys)
             if(present(sorted_levels).and.associated(sorted_levels))deallocate(sorted_levels)
             if(present(sorted_ip1s).and.associated(sorted_ip1s))deallocate(sorted_ip1s)
             return
          endif
       end do
       nkns=count
    endif
    sort: do k=1,nkns-1
       iptr=k
       do i=k+1,nkns
          if(hyb(i).lt.hyb(iptr))iptr=i
       enddo
       swap: if(k.ne.iptr)then
          tempo=hyb(k)
          hyb(k)=hyb(iptr)            
          hyb(iptr)=tempo
          itempo=ip1s(k)
          ip1s(k)=ip1s(iptr)            
          ip1s(iptr)=itempo
          itempo=keys(k)
          keys(k)=keys(iptr)
          keys(iptr)=itempo      
       endif swap
    enddo sort
    ! Identify duplicate by -1 in hyb and count distinct levels
    dup: if(.not.l_remove_duplicate_L)then
       nks=nkns
    else
       nks=1
       do k=1,nkns-1
          if(hyb(k).eq.-1)cycle
          do i=k+1,nkns
             if(hyb(i).eq.hyb(k))then
                hyb(i)=-1
                ip1s(i)=-1
             else
                nks=nks+1
                exit
             endif
          enddo
       enddo
    endif dup
    if (associated(sorted_keys))deallocate(sorted_keys)
    allocate(sorted_keys(nks),stat=ier)
    if (ier /= 0)then
       write(app_msg,*) 'cpg_sort_key_by_levels: allocation problem 2 in cpg_sort_key_by_levels'
       call app_log(APP_ERROR,app_msg)
         return         
    endif
    if(present(sorted_levels))then
       if (associated(sorted_levels))deallocate(sorted_levels)
       allocate(sorted_levels(nks),stat=ier)
       if (ier /= 0)then
          call app_log(APP_ERROR,'ERROR in cpg_sort_key_by_levels: allocation problem 3 in cpg_sort_key_by_levels')
          return         
       endif
    endif
    if(present(sorted_ip1s))then
       if (associated(sorted_ip1s))deallocate(sorted_ip1s)
       allocate(sorted_ip1s(nks),stat=ier)
       if (ier /= 0)then
          call app_log(APP_ERROR,'ERROR in cpg_sort_key_by_levels: allocation problem 4 in cpg_sort_key_by_levels')
               return         
       endif
    endif
    i=1
    do k=1,nkns
       if(hyb(k).ne.-1)then
          sorted_keys(i)=keys(k)
          if(present(sorted_levels))sorted_levels(i)=hyb(k)
          if(present(sorted_ip1s))sorted_ip1s(i)=ip1s(k)
          i=i+1
       endif
    enddo
    !
    ! Reverse order for certain kind
    ! 0       : m  [metres] (height with respect to sea level)
    ! 1       : sg [sigma] (0.0->1.0)
    ! 2       : mb [mbars] (pressure in millibars)
    ! 3       :    [others] (arbitrary code)
    ! 4       : M  [metres] (height with respect to ground level)
    ! 5       : hy [hybrid] (0.0->1.0)
    ! 6       : th [theta]
    if(l_sort_as_kind.ne.-1)kind=l_sort_as_kind
    select case (kind)
    case (1,2,3,5) 
       reverse_L=.false.
    case (0,4,6)       
       reverse_L=.true.
       keys(1:size(sorted_keys)) = sorted_keys(size(sorted_keys):1:-1)
       sorted_keys=keys(1:size(sorted_keys))         
       if(present(sorted_levels))sorted_levels(1:size(sorted_keys))=hyb(size(sorted_keys):1:-1)         
       if(present(sorted_ip1s))sorted_ip1s(1:size(sorted_keys))=ip1s(size(sorted_keys):1:-1)   
    case DEFAULT
       write(app_msg,*) 'cpg_sort_key_by_levels: unsupported kind=',kind
       call app_log(APP_ERROR,app_msg)
         return    
    end select
    !
    if(app_loglevel('').gt.1)then
       write(app_msg,*) 'There are ',nks,' sorted levels'
       call app_log(APP_INFO,app_msg)
       if(l_sort_as_kind.ne.-1) then
         write(app_msg,*) 'Levels were sorted as kind =',l_sort_as_kind
         call app_log(APP_INFO,app_msg)
       endif
       if(l_remove_duplicate_L) then
         write(app_msg,*) 'There was',nkns-nks,' duplicated levels'
         call app_log(APP_ERROR,app_msg)
       endif
       do k=1,nks
          ier=cpg_fstprm(sorted_keys(k),prm)
          if (ier /= 0)then
             return         
          endif
          call convip(prm%ip1,hyb(k),kind,-1,dummy_S,.false.)      
          write(app_msg,*) 'Sorted levels, k=',k,hyb(k),ip1s(k)
          call app_log(APP_INFO,app_msg)
       enddo
    endif
    deallocate(keys,hyb,ip1s)
    status=VGD_OK
  end function cpg_sort_key_by_levels
  !===================================================================================
  integer function cpg_get_px_gz() result(stat)
    use app
    use vgrid_descriptors, only: vgd_get
    implicit none
    ! Local variables
    integer :: vcode
    stat=VGD_ERROR
    if(vgd_get(cpg_vgd,"VCOD",vcode) == VGD_ERROR)then
       return
    endif
    select case ( vcode )
    case (1001,1002,1003,5001)
       ! Computation of gz like it was done in GEM 3.3.3 with p0vt2gz_hyb
       ! It doesn't give better results than simple average below.
       if( cpg_get_gz_1002_5001_simple_avg(vcode) == VGD_ERROR )then
          stat=app_end(-1)
          error stop
       endif
    case (5002)
       if( cpg_get_gz_5002_5005_5100(vcode) == VGD_ERROR )then
          stat=app_end(-1)
          error stop
       endif
    case (5005)     
       if( cpg_get_gz_5002_5005_5100(vcode) == VGD_ERROR )then
          stat=app_end(-1)
          error stop
       endif
    case (5100)     
       if( cpg_get_gz_5002_5005_5100(vcode) == VGD_ERROR )then
          stat=app_end(-1)
          error stop
       endif
    case (21001)
       if( cpg_get_px_21001() == VGD_ERROR )then
          stat=app_end(-1)
          error stop
       endif
    case DEFAULT
      write(app_msg,*) 'Vcode ',vcode,' not supported, please contact developer to add this'
      call app_log(APP_ERROR,app_msg)
       stat=app_end(-1)
       error stop
    end select
    stat=VGD_OK
  end function cpg_get_px_gz
  !===================================================================================  
  integer function cpg_get_gz_1002_5001_simple_avg(F_vcode) result(status)
    use vgrid_descriptors, only: vgrid_descriptor, vgd_new, vgd_get, vgd_levels
    use tdpack, only: grav_8, rgasd_8, tcdk, delta
    implicit none
    integer, intent(in) :: F_vcode
    ! local variables
    integer, parameter :: nmax=10000
    integer :: nij, ier, nk, i,j,k, infon, fstinl,ig1,ig2,kind
    integer, dimension(:), pointer :: ip1s_m ,sorted_keys     
    integer, dimension(nmax) :: unsorted_keys
    real, dimension(:,:), pointer :: w2a,w2b
    real, dimension(:,:,:), pointer :: w3a,w3b
    real, dimension(:,:), pointer :: ptr2_t1, ptr2_t2, ptr2_tmp
    real, dimension(:,:,:), pointer :: ptr3_p1, ptr3_p2, ptr3_tmp
    real :: aaa
    logical :: vt_L
    type(cpg_rpn) :: p0,gz,tt,hu,prm 
    external :: mfotvt
    nullify(w2a,w2b,w3a,w3b,ptr3_p1,ptr3_p2,ptr2_t1,ptr2_t2,ip1s_m, &
         p0%data,gz%data,tt%data,hu%data,prm%data,sorted_keys)
    status = VGD_ERROR
    p0%nomvar="p0"
    if( cpg_get_rec(p0, cpg_lui) == VGD_ERROR ) return      
    p0%data=p0%data*100.
    ig1 = p0%ig1; ig2 = p0%ig2;
    nij =  p0%ni*p0%nj
    allocate(w3a(p0%ni,p0%nj,1),w3b(p0%ni,p0%nj,1),&
         w2a(p0%ni,p0%nj),w2b(p0%ni,p0%nj),prm%data(p0%ni,p0%nj),stat=ier)
    if( ier /= 0 )then
      write(app_msg,*) 'cpg_get_gz_1002_5001_simple_avg: cannot allocate w* of size ->',prm%ni,' x',prm%nj, 'x 1'
      call app_log(APP_ERROR,app_msg)
      return
    end if
    ! get all VT, or TT ip1s
    ier = fstinl(cpg_lui,tt%ni,tt%ni,tt%nk,p0%datev,p0%etiket,-1,p0%ip2,-1,'','VT',&
         unsorted_keys, infon, nmax)
    if(ier < 0 )then
      write(app_msg,*) 'cpg_get_gz_1002_5001_simple_avg: fstinl on VT'
      call app_log(APP_ERROR,app_msg)
       return
    endif
    vt_L = .false.
    if( infon > 0 )then
       vt_L = .true.
    else
       ier = fstinl(cpg_lui,tt%ni,tt%ni,tt%nk,p0%datev,p0%etiket,-1,p0%ip2,-1,'','TT',&
            unsorted_keys, infon, nmax)
       if(ier < 0 )then
         write(app_msg,*) 'cpg_get_gz_1002_5001_simple_avg: fstinl on TT'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       if( infon == 0 )then
         write(app_msg,*) 'There is no VT or TT record in the input file'
         call app_log(APP_ERROR,app_msg)
          return
       endif
    endif
    ier = vgd_get(cpg_vgd,"KIND",kind)
    if(cpg_sort_key_by_levels(sorted_keys,unsorted_keys(1:infon),&
         remove_duplicate_L=.true.,&
         only_kind=kind,sorted_ip1s=ip1s_m) == VGD_ERROR)return
    write(app_msg,*) 'List of momentum ip1s',ip1s_m
    call app_log(APP_INFO,app_msg)

    nk = size(ip1s_m)
    write(app_msg,*) 'computing 3d pressure'
    call app_log(APP_INFO,app_msg)

    gz%nomvar="gz"
    ier = cpg_cp_params(prm,p0)
    prm%ip1=12000
    if(kind == 5)prm%ip1=93423264
    if( cpg_get_rec(gz, cpg_lui, match_prm=prm) == VGD_ERROR ) return
    ier = cpg_cp_params(prm,gz)
    if(cpg_zap_etiket_L)prm%etiket=cpg_etiket_S
    if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=gz%data) == VGD_ERROR)return
    gz%data=gz%data*10.
    ! start integration from surface up
    aaa=rgasd_8/grav_8
    w3a(:,:,1)=log(p0%data(:,:))
    ptr3_p1 => w3a
    ptr3_p2 => w3b
    ier = cpg_cp_params(tt,gz)
    if(vt_L)then
       tt%nomvar="vt"
    else
       tt%nomvar="tt"
       ier = cpg_cp_params(hu,gz)
       hu%nomvar="hu"
    endif
    if( cpg_get_rec(tt, cpg_lui, match_prm=tt) == VGD_ERROR ) return
    ptr2_t1 => w2a
    ptr2_t2 => w2b
    if(vt_L)then
       ptr2_t1 = tt%data + tcdk
    else
       tt%data = tt%data + tcdk
       if( cpg_get_rec(hu, cpg_lui, match_prm=hu) == VGD_ERROR ) return
       call mfotvt(ptr2_t1,tt%data,hu%data,nij,1,nij)
    endif
    ier = cpg_cp_params(prm,p0)
    prm%nomvar='GZ'
    do k=nk-1,1,-1
       tt%ip1=ip1s_m(k);
       if( cpg_get_rec(tt, cpg_lui, match_prm=tt) == VGD_ERROR ) return
       if(vt_L)then
          ptr2_t2 = tt%data + tcdk
       else
          tt%data = tt%data + tcdk
          hu%ip1=ip1s_m(k);
          if( cpg_get_rec(hu, cpg_lui, match_prm=hu) == VGD_ERROR ) return
          call mfotvt(ptr2_t2,tt%data,hu%data,nij,1,nij)
       endif
       if( vgd_levels(cpg_vgd,sfc_field=p0%data,ip1_list=ip1s_m(k:k),levels=ptr3_p2,in_log=.true.) == VGD_ERROR)return
       do j=1,p0%nj
          do i=1,p0%ni
             gz%data(i,j) = gz%data(i,j) &
                  - aaa * .5*(ptr2_t2(i,j)+ptr2_t1(i,j)) &
                  * (ptr3_p2(i,j,1) - ptr3_p1(i,j,1))
          end do
       end do
       prm%ip1=ip1s_m(k);
       prm%data=gz%data/10.
       if(cpg_zap_etiket_L)prm%etiket=cpg_etiket_S
       if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=prm%data) == VGD_ERROR)return
       ! Swap pointers
       ptr2_tmp => ptr2_t1
       ptr2_t1 => ptr2_t2
       ptr2_t2 => ptr2_tmp
       ptr3_tmp => ptr3_p1
       ptr3_p1 => ptr3_p2
       ptr3_p2 => ptr3_tmp
    end do
    status = VGD_OK   
  end function cpg_get_gz_1002_5001_simple_avg
  !===================================================================================  
  integer function cpg_get_gz_5002_5005_5100(F_vcode) result(status)
    use vgrid_descriptors, only: vgrid_descriptor, vgd_new, vgd_get, &
         vgd_levels
    use, intrinsic :: iso_fortran_env, only: REAL64
    use tdpack, only: rgasd,grav,tcdk,delta
    implicit none
    integer, intent(in) :: F_vcode
    ! local variables
    type(cpg_rpn) :: p0,p0ls,gz,prm,tt,hu,vt
    real, dimension(:,:), pointer :: gzt
    real, dimension(:,:,:), pointer :: wa,wb
    real, dimension(:,:,:), pointer :: ptr_p1, ptr_p2, ptr_tempo
    real :: work, inc
    real(kind=REAL64) :: aaa
    integer, dimension(:), pointer :: ip1s_m, ip1s_t
    integer :: i, j, ier, k, nk, nij, kp, ig1, ig2,vcode
    logical :: thermo_L, momentum_L, vt_L
    external :: mfotvt
    nullify(gzt, wa, wb, ptr_p1, ptr_p2, ptr_tempo, ip1s_m, ip1s_t, p0%data, p0ls%data, gz%data, prm%data, tt%data, hu%data, vt%data)
    status = VGD_ERROR
    if(vgd_get(cpg_vgd,"VCOD",vcode) == VGD_ERROR)then
       return
    endif
    thermo_L=.false.
    momentum_L=.false.
    if(trim(cpg_levels_S) == 'THERMO')then
       thermo_L=.true.
    elseif(trim(cpg_levels_S) == 'MOMENTUM')then
       momentum_L=.true.
    elseif(trim(cpg_levels_S) == 'ALL_LEVELS')then
       thermo_L=.true.
       momentum_L=.true.
    else
      write(app_msg,*) 'ERROR with option -levels expected THERMO, MOMENTUM or ALL_LEVELS but got ',trim(cpg_levels_S)
      call app_log(APP_ERROR,app_msg)
      write(app_msg,*) 'NOTE: passing a nomvar to option -levels is not implemented for vcode ',vcode
      call app_log(APP_VERBATIM,app_msg)
       return
    endif

    ! get all momentum ip1s
    if( vgd_get(cpg_vgd,"vipm - level ip1 list (m)", ip1s_m) == vgd_error)return
    write(app_msg,*) 'List of momentum ip1s',ip1s_m
    call app_log(APP_INFO,app_msg)

    ! get all thermo ip1s
    if( vgd_get(cpg_vgd,"vipt - level ip1 list (t)", ip1s_t) == vgd_error)return
    write(app_msg,*) 'List of thermo ip1s_t',ip1s_t
    call app_log(APP_INFO,app_msg)

    ! For vcode 5002, there are nk+t dynamic temperature, there are nk+2 ip1s_t counting 
    !                 the diag level. This last level is not
    !                 used for computing gz.
    !                 There are nk momentum levels, there are nk+1 ip1s_m counting
    !                 the diag level.
    ! For vcode 5005
    ! and vcode 5100, there are nk dynamic temperature, 
    !                 there are nk+2 ip1s_t counting hyb=1.0 and diag levels, these 
    !                 last 2 levels are not used for computing gz.
    if( F_vcode == 5002 )then
       nk = size(ip1s_m)-1
       kp = 1
    else if(F_vcode == 5005 .or. F_vcode == 5100 )then
       nk = size(ip1s_m)-2
       kp = 0
    else
      write(app_msg,*) 'invalid Vcode in call to get_gz_5002_5005_5100, got ',F_vcode
      call app_log(APP_ERROR,app_msg)
       return
    endif
    p0%nomvar="p0"
    if( cpg_get_rec(p0, cpg_lui) == VGD_ERROR ) return
    p0%data=p0%data*100.
    if(F_vcode == 5100)then
       p0ls%nomvar="p0ls"
       if( cpg_get_rec(p0ls, cpg_lui) == VGD_ERROR ) return
       p0ls%data=p0ls%data*100.
    endif
    ig1 = p0%ig1; ig2 = p0%ig2;
    nij =  p0%ni*p0%nj    
    allocate(wa(p0%ni,p0%nj,1),wb(p0%ni,p0%nj,1),gzt(p0%ni,p0%nj),stat=ier)
    if( ier /= 0 )then
      write(app_msg,*) 'error in get_gz_5002_5005_5100, cannot allocate wa, wb and gztof size ->',prm%ni,' x',prm%nj
      call app_log(APP_ERROR,app_msg)
    end if
    write(app_msg,*) 'computing 3d pressure'
    call app_log(APP_INFO,app_msg)

    gz%nomvar="gz"
    ier = cpg_cp_params(prm,p0)    
    prm%ip1=93423264
    !ier=cpg_print_prm(prm)
    if( cpg_get_rec(gz, cpg_lui, match_prm=prm) == VGD_ERROR ) return
    ier = cpg_cp_params(prm,gz)
    if( cpg_fstecr(cpg_luo, gz, F_rewrit=.true., data=gz%data) == VGD_ERROR)return
    gz%data=gz%data*10.
    ! start integration from surface up
    aaa=rgasd/grav
    wa(:,:,1)=log(p0%data(:,:))
    ptr_p1 => wa
    ptr_p2 => wb
    vt%nomvar="vt"
    tt%nomvar="tt"
    hu%nomvar="hu"

    ! Try to get one VT record to set vt_L flag
    vt_L = .true.
    ier = cpg_cp_params(prm,p0)
    prm%ip1=ip1s_t(1+kp);
    if( cpg_get_rec(vt, cpg_lui, match_prm=prm) == VGD_ERROR )then
       vt_L = .false.
    endif
    
    do k=nk,1,-1
       ier = cpg_cp_params(prm,p0)
       prm%ip1=ip1s_t(k+kp);
       if( vt_L )then
          if( cpg_get_rec(vt, cpg_lui, match_prm=prm) == VGD_ERROR )return
          vt%data(:,:) = vt%data(:,:) + tcdk
          ier = cpg_cp_params(prm,vt)
       else
          if( cpg_get_rec(tt, cpg_lui, match_prm=prm) == VGD_ERROR )return
          tt%data(:,:) = tt%data(:,:) + tcdk
          if( cpg_get_rec(hu, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          if(.not.associated(vt%data))then
             allocate(vt%data(prm%ni,prm%nj),stat=ier)
             if( ier /= 0 )then
               write(app_msg,*) 'get_gz_5002_5005: cannot allocate vt%data of size ->',prm%ni,' x',prm%nj
               call app_log(APP_ERROR,app_msg)
               return
             end if
          endif
          call mfotvt(vt%data,tt%data,hu%data,nij,1,nij)
          ier = cpg_cp_params(prm,tt)
       endif
       if(F_vcode == 5100)then
          if( vgd_levels(cpg_vgd,sfc_field=p0%data,sfc_field_ls=p0ls%data,ip1_list=ip1s_m(k:k),levels=ptr_p2,in_log=.true.) &
               == VGD_ERROR)return
       else
          if( vgd_levels(cpg_vgd,sfc_field=p0%data,ip1_list=ip1s_m(k:k),levels=ptr_p2,in_log=.true.) == VGD_ERROR)return
       endif
       do j=1,p0%nj
          do i=1,p0%ni
             inc = - aaa * vt%data(i,j) * ( ptr_p2(i,j,1) - ptr_p1(i,j,1) )
             ! GZ thermo is at the center of the layer for Vcode 5005
             gzt(i,j) = (gz%data(i,j) + .5*inc)/10.             
             gz%data(i,j) = (gz%data(i,j) + inc)/10.
          end do
       end do
       if( cpg_zap_etiket_L ) prm%etiket = cpg_etiket_S
       prm%nomvar='GZ'
       if(thermo_L)then
          if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=gzt) == VGD_ERROR)return
       endif
       prm%ip1=ip1s_m(k);
       if(momentum_L)then
          if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=gz%data) == VGD_ERROR)return
       endif
       gz%data=gz%data*10.
       ptr_tempo => ptr_p1
       ptr_p1 => ptr_p2
       ptr_p2 => ptr_tempo
    end do
    if( F_Vcode == 5002 .and. thermo_L)then
       ! Additional thermo level above the firt momentum level.
       prm%ip1=ip1s_t(1);
       if( vt_L )then
          if( cpg_get_rec(vt, cpg_lui, match_prm=prm) == VGD_ERROR )return
          vt%data(:,:) = vt%data(:,:) + tcdk
          ier = cpg_cp_params(prm,tt)
       else
          if( cpg_get_rec(tt, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          tt%data(:,:) = tt%data(:,:) + tcdk
          if( cpg_get_rec(hu, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          call mfotvt(vt%data,tt%data,hu%data,nij,1,nij)
          ier = cpg_cp_params(prm,tt)
       endif
       if( vgd_levels(cpg_vgd,sfc_field=p0%data,ip1_list=ip1s_t(1:1),levels=ptr_p2,in_log=.true.) == VGD_ERROR)return
       do j=1,p0%nj
          do i=1,p0%ni
             gzt(i,j) = (gz%data(i,j) - aaa * vt%data(i,j) * ( ptr_p2(i,j,1) - ptr_p1(i,j,1)))/10.
          end do
       end do
       if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=gzt) == VGD_ERROR)return
    endif
    deallocate(wa, wb, gzt, ip1s_m, ip1s_t, vt%data, p0%data, gz%data)
    if(associated(p0ls%data))deallocate(p0ls%data)
    if(associated(tt%data))deallocate(tt%data)
    if(associated(hu%data))deallocate(hu%data)
    status = VGD_OK
  end function cpg_get_gz_5002_5005_5100
  !===================================================================================
  integer function cpg_get_vgd_levels() result(stat)

    use vGrid_Descriptors, only: vgd_new,vgd_levels,vgd_get, &
           vgrid_descriptor,vgd_print,vgd_write,&
           VGD_LEN_RFLD, VGD_LEN_RFLS, VGD_NO_REF_NOMVAR
    implicit none

    !Local variable
    
    type(cpg_rpn) :: record,record2
    integer, parameter :: nmax=1000
    integer, dimension(nmax) :: liste,liste2
    integer :: nk2,nk3,ier,fstinl,fstluk,fstecr,fstinf,ni,nj,nk,infon,n_ip1,i,j,key,kind,version
    integer, dimension(:), pointer :: ip1_list,ip1_list2,ip1_list3,ip1_single
    character(len=VGD_LEN_RFLD) :: nomvar
    character(len=VGD_LEN_RFLS) :: nomvar_ls
    real :: dummy
    real, dimension(:,:,:), pointer :: pres
    real, dimension(:,:), pointer :: work,p0ls
    real :: pp
    logical :: OK_2001_L, final_ip1_list_L
    
    stat=VGD_ERROR

    nullify(ip1_list,ip1_list2,ip1_list3,ip1_single,pres,work,p0ls,record%data)
    if(vgd_get(cpg_vgd,'KIND - vertical coordinate ip1 kind', kind) == VGD_ERROR)then
       return
    endif
    allocate(ip1_single(1))
    OK_2001_L=.false.
    final_ip1_list_L=.true.
    if(trim(cpg_levels_S).eq.'THERMO')then
       ier=vgd_get(cpg_vgd,'VIPT - level ip1 list (t)',value=ip1_list)
       if(ier.ne.VGD_OK)then
         write(app_msg,*) 'vgd_get on VIPT for var (group) ',cpg_levels_S
         call app_log(APP_ERROR,app_msg)
          return
       endif
    else if(trim(cpg_levels_S).eq.'MOMENTUM')then
       ier=vgd_get(cpg_vgd,'VIPM - level ip1 list (m)',value=ip1_list)
       if(ier.ne.VGD_OK)then
         write(app_msg,*) 'vgd_get on VIPM for var (group) ',cpg_levels_S
         call app_log(APP_ERROR,app_msg)
          return
       endif
    else if(trim(cpg_levels_S).eq.'VERTICAL_VELOCITY')then
       ier=vgd_get(cpg_vgd,'VIPW - level ip1 list (w)',value=ip1_list)
       if(ier.ne.VGD_OK)then
         write(app_msg,*) 'vgd_get on VIPW for var (group) ',cpg_levels_S
         call app_log(APP_ERROR,app_msg)
          return
       endif
    else if (trim(cpg_levels_S).eq.'ALL_LEVELS')then
       ier=vgd_get(cpg_vgd,'VIPT - level ip1 list (t)',value=ip1_list2)
       if(ier.ne.VGD_OK)then
         write(app_msg,*) 'vgd_get on VIPT for var (group) ',cpg_levels_S
         call app_log(APP_ERROR,app_msg)
          return
       endif
       ier=vgd_get(cpg_vgd,'VIPM - level ip1 list (m)',value=ip1_list3)
       if(ier.ne.VGD_OK)then
         write(app_msg,*) 'vgd_get on VIPM for var (group) ',cpg_levels_S
         call app_log(APP_ERROR,app_msg)
          return
       endif
       !print*,"TODO vgd_get(cpg_vgd,'VIPW ..."
       !error stop
       nk2=size(ip1_list2)
       nk3=size(ip1_list3)
       allocate(ip1_list(nk2+nk3))
       ip1_list(1:nk2)=ip1_list2
       ip1_list(nk2+1:nk2+nk3)=ip1_list3
    else
       OK_2001_L=.true.
       final_ip1_list_L=.false.
       ! If there is more than one P0 in file, the ip1_list will contain all 
       !    ip1 for nomvar cpg_levels_S and there will be repetitions.
       !    The ip1_list will have to be rebuild once the datev of P0 is known in P0 loop below.
       ier=fstinl(cpg_lui,ni,nj,nk,-1,' ',-1,-1,-1,' ',cpg_levels_S,liste,n_ip1,nmax)           
       if(ier.lt.0)then
          return
       endif
       if(n_ip1.eq.0)then
         write(app_msg,*) 'No record of nombav '//trim(cpg_levels_S)//' in input file'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       allocate(ip1_list(n_ip1),stat=ier)
       if(ier.ne.0)then
         write(app_msg,*) 'Problem in allocate ip1_list(n_ip1)'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       do i=1,n_ip1
          if(cpg_fstprm(liste(i),record) == VGD_ERROR)then
                return
          endif
          ip1_list(i)=record%ip1
       enddo
       kind = cpg_kind
       if ( kind == -1 )then
          ! Loop on ip1_list to find kind other then 4 (diag level)
          do i=1,n_ip1
             call convip(ip1_list(i), pp, kind, -1,"",.false.)
             if(kind.ne.4)exit
          enddo
       endif
    endif
    n_ip1=size(ip1_list)
    ier=vgd_get(cpg_vgd,'KIND - vertical coordinate ip1 kind', kind)
    if(ier.ne.VGD_OK)then
       return
    endif
    ier=vgd_get(cpg_vgd,'VERS - vertical coordinate version' , version)
    if(ier.ne.VGD_OK)then
       return
    endif
    if(kind*1000+version == 2001) then
       if(.not. OK_2001_L)then
         write(app_msg,*) 'levels trim(cpg_levels_S) not allowed with pressure level file. &
      & Use a nomvar present in the input file e.g. -levels TT'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       ! Get the list of field at different time, this will be use has P0 but 
       ! will never be used in computation since 2001 is pressure levels         
         
       key=fstinf(cpg_lui,ni,nj,nk,-1,' ',record%ip1,-1,-1,' ',cpg_levels_S)
       if(key.lt.0)then
         write(app_msg,*) 'cannot find variable ',nomvar,' in input file'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       nomvar=cpg_levels_S
    else
       ier=vgd_get(cpg_vgd,'RFLD - reference field name',value=nomvar)
       if(ier.ne.VGD_OK)then
          return
       endif
       ! Use wild card to find RFLD
       record%ip1=-1
    endif
    ier=fstinl(cpg_lui,ni,nj,nk,-1,' ',record%ip1,-1,-1,' ',nomvar,liste,infon,nmax)
    if(ier.lt.0)then
      write(app_msg,*) 'fstinl on ',nomvar
      call app_log(APP_ERROR,app_msg)
       return
    endif
    if(infon.eq.0)then
      write(app_msg,*) 'No record ',nomvar,' in input file'
      call app_log(APP_ERROR,app_msg)
       return
    endif
    ! Loop on all p0
    LOOP_ON_P0: do i=1,infon
       if(cpg_fstprm(liste(i),record) == VGD_ERROR)then
          return
       endif       
       ier=vgd_write(cpg_vgd,cpg_luo,'fst')
       if(record%ni.ne.ni.or. &
            record%nj.ne.nj.or. &
            record%nk.ne.nk)then
         write(app_msg,*) 'Size of record ',nomvar,' inconsistant with previous one (ni:',record%ni, &
            'vs',ni,'nj:',record%nj,'vs',nj,'nk:',record%nk,'vs',nk,')'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       if(associated(record%data))deallocate(record%data)
       allocate(record%data(ni,nj),stat=ier)
       if(ier.ne.0)then
         write(app_msg,*) 'Problem in allocate record%data'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       if(associated(work))deallocate(work)
       allocate(work(ni,nj),stat=ier)
       if(ier.ne.0)then
         write(app_msg,*) 'Problem in allocate work'
         call app_log(APP_ERROR,app_msg)
          return
       endif
       ier=fstluk(record%data,liste(i),ni,nj,nk)
       if(ier.lt.0)then
          return
       endif
       work=record%data
       if(trim(nomvar) == "P0")record%data=record%data*100.       
       ier=vgd_get(cpg_vgd,'RFLS - large scale reference field name',value=nomvar_ls)         
       if(nomvar_ls == VGD_NO_REF_NOMVAR)then
          if(associated(p0ls))deallocate(p0ls)
       else
          key = fstinf(cpg_lui,ni,nj,nk,record%datev,record%etiket,record%ip1,record%ip2,-1,' ',nomvar_ls)
          if(key < 0)then
            write(app_msg,*) 'cannot find large scale reference field ',nomvar_ls,' with the following parameters'
            call app_log(APP_ERROR,app_msg)
             write(6,'("datev: ",i12,", Etiket: ",a,", ip1: ",i4,", ip2: ",i4)') &
                  record%datev,record%etiket,record%ip1,record%ip2
             return
          endif
          if(associated(p0ls))deallocate(p0ls)
          allocate(p0ls(ni,nj),stat=ier)
          if(ier.ne.0)then
            write(app_msg,*) 'Problem in allocate p0ls'
            call app_log(APP_ERROR,app_msg)
             return
          endif
          ier=fstluk(p0ls,key,ni,nj,nk)
          if(ier.lt.0)then
            write(app_msg,*) 'Problem with fstluk on ',nomvar_ls
            call app_log(APP_ERROR,app_msg)
             return
          endif
          if(trim(nomvar_ls) == "P0LS")p0ls=p0ls*100.
       endif       
       ! Loop on ip1_list was added since memory allocation for the cube was too large with
       ! some global grids.               
       if(.not. final_ip1_list_L)then
          ! Get ip1_list for current P0 datev
          ier=fstinl(cpg_lui,ni,nj,nk,record%datev,' ',-1,-1,-1,' ',cpg_levels_S,liste2,n_ip1,nmax)           
          if(ier.lt.0)then
            write(app_msg,*) 'Problem with fstinl 2 on ',cpg_levels_S,' datev = ',record%datev
            call app_log(APP_ERROR,app_msg)
             return
          endif
          if(n_ip1.eq.0)then
            write(app_msg,*) 'No record of nomvar ',trim(cpg_levels_S),' datev = ',record%datev,' in input file'
            call app_log(APP_ERROR,app_msg)
             return
          endif
          ! Note we got a subset of ip1s, therefore ip1_list is long enough.
          do j=1,n_ip1
             if(cpg_fstprm(liste2(j),record2) == VGD_ERROR)then
                return
             endif
             ip1_list(j)=record2%ip1
          enddo
       endif       
       LOOP_ON_IP1_LIST: do j=1,n_ip1
          ip1_single=ip1_list(j)
          if(associated(p0ls))then
             ier=vgd_levels(cpg_vgd,sfc_field=record%data,sfc_field_ls=P0ls,ip1_list=ip1_single,levels=pres)            
          else
             ier=vgd_levels(cpg_vgd,sfc_field=record%data,ip1_list=ip1_single,levels=pres)            
          endif
          if(ier.ne.VGD_OK)then
             return
          endif
          if(cpg_is_pressure_L)then
             pres=pres/100.
             record%nomvar="PX"
          else
             pres=pres/10.
             record%nomvar="GZ"
          endif
          if(cpg_zap_etiket_L)record%etiket=cpg_etiket_s
          ier=fstecr(pres,dummy,-record%nbits,cpg_luo,record%dateo,record%deet,record%npas, &
               record%ni,record%nj,record%nk,ip1_single,record%ip2,record%ip3, &
               record%typvar,record%nomvar,record%etiket,record%grtyp, &
               record%ig1,record%ig2,record%ig3,record%ig4,record%datyp,.true.)
          if(ier.lt.0)then
            write(app_msg,*) 'Problem with fstecr for var (group) ',cpg_levels_S
            call app_log(APP_ERROR,app_msg)
             return
          endif
       enddo LOOP_ON_IP1_LIST
    enddo LOOP_ON_P0
    deallocate(record%data,work,ip1_single)
    if(associated(pres )    )deallocate(pres)
    if(associated(ip1_list) )deallocate(ip1_list)
    if(associated(ip1_list2))deallocate(ip1_list2)
    if(associated(ip1_list3))deallocate(ip1_list3)
    if(associated(work)     )deallocate(work)
    if(associated(p0ls)     )deallocate(p0ls)
    stat=VGD_OK
    return
  end function cpg_get_vgd_levels
  !===================================================================================
  integer function cpg_get_px_21001() result(status)
    use vgrid_descriptors, only: vgrid_descriptor, vgd_new, vgd_get,vgd_levels,&
         VGD_LEN_RFLD,VGD_NO_REF_NOMVAR
    use, intrinsic :: iso_fortran_env, only: REAL64
    use tdpack, only: rgasd,grav,tcdk,delta
    implicit none
    ! local variables
    type(cpg_rpn) :: me,mels,px,gz,prm,tt,hu,vt
    real, dimension(:,:), pointer :: pxt
    real, dimension(:,:,:), pointer :: wa,wb
    real, dimension(:,:,:), pointer :: ptr_p1, ptr_p2, ptr_tempo
    real :: work, inc
    real(kind=REAL64) :: aaa,fac_8
    integer, dimension(:), pointer :: ip1s_m, ip1s_t
    integer :: i, j, ier, k, ni,nj,nk, nij,vcode
    integer, parameter :: ip1_sfc=93423364
    logical :: thermo_L, momentum_L, sleve_L, vt_L
    character(len=VGD_LEN_RFLD) :: refls_name
    integer, external :: fstinf
    external :: mfotvt
    nullify(me%data,mels%data, pxt, wa,wb, ptr_p1, ptr_p2, ptr_tempo, ip1s_m, ip1s_t, px%data, &
         gz%data, prm%data,tt%data, hu%data, vt%data)
    status = VGD_ERROR
    if(vgd_get(cpg_vgd,"VCOD",vcode) == VGD_ERROR)then
       return
    endif
    if(vcode /= 21001)then
      write(app_msg,*) 'cpg_get_px_21001 called with wrong Vcode, expected 21001 got',vcode
      call app_log(APP_ERROR,app_msg)
       return
    endif
    ier=vgd_get(cpg_vgd,"RFLS",refls_name)
    sleve_L = .not. (trim(refls_name) == trim(VGD_NO_REF_NOMVAR))
    me%nomvar="ME"
    if( cpg_get_rec(me, cpg_lui) == VGD_ERROR ) return
    if(sleve_L )then
       mels%nomvar=refls_name
       if( cpg_get_rec(mels, cpg_lui) == VGD_ERROR ) return
    endif
    thermo_L=.false.
    momentum_L=.false.
    if(trim(cpg_levels_S) == 'THERMO')then
       thermo_L=.true.
    elseif(trim(cpg_levels_S) == 'MOMENTUM')then
       momentum_L=.true.
    elseif(trim(cpg_levels_S) == 'ALL_LEVELS')then
       thermo_L=.true.
       momentum_L=.true.
    else
      write(app_msg,*) 'option -levels expected THERMO, MOMENTUM or ALL_LEVELS but got ',trim(cpg_levels_S)
      call app_log(APP_ERROR,app_msg)
       write(app_msg,*) 'NOTE: passing a nomvar to option -levels is not implemented for vcode ',vcode
       call app_log(APP_VERBATIM,app_msg)
     return
    endif
    ! get all momentum ip1s
    if( vgd_get(cpg_vgd,"vipm - level ip1 list (m)", ip1s_m) == vgd_error)return
    write(app_msg,*) 'List of momentum ip1s',ip1s_m
    call app_log(APP_INFO,app_msg)
    
    ! get all thermo ip1s
    if( vgd_get(cpg_vgd,"vipt - level ip1 list (t)", ip1s_t) == vgd_error)return
    write(app_msg,*) 'List of thermo ip1s_t',ip1s_t
    call app_log(APP_INFO,app_msg)

    ! For vcode 21001, there are nk dynamic temperature, 
    !                  there are nk+2 ip1s_t counting hyb=1.0 and diag levels, these 
    !                  last 2 levels are not used for computing gz.
    nk = size(ip1s_m)-2
    px%nomvar="p0"
    if( cpg_get_rec(px, cpg_lui) == VGD_ERROR ) return
    px%nomvar='PX'; px%ip1=ip1_sfc
    ier = cpg_cp_params(prm,px)    
    if( cpg_zap_etiket_L ) px%etiket = cpg_etiket_S
    if( cpg_fstecr(cpg_luo, px, F_rewrit=.true., data=px%data) == VGD_ERROR)return
    nij = px%ni*px%nj
    allocate(wa(px%ni,px%nj,1),wb(px%ni,px%nj,1),vt%data(px%ni,px%nj),pxt(px%ni,px%nj),stat=ier)
    if( ier /= 0 )then
      write(app_msg,*) 'cpg_get_px_21001: cannot allocate wa, wb, vt%data, pxtsize ->',prm%ni,' x',prm%nj, 'x 1'
      call app_log(APP_ERROR,app_msg)
    endif

    write(app_msg,*) 'Computing 3d pressure for Vcode ',vcode
    call app_log(APP_INFO,app_msg)

    gz%nomvar="gz"
    prm%ip1=ip1_sfc
    !ier=cpg_print_prm(prm)
    if( cpg_get_rec(gz, cpg_lui, match_prm=prm) == VGD_ERROR ) return
    ier = cpg_cp_params(prm,gz)
    wa(:,:,1)=gz%data(:,:)*10.
    ptr_p1 => wa
    ptr_p2 => wb
    aaa=grav/rgasd    
    tt%nomvar="tt"
    vt%nomvar="vt"
    ! Check if VT present
    vt_L = fstinf(cpg_lui,vt%ni,vt%nj,vt%nk,-1,' ',-1,-1,-1,' ',vt%nomvar) >= 0
    hu%nomvar="hu"
    ! start computation from surface up
    do k=nk,1,-1
       !ier = cpg_cp_params(prm,px)
       prm%ip1=ip1s_t(k);
       if(vt_L)then
          if( cpg_get_rec(vt, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          vt%data(:,:) = vt%data(:,:) + tcdk
       else
          if( cpg_get_rec(tt, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          tt%data(:,:) = tt%data(:,:) + tcdk
          if( cpg_get_rec(hu, cpg_lui, match_prm=prm) == VGD_ERROR ) return
          call mfotvt(vt%data,tt%data,hu%data,nij,1,nij)
       endif
       if(sleve_L)then
          if( vgd_levels(cpg_vgd,sfc_field=me%data,sfc_field_ls=mels%data,&
               ip1_list=ip1s_m(k:k),levels=ptr_p2,in_log=.true.) == VGD_ERROR)return
       else
          if( vgd_levels(cpg_vgd,sfc_field=me%data,&
               ip1_list=ip1s_m(k:k),levels=ptr_p2,in_log=.true.) == VGD_ERROR)return
       endif
       do j=1,px%nj
          do i=1,px%ni
             pxt(i,j)=px%data(i,j)
             px%data(i,j)=px%data(i,j)*exp(-aaa/vt%data(i,j)*(ptr_p2(i,j,1)-ptr_p1(i,j,1)))
             ! Geometric average to get PX Thermo which is analitically equivalent to
             !pxt(i,j)=pxt(i,j)*exp(-aaa/vt%data(i,j)*0.5*(ptr_p2(i,j,1)-ptr_p1(i,j,1)))
             if(thermo_L)pxt(i,j)=sqrt(pxt(i,j)*px%data(i,j))
          end do
       end do
       if(vt_L)then
          ier = cpg_cp_params(prm,vt)
       else
          ier = cpg_cp_params(prm,tt)
       endif
       if( cpg_zap_etiket_L ) prm%etiket = cpg_etiket_S
       prm%nomvar='PX'
       if(momentum_L)then
          prm%ip1=ip1s_m(k);
          if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=px%data) == VGD_ERROR)return
       endif
       if(thermo_L)then
          prm%ip1=ip1s_t(k);                 
          if( cpg_fstecr(cpg_luo, prm, F_rewrit=.true., data=pxt) == VGD_ERROR)return
       endif
       if( cpg_zap_etiket_L )prm%etiket=tt%etiket
       ptr_tempo => ptr_p1
       ptr_p1 => ptr_p2
       ptr_p2 => ptr_tempo
    end do    
    deallocate(wa, vt%data, pxt, ip1s_m, ip1s_t, me%data, px%data, gz%data)
    if(associated(tt%data))deallocate(tt%data)
    if(associated(hu%data))deallocate(hu%data)
    if(associated(mels%data))deallocate(mels%data)
    status = VGD_OK
  end function cpg_get_px_21001
  !===================================================================================
  integer function cpg_fstprm(fstkey,record) result(status)
    implicit none
    integer, intent(in) :: fstkey
    type(cpg_rpn) :: record
    ! Local variables
    integer :: error
    real(kind=8) :: nhours
    integer, external :: fstprm
    status = VGD_ERROR
    error=fstprm(fstkey,record%dateo,record%deet,record%npas, &
         record%ni,record%nj,record%nk,record%nbits,record%datyp,record%ip1,record%ip2, &
         record%ip3,record%typvar,record%nomvar,record%etiket,record%grtyp, &
         record%ig1,record%ig2,record%ig3,record%ig4,record%swa, &
         record%lng,record%dltf,record%ubc,record%extra1,record%extra2, &
         record%extra3)
    if (error < 0) then
       write(6,*) 'ERROR: in cpg_fstprm, cannot fstprm for fstkey ',fstkey
       return
    end if
    nhours=record%deet*record%npas/3600.d0
    call incdatr(record%datev,record%dateo,nhours)
    status = VGD_OK
  end function cpg_fstprm
end module mod_comp_pres_gz

!========================================================================

program compute_pressure_gz
   use app
   use vGrid_Descriptors, only: vgd_new,vgd_print,vgd_get,VGD_PRES_TYPE,VGD_HEIGHT_TYPE,VGD_OK,VGD_ERROR
  use mod_comp_pres_gz, only: cpg_process_arguments,cpg_lui,cpg_luo,cpg_samefile_L,cpg_compute_gz_L,&
       cpg_get_vgd_levels,cpg_vgd,cpg_kind,cpg_version,cpg_levels_S,cpg_is_pressure_L,&
       cpg_get_px_gz
  implicit none
  integer :: stat,coor_type, fstfrm
  character(len=12), parameter :: version='1.0.0'

!==========================================================================
app_ptr=app_init(0,'r.compute_pressure_gz',version,'',BUILD_TIMESTAMP)
call app_start()
!==========================================================================

  if(cpg_process_arguments() == VGD_ERROR)then
   stat=app_end(-1)
   error stop
  endif

  ! Determine the vertical coordinate type (pressure or height) in order
  ! to call the appropriate computation function
  if(vgd_new(cpg_vgd,cpg_lui,'fst',kind=cpg_kind,version=cpg_version) == VGD_ERROR)then
     write(app_msg,*) 'vgd_new for levels',cpg_levels_S
     call app_log(APP_ERROR,app_msg)
     stat=app_end(-1)
     error stop
  endif
  if(app_loglevel('').gt.1)then
     if(vgd_print(cpg_vgd) == VGD_ERROR)then
        stat=app_end(-1)
        error stop
     endif
  endif
  if(vgd_get(cpg_vgd,"TYPE",coor_type) == VGD_ERROR)then
     write(app_msg,*) 'vgd_get on key "TYPE"'
     call app_log(APP_ERROR,app_msg)
     stat=app_end(-1)
     error stop
  endif
  if(coor_type == VGD_PRES_TYPE)then
     cpg_is_pressure_L=.true.
     if(cpg_compute_gz_L)then
        if(cpg_get_px_gz() == VGD_ERROR)then
           stat=app_end(-1)
           error stop
        endif
     else
        if(cpg_get_vgd_levels() == VGD_ERROR)then
         stat=app_end(-1)
         error stop
        endif
     endif
  elseif(coor_type == VGD_HEIGHT_TYPE)then
     cpg_is_pressure_L=.false.
     if(cpg_compute_gz_L)then
        if(cpg_get_vgd_levels() == VGD_ERROR) then
         stat=app_end(-1)
         error stop  
        endif      
     else
        if(cpg_get_px_gz() == VGD_ERROR)then
         stat=app_end(-1)
         error stop
        endif
     endif
  else
     write(app_msg,*) 'unsuported vertical coordinate type: ',coor_type
     call app_log(APP_ERROR,app_msg)
     stat=app_end(-1)
     error stop
  endif
  
  stat=fstfrm(cpg_lui)
  if(.not.cpg_samefile_L)stat=fstfrm(cpg_luo)
  !
  stat=app_end(-1)
  !
end program compute_pressure_gz
