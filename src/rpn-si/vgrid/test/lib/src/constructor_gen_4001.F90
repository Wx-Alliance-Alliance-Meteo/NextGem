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
module mod_constructor_4001
  
  private
  public :: FSTD_ext, my_fstprm

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
 end module mod_constructor_4001

program constructor

  use vGrid_Descriptors, only: vgrid_descriptor,vgd_new,vgd_free,vgd_get,vgd_write,VGD_ERROR
  use Unit_Testing, only: ut_report
  use mod_constructor_4001, only: FSTD_ext, my_fstprm
  

  implicit none

  type(vgrid_descriptor) :: vgd
  integer, parameter :: nk=10
  integer :: stat,k,nl_m=-1,nl_t=-1,nl_w=-1
  integer :: fnom,fstouv,fstfrm,fstluk,fstecr,ier,lu1=11,lu2=12,key,fstinf,ni,nj,nkk
  integer, dimension(:), pointer :: vipm,vipt,vipw,work_i
  real :: my_epsilon=1.e-6, dummy
  real, dimension(nk) :: hgts,hgts2
  real, dimension(:), pointer :: vcdm,vcdt,vcdw,work
  real, dimension(:,:), allocatable :: ff
  real(kind=8), dimension(:), pointer :: a_m_8,b_m_8,c_m_8,a_t_8,b_t_8,c_t_8,a_w_8,b_w_8,c_w_8,work_8
  logical :: OK=.true.
  logical, parameter :: write_control_L=.false.
  character (len=256) :: file
  type(FSTD_ext) :: fst

  nullify(vipm,vipt,vipw,work_i,vcdm,vcdt,vcdw,work,a_m_8,b_m_8,c_m_8,a_t_8,b_t_8,c_t_8,a_w_8,b_w_8,c_w_8,work_8)

  hgts=(/0.,20.,25.,40.,50.,80.,100.,120.,150.,200./)

  if( vgd_new(vgd,kind=4,version=1,hyb=hgts) == VGD_ERROR)error stop 1

  if( vgd_get(vgd,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='CA_T - vertical A coefficient (t)'   ,value=a_t_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='CA_W - vertical A coefficient (w)'   ,value=a_w_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='CB_M - vertical B coefficient (m)'   ,value=b_m_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='CB_T - vertical B coefficient (t)'   ,value=b_t_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='CB_W - vertical B coefficient (w)'   ,value=b_w_8) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VIPM - level ip1 list (m)'           ,value=vipm) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VIPT - level ip1 list (t)'           ,value=vipt) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VIPW - level ip1 list (w)'           ,value=vipw) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VCDM - vertical coordinate (m)'      ,value=vcdm) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VCDT - vertical coordinate (t)'      ,value=vcdt) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='VCDW - vertical coordinate (w)'      ,value=vcdw) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='NL_M - Number of vertical levels (m)',value=nl_m) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='NL_T - Number of vertical levels (t)',value=nl_t) == VGD_ERROR )error stop 1
  if( vgd_get(vgd,key='NL_W - Number of vertical levels (w)',value=nl_w) == VGD_ERROR )error stop 1

  file='data/data_constructor_gen_4001.txt'
  if(write_control_L)then
     open(unit=10,file=file)
     write(10,*)a_m_8
     write(10,*)vipm
     write(10,*)vcdm
     close(10)
  endif

  if(nl_m.ne.size(a_m_8))then
     OK=.false.
     print*,'Problem with NL_M 4001 should be',size(b_m_8),' got',nl_m
  endif
  if(nl_t.ne.size(b_m_8))then
     OK=.false.
     print*,'Problem with NL_T 4001 should be',size(b_m_8),' got',nl_t
  endif
  if(nl_w.ne.size(b_m_8))then
     OK=.false.
     print*,'Problem with NL_W 4001 should be',size(b_m_8),' got',nl_w
  endif

  ! Check 1
  open(unit=10,file=file,ACTION='READ')
  allocate(work_8(size(a_m_8)))
  print*,'Reading A'
  read(10,*)work_8
  do k=1,nk
     if( abs(work_8(k)) < my_epsilon .or. abs(a_m_8(k)) < my_epsilon )then
        ! work_8(k) is probably zero
        if( abs(work_8(k)-a_m_8(k)) > my_epsilon)then
           OK=.false.
           print*,'Probleme avec A, pas dans les limites tollerees'
           print*,work_8(k),'vs',a_m_8(k)
        endif
     else
        if(abs(work_8(k)-a_m_8(k))/a_m_8(k)>my_epsilon)then
           OK=.false.
           print*,'Probleme avec A, pas dans les limites tollerees'
           print*,work_8(k),'vs',a_m_8(k)
        endif
     endif
  enddo

  ! Check B
  if(any(b_m_8.ne.0.d0))then
     OK=.false.
     print*,'Probleme avec B M, doit egaler 0.0d :',b_m_8
  endif
  if(any(b_t_8.ne.0.d0))then
     OK=.false.
     print*,'Probleme avec B T, doit egaler 0.0d :',b_t_8
  endif
  if(any(b_w_8.ne.0.d0))then
     OK=.false.
     print*,'Probleme avec B W, doit egaler 0.0d :',b_w_8
  endif

  if(any(a_m_8.ne.a_t_8))then
     OK=.false.
     print*,'Probleme avec A T, dot egaler A M :'
     print*,a_t_8
     print*,a_m_8
  endif

  if(any(a_w_8.ne.a_t_8))then
     OK=.false.
     print*,'Probleme avec A W, dot egaler A M :'
     print*,a_t_8
     print*,a_w_8
  endif

  ! Check IP
  allocate(work_i(size(vipm)))
  print*,'Reading IP'
  read(10,*)work_i
  do k=1,nk
     if(work_i(k).ne.vipm(k))then
        OK=.false.
        print*,'Probleme avec IP :'
        print*,work_i(k),'vs',vipm(k)
     endif
  enddo
  if(any(vipt.ne.vipm))then
     OK=.false.
     print*,'Probleme avec ip1 T , doit egaler ip1 M :'
     print*,vipt
     print*,vipm
  endif
  if(any(vipw.ne.vipm))then
     OK=.false.
     print*,'Probleme avec ip1 W , doit egaler ip1 M :'
     print*,vipt
     print*,vipw
  endif

  ! Check vcdm
  allocate(work(size(vcdm)))
  print*,'Reading VCDM'
  read(10,*)work
  do k=1,nk
     if(abs(work(k)) < my_epsilon .or. abs(vcdm(k)) < my_epsilon)then
        if( abs(work(k)-vcdm(k)) > my_epsilon)then
           OK=.false.
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs',vcdm(k)
        endif
     else
        if(abs(work(k)-vcdm(k))/vcdm(k)>my_epsilon)then
           OK=.false.
           print*,'Probleme avec vcdm, pas dans les limites tollerees'
           print*,work(k),'vs',vcdm(k)
        endif
     endif
  enddo  
  if(any(vcdt.ne.vcdm))then
     OK=.false.
     print*,'Probleme avec VCDT, doit egaler VCDM'
     print*,vcdt
     print*,vcdm
  endif
  if(any(vcdw.ne.vcdm))then
     OK=.false.
     print*,'Probleme avec VCDW, doit egaler VCDM'
     print*,vcdw
     print*,vcdm
  endif
  
  print*,'The following error is normal'
  print*,'Testing height not positive'
  stat = vgd_new(vgd,kind=4,version=1,hyb=-hgts)
  if(stat.ne.VGD_ERROR)OK=.false.

  print*,'The following error is normal'
  print*,'Testing sig levels not in order'
  hgts2=hgts
  hgts2(size(hgts)-3)=hgts2(size(hgts)-2)
  stat = vgd_new(vgd,kind=1,version=1,hyb=hgts2)
  if(stat.ne.VGD_ERROR)OK=.false.
  deallocate(a_m_8,a_t_8,b_m_8,b_t_8,vipm,vipt,vcdm)
  stat = vgd_free(vgd)

  if(.true.) then
     ! Generate file to copy in data/dm_4001_from_model_run
     if( vgd_new(vgd,kind=4,version=1,hyb=hgts) == VGD_ERROR)error stop 1
     if( vgd_get(vgd,key='NL_M - Number of vertical levels (m)',value=nl_m) == VGD_ERROR )error stop 1
     if( vgd_get(vgd,key='CA_M - vertical A coefficient (m)'   ,value=a_m_8) == VGD_ERROR )error stop 1
      if( vgd_get(vgd,key='VIPM - level ip1 list (m)'           ,value=vipm) == VGD_ERROR )error stop 1
     file="data_out/to_copy_in_dm_4001_from_model_run"
     if( fnom(lu1,file,"RND",0) < 0 )then
        print*,'ERROR with fnom on file ',trim(file)
        error stop 1
     endif
     ier = fstouv(lu1,'RND')
     ! Use the following file only the get a record info
     file="data/dm_21002_from_model_run_NON_SLEVE"
     if( fnom(lu2,file,"RND",0) < 0 )then
        print*,'ERROR with fnom on file ',trim(file)
        error stop 1
     endif
     if( fstouv(lu2,'RND') <= 0 )then
        print*,'file empty',trim(file)
     endif
     key = fstinf(lu2,ni,nj,nkk,-1,' ',-1,-1,-1,' ','GZ')
     if( key < 0 )then
        print*,'ERROR cannot find GZ at 93423264, key=',key
        error stop 1
     endif
     ier = my_fstprm(key,fst)
     allocate(ff(fst%ni,fst%ni))
     if( fstluk(ff,key,ni,nj,nkk) < 0 )error stop 1     
     do k=1,nl_m
        ff=a_m_8(k)/10.
        ier=fstecr(ff,dummy,-32,lu1,fst%dateo,fst%deet,fst%npas, &
             ni,nj,nkk,vipm(k),fst%ip2,fst%ip3,fst%typvar,fst%nomvar,fst%etiket,fst%grtyp, &
             fst%ig1,fst%ig2,fst%ig3,fst%ig4,fst%datyp,.true.)
     enddo
     if( vgd_write(vgd,lu1,format='fst') == VGD_ERROR )error stop 1
     ier = fstfrm(lu1)
     ier = fstfrm(lu2)
  endif
  
  stat = vgd_free(vgd)

  call ut_report(OK,'Grid_Descriptors::vgd_new vertical generate initializer (4001) value')

end program constructor
