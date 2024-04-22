!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

module INs_base
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use IOs
      use INs
      use tdpack
      implicit none
#include <rmnlib_basics.hf>
      
contains

!**   s/r INs_read_mt - Read variable F_var_S and perform horizontal
!                       interpolation to F_nd Arakawa grid destinations

      integer function INs_read_mt ( F_var_S, F_hgrid_S, F_nomvar_S,&
                    F_dest, F_nd, F_ip1, F_nka, F_hint_S, F_quiet_L )
      implicit none
      character(len=*)          ,intent(in)  :: F_var_S
      character(len=*), dimension(*),intent(in) :: F_hgrid_S
      character(len=4), intent(OUT)          :: F_nomvar_S
      character(len=*), optional,intent(in)  :: F_hint_S
      logical         , optional,intent(in)  :: F_quiet_L
      integer                   ,intent(in ) :: F_nd
      integer                   ,intent(out) :: F_nka
      integer, intent(OUT) :: F_ip1(*)
      real, intent(OUT) :: F_dest(G_ni+2*G_halox,G_nj+2*G_haloy,*)

      integer, external :: samegrid_gid, samegrid_rot
      character(len=1) typ,grd
      character(len=4) nomvar,var,dumc
      character(len=12) lab,interp_S
      logical :: quiet_L
      integer, parameter :: nlis = 1024
      integer i, k, idst, err, nz, n1,n2,n3, nrec, liste(nlis),&
              liste_sorted(nlis),lislon,maxdim_wk2
      integer subid,nicore,njcore,datev,ni_dest,nj_dest
      integer mpx,local_nk,irest,kstart, dstf_gid,src_gid, ip1
      integer dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
              dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      real :: surface_level
      real, dimension(:  ), allocatable :: wk1
      real, dimension(:  ), pointer     :: posx,posy
      real(kind=REAL64) add, mult
!
!---------------------------------------------------------------------
!
      INs_read_mt= -1
      F_nka= -1 ; local_nk= 0
      add= 0.d0 ; mult= 1.d0
      quiet_L=.false.
      if (present(F_quiet_L)) quiet_L= F_quiet_L

      nomvar = F_var_S ; ip1= -1
      select case (F_var_S)
         case ('OROGRAPHY')
            if (Inp_kind == 2  ) then
               nomvar= '@NUL'
!!$               if (Inp_src_PX_L) then
!!$                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
!!$                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
!!$               endif
            endif
!!$            if (Inp_kind == 1 ) then
!!$               nomvar= 'GZ' ; ip1= 12000
!!$               if (Inp_src_PX_L) then
!!$                  nomvar= 'GZ' ; surface_level= 1. ; p1=5
!!$                  call convip ( ip1, surface_level,p1,1,dumc,.false. )
!!$               endif
!!$            endif
            if (Inp_kind == 5 ) then
               nomvar= 'GZ' ; surface_level= 1.
               call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( Inp_src_hauteur_L ) then
              ! nomvar= 'GZ'
               nomvar= 'ME'
               ip1=0
              ! if (Inp_kind==21) surface_level= 0.
              ! if (Inp_kind==5 ) surface_level= 1.
              ! call convip ( ip1, surface_level,Inp_kind,1,dumc,.false. )
            endif
            if ( nomvar == 'GZ' ) mult= 10.d0 * grav_8
         case ('SFCPRES')
            nomvar= 'P0'
            if (Inp_kind == 2  ) nomvar= '@NUL'
            !if (Inp_kind == 1  ) nomvar= 'P0'
            !if (Inp_kind == 5  ) nomvar= 'P0'
            !if (Inp_src_hauteur_L ) nomvar= 'P0'
            if ( nomvar == 'P0' ) mult= 100.d0
         case ('TEMPERATURE')
            nomvar= 'TT'
            !if (Inp_kind == 2  ) nomvar= 'TT'
            !if (Inp_kind == 1  ) nomvar= 'TT'
            !if (Inp_kind == 5  ) nomvar= 'TT'
            !if (Inp_src_hauteur_L ) nomvar= 'TT'
            if ( nomvar == 'TT' ) add= tcdk_8
         case ('GEOPOTENTIAL')
            nomvar= 'GZ' ; mult= 10.d0
         case ('PX')
            mult= 100.d0
         case ('URT1')
            mult= knams_8
         case ('VRT1')
            mult= knams_8
      end select

      datev= Inp_cmcdate
      if ( F_var_S(1:min(3,len_trim(F_var_S))) == 'TR/' ) then
         nomvar= F_var_S(4:)
         if (Tr3d_anydate_L) datev= -1
      end if

      F_nomvar_S = trim(nomvar)
      if ( nomvar == '@NUL' ) return

      nrec= fstinl (Inp_handle, n1,n2,n3, datev,' ', &
                    ip1,-1,-1,' ', nomvar,liste,lislon,nlis)
      if (lislon == 0) goto 999

      err= fstprm (liste(1), DTE, DET, IPAS, n1, n2, n3    ,&
                   BIT, DTY, P1, P2, P3, TYP, VAR, LAB, GRD,&
                   G1,G2,G3,G4,SWA,LNG,DLF,UBC,EX1,EX2,EX3)

      src_gid= ezqkdef (n1, n2, GRD, g1, g2, g3, g4, Inp_handle)

      if ((trim(nomvar) == 'URT1').or.(trim(nomvar) == 'VRT1').or.&
          (trim(nomvar) == 'UT1' ).or.(trim(nomvar) == 'VT1' )) then
         err= samegrid_rot ( src_gid, &
         Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4)
         if (err < 0) then
            lislon= 0
            goto 999
         end if
      end if

      call sort_ip1 (liste,liste_sorted,lislon)

      if (lislon > 1) then
         F_ip1(1:lislon) = liste_sorted(1:lislon)
      else
         F_ip1(1) = p1
      end if

      F_nka= lislon
      call splitW ( myproc_IOS, numproc_IOS,lislon,1,&
                    local_nk, kstart, i)
      if (kstart<0) goto 999

      allocate (wk1(n1*n2))
      interp_S= 'CUBIC'
      if (present(F_hint_S)) interp_S= F_hint_S
      
      do idst= 1, F_nd
          
          if (local_nk > 0) then
             if (F_hgrid_S(idst) == 'Q') then
                posx => geomh_longs
                posy => geomh_latgs
             end if
             if (F_hgrid_S(idst) == 'U') then
                posx => geomh_longu
                posy => geomh_latgs
             end if
             if (F_hgrid_S(idst) == 'V') then
                posx => geomh_longs
                posy => geomh_latgv
             end if
             if (F_hgrid_S(idst) == 'F') then
                posx => geomh_longu
                posy => geomh_latgv
             end if
             ni_dest= G_ni+2*G_halox
             nj_dest= G_nj+2*G_haloy
             dstf_gid = ezgdef_fmem (ni_dest, nj_dest, 'Z', 'E', &
                             Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, &
                             posx, posy)

             if ( GRD == 'U' ) then
                nicore = G_ni-Glb_pil_w-Glb_pil_e
                njcore = G_nj-Glb_pil_s-Glb_pil_n
                if (n1 >= nicore .and. n2/2 >= njcore) then
                   subid= samegrid_gid ( &
                   src_gid, Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4,&
                   posx(1+Glb_pil_w), posy(1+Glb_pil_s), nicore,njcore )
                else
                   subid=-1
                end if
                if (subid >= 0) then
                   interp_S = 'NEAREST'
                   err = ezsetopt ('USE_1SUBGRID', 'YES')
                   err = ezsetival('SUBGRIDID', subid)
                else
                   err = ezsetopt ('USE_1SUBGRID', 'NO')
                end if
             end if

             err = ezdefset ( dstf_gid , src_gid )
             err = ezsetopt ('INTERP_DEGREE', interp_S)
             if (lun_out>0) write(lun_out,1001) &
                'Interpolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
                lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),&
                'grid, levels:',kstart,kstart+local_nk-1
          end if

          do i=1,local_nk
             err = fstluk (wk1, liste(kstart+i-1), n1,n2,n3)
             k= (idst-1)*lislon+kstart+i-1
             err = ezsint (F_dest(1,1,k),wk1)
             F_dest(:,:,k)= F_dest(:,:,k)*mult + add
          end do
          if ((err==2).and.(lun_out>0)) &
             write(lun_out,1001) &
             'EXTRApolating: ',trim(F_var_S),trim(nomvar),', nka= ',&
             lislon,',valid: ',Inp_datev,' on ',F_hgrid_S(idst),' grid'
          err= ezsetopt ( 'USE_1SUBGRID', 'NO' )
      end do
      deallocate (wk1)

 999  if (lislon > 0) then
         INs_read_mt= 0
      else
         if ((.not.quiet_L).and.(lun_out>0)) write(lun_out,'(7a)') &
              ' FIELD: ',trim(F_var_S),':',trim(nomvar),' valid: ',&
              Inp_datev, 'NOT FOUND'
      end if
      call MPI_barrier (MY_WORLD_COMM,err)

 1001 format (2a,':',2a,i3,5a,2i3)
!
!---------------------------------------------------------------------
!
      return
      end function INs_read_mt
    
!**s/r INs_read_uv - Read UU and VV and perform horizontal
!                    interpolation to U and V points respectively

      integer function INs_read_uv (F_dest, F_ip1, F_nka)
      implicit none

      integer                  ,intent(out) :: F_nka
      integer, intent(OUT) :: F_ip1(*)
      real, intent(OUT) :: F_dest(G_ni+2*G_halox,G_nj+2*G_haloy,*)

      character(len=1) typ,grd
      character(len=4) var,dumc
      character(len=12) lab
      integer, parameter :: nlis = 1024
      integer :: i, k, idst, err, nz, n1,n2,n3, nrec
      integer :: liste_u(nlis),liste_v(nlis),liste_sorted(nlis)
      integer :: datev,local_nk,src_gid,nku,nkv
      integer :: dstu_gid,dstv_gid,kstart,erru,errv
      integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
                 dty, swa, lng, dlf, ubc, ex1, ex2, ex3
      real, dimension(:), pointer     :: posxu,posyu,posxv,posyv
      real, dimension(:), allocatable :: uv,u,v
!
!---------------------------------------------------------------------
!
      INs_read_uv= -1
      F_nka= -1 ; local_nk= 0

      datev= Inp_cmcdate
      nrec= fstinl (Inp_handle,n1,n2,n3,Inp_cmcdate,' ',-1,-1,-1,' ',&
                    'UU',liste_u,nku,nlis)
      nrec= fstinl (Inp_handle,n1,n2,n3,Inp_cmcdate,' ',-1,-1,-1,' ',&
                    'VV',liste_v,nkv,nlis)
      if ((nku/=nkv).or.(nku<3)) goto 999

      err= fstprm (liste_u(1), DTE, DET, IPAS, n1, n2, n3  ,&
                   BIT, DTY, P1, P2, P3, TYP, VAR, LAB, GRD,&
                   G1,G2,G3,G4,SWA,LNG,DLF,UBC,EX1,EX2,EX3)

      src_gid= ezqkdef (n1, n2, GRD, g1, g2, g3, g4, Inp_handle)

      call sort_ip1 (liste_u,liste_sorted,nku)
      call sort_ip1 (liste_v,liste_sorted,nkv)

      F_ip1(1:nku) = liste_sorted(1:nku)
      F_ip1(nku+1:2*nku) = F_ip1(1:nku)
      
      F_nka= 2*nku
      call splitW ( myproc_IOS, numproc_IOS,nku,1,&
                    local_nk, kstart, i)
      if (kstart<0) goto 999

      allocate (u(n1*n2), v(n1*n2), uv(INs_nid*INs_njd))
      err = ezsetopt ('INTERP_DEGREE', 'CUBIC')
      posxu => geomh_longu
      posyu => geomh_latgs
      posxv => geomh_longs
      posyv => geomh_latgv

      if (local_nk > 0) then
         if (lun_out>0) write(Lun_out,1001) 'Interpolating: UU, nka= ',&
                    nku,', valid: ',Inp_datev,' on U grid'
         dstu_gid = ezgdef_fmem ( INs_nid, INs_njd, 'Z', 'E', &
                           Rot_ig1,Rot_ig2, Rot_ig3, Rot_ig4, &
                                                 posxu, posyu )
         if (lun_out>0) write(Lun_out,1001) 'Interpolating: VV, nka= ',&
                    nkv,', valid: ',Inp_datev,' on V grid'
         dstv_gid = ezgdef_fmem ( INs_nid, INs_njd, 'Z', 'E', &
                          Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, &
                                                 posxv, posyv )
         do i=1,local_nk
            k= kstart+i-1
            err= fstluk ( u, liste_u(k), n1,n2,n3)
            err= fstluk ( v, liste_v(k), n1,n2,n3)
            err = ezdefset ( dstu_gid , src_gid )
            erru= ezuvint  ( F_dest(1,1,k),uv, u,v )
            F_dest(:,:,k)= F_dest(:,:,k) * knams_8
            err = ezdefset ( dstv_gid , src_gid )
            errv= ezuvint  ( uv,F_dest(1,1,nku+k), u,v )
            F_dest(:,:,nku+k)= F_dest(:,:,nku+k) * knams_8
         end do
         if ((erru==2).and.(Lun_out>0)) &
            write(Lun_out,1002) 'EXTRApolating: UU, nka= ',&
                      nku,', valid: ',Inp_datev,' on U grid'
         if ((errv==2).and.(Lun_out>0)) &
            write(Lun_out,1002) 'EXTRApolating: VV, nka= ',&
                      nkv,', valid: ',Inp_datev,' on V grid'
      endif
      deallocate (u,v,uv)

 999  if (nku > 0) then
         INs_read_uv= 0
      else
         if (Lun_out>0) write(Lun_out,'(3a)') &
         'Variable: UU,VV valid: ',Inp_datev, 'NOT FOUND'
      end if
      call MPI_barrier (MY_WORLD_COMM,err)

 1001 format (a,i3,3a)
 1002 format (a,i3,3a)
!
!---------------------------------------------------------------------
!
      return
      end function INs_read_uv

!**s/r INs_hwnd - Read and interpolate horizontal winds UU,VV

      integer function INs_hwnd (F_vname, F_nk, F_deb)
      implicit none

      character(len=*), intent (OUT) :: F_vname
      integer         , intent (OUT) :: F_nk, F_deb
      
      character(len=4) vname
      integer nkau,nkav,err,deb1,deb2,initial
!
!---------------------------------------------------------------------
!
      INs_hwnd= -1 ; initial= INs_nplans
      F_vname= '' ; F_nk=0
      deb1= INs_nplans+1
      deb2= INs_nplans*INs_hord+1
      F_deb= deb1
      err= INs_read_mt ('URT1','U',vname,ND(deb2:),1,DIP1(deb1),nkau)
      if (nkau>0) then
         INs_nplans= INs_nplans+nkau
         deb1= INs_nplans+1
         deb2= INs_nplans*INs_hord+1
         err= INs_read_mt ('VRT1','V',vname,ND(deb2:),1,DIP1(deb1),nkav)
         if (nkav/=nkau) then
            print*, 'URT1 / VRT1 are unusable'
            INs_nplans= initial
            F_deb=-1
         else
            INs_nplans= INs_nplans+nkav
            INs_hwnd= 0
            F_vname= 'UVRT1' ; F_nk=2*nkau
         endif
      endif

      if (INs_hwnd/=0) then
         deb1= INs_nplans+1
         deb2= INs_nplans*INs_hord+1
         F_deb= deb1
         err= INs_read_uv (ND(deb2:),DIP1(deb1),nkau)
         if (nkau>0) then
            INs_nplans= INs_nplans+nkau
            INs_hwnd=0
            F_vname= 'UV' ; F_nk=nkau
         endif
      end if
!     
!---------------------------------------------------------------------
!
      return
      end function INs_hwnd

end module INs_base
