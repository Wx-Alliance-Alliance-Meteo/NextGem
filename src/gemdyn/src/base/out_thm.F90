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

!**s/r out_thm - output  temperature, humidity and mass fields

      subroutine out_thm ( levset, set )
      use dynkernel_options
      use metric
      use vertical_interpolation
      use vGrid_Descriptors
      use vgrid_wb
      use geomh
      use gmm_vt1
      use gmm_pw
      use gmm_geof
      use VERgrid_options
      use gem_options
      use out_options
      use step_options
      use dyn_fisl_options
      use tdpack
      use glb_ld
      use cstv
      use out_meta
      use out_mod
      use out3
      use levels
      use outp
      use outd
      use ver
      use mem_tstp
      use rmn_gmm
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: levset, set

      type :: stg_i
         integer :: t,m,p
      end type stg_i

      type(vgrid_descriptor) :: vcoord

      logical,save :: done_L= .false.
      logical      :: write_diag_lev,near_sfc_L
      logical      :: satues_L= .false.

      integer,save :: lastdt= -1
      integer i,j,k,ii,l_ninj,nko,istat,nk_under,nk_src,knd,nkindx
      integer pngz,pnvt,pntt,pnes,pntd, pnhr,pnpx,pntw,pnwe,pnww,&
              pnzz,pnth,pnpn,pnp0,psum,pnla,pnlo,pnme,pnmx
      integer psl1,psl2,psl3,&
              psl4,psl5,psl6,&
              psl7,psl8,psl9,&
              psl10,psl11,psl12
      integer, dimension(:), allocatable :: indo
      integer, dimension(:), pointer     :: ip1m

      real  , parameter :: theta_p0 = 100000.
      real  , parameter :: ES_MAX   = 30.
      real(kind=REAL64), parameter :: ZERO_8   = 0.0

      real w1(l_minx:l_maxx,l_miny:l_maxy), w2(l_minx:l_maxx,l_miny:l_maxy),&
           p0(l_minx:l_maxx,l_miny:l_maxy), deg2rad, zd2etad

      real, dimension(:,:,:), pointer    :: gmm_hut1, wlnph_m, wlnph_ta,www
      real ,dimension(:,:,:), allocatable:: px_pres,hu_pres,td_pres    ,&
                                            tt_pres,vt_pres,w5,w6,cible,&
                                            gzm,gzt,ttx,htx,ffwe       ,&
                                            px_ta,px_m,th,t8,myomega   
      real, dimension(:,:  ), allocatable:: wlao
      real ,dimension(:    ), allocatable:: prprlvl,rf
      real, dimension(:    ), pointer    :: hybm,hybt,hybt_w
      integer ind0(1) ! One level output
      real hyb0(1),hybt_gnk1(1),hybt_gnk2(1) ! One level output

      real tt (l_minx:l_maxx,l_miny:l_maxy,G_nk+1),&
           hu (l_minx:l_maxx,l_miny:l_maxy,G_nk+1),&
           vt (l_minx:l_maxx,l_miny:l_maxy,G_nk+1)

      save gzm,gzt,htx,ttx,wlao,hybm,hybt,hybt_w
!
!-------------------------------------------------------------------
!
      l_ninj= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)

      pnpn=0 ; pnp0=0 ; pnla=0 ; pnlo=0 ; pnme=0 ; pnmx=0
      pngz=0 ; pnvt=0 ; pntt=0 ; pnes=0 ; pntd=0 ; pnhr=0 ; pnpx=0
      pntw=0 ; pnwe=0 ; pnww=0 ; pnzz=0 ; pnth=0 ; psl1=0 ; psl2=0
      psl3=0 ; psl4=0 ; psl5=0 ; psl6=0 ; psl7=0 ; psl8=0 ; psl9=0
      psl10=0; psl11=0; psl12=0

      hyb0(1)=0.0
      ind0(1)=1

      do ii=1,Outd_var_max(set)
         if (Outd_var_S(ii,set) == 'PN') pnpn=ii
         if (Outd_var_S(ii,set) == 'P0') pnp0=ii
         if (Outd_var_S(ii,set) == 'LA') pnla=ii
         if (Outd_var_S(ii,set) == 'LO') pnlo=ii
         if (Outd_var_S(ii,set) == 'ME') pnme=ii
         if (Outd_var_S(ii,set) == 'MX') pnmx=ii
         if (Outd_var_S(ii,set) == 'GZ') pngz=ii
         if (Outd_var_S(ii,set) == 'VT') pnvt=ii
         if (Outd_var_S(ii,set) == 'TT') pntt=ii
         if (Outd_var_S(ii,set) == 'ES') pnes=ii
         if (Outd_var_S(ii,set) == 'TD') pntd=ii
         if (Outd_var_S(ii,set) == 'HR') pnhr=ii
         if (Outd_var_S(ii,set) == 'PX') pnpx=ii
         if (Outd_var_S(ii,set) == 'TW') pntw=ii
         if (Outd_var_S(ii,set) == 'WE') pnwe=ii
         if (Outd_var_S(ii,set) == 'WW') pnww=ii
         if (Outd_var_S(ii,set) == 'ZZ') pnzz=ii
         if (Outd_var_S(ii,set) == 'TH') pnth=ii
         if (Outd_var_S(ii,set) == 'S1') psl1=ii
         if (Outd_var_S(ii,set) == 'S2') psl2=ii
         if (Outd_var_S(ii,set) == 'S3') psl3=ii
         if (Outd_var_S(ii,set) == 'S4') psl4=ii
         if (Outd_var_S(ii,set) == 'S5') psl5=ii
         if (Outd_var_S(ii,set) == 'S6') psl6=ii
         if (Outd_var_S(ii,set) == 'S7') psl7=ii
         if (Outd_var_S(ii,set) == 'S8') psl8=ii
         if (Outd_var_S(ii,set) == 'S9') psl9=ii
         if (Outd_var_S(ii,set) == 'S10') psl10=ii
         if (Outd_var_S(ii,set) == 'S11') psl11=ii
         if (Outd_var_S(ii,set) == 'S12') psl12=ii
      end do

      psum=pnpn+pnp0+pnla+pnlo+pnme+pnmx
      psum=psum +  &
           pngz+pnvt+pntt+pnes+pntd+pnhr+pnpx+ &
           pntw+pnwe+pnww+pnzz+pnth+psl1+psl2+psl3+&
           psl4+psl5+psl6+psl7+psl8+psl9+psl10+psl11+psl12

      if (psum == 0) return

      if (pnww /= 0) allocate ( myomega(l_minx:l_maxx,l_miny:l_maxy,G_nk  ) )
      if (pnth /= 0) allocate ( th   (l_minx:l_maxx,l_miny:l_maxy,G_nk+1) )

!     Obtain humidity HUT1 and other GMM variables
      nullify (gmm_hut1,wlnph_m,wlnph_ta)
      istat= gmm_get('TR/'//'HU'//':P', gmm_hut1  )
      istat= gmm_get(gmmk_pw_log_pm_s , wlnph_m   )
      istat= gmm_get(gmmk_pw_log_pt_s , wlnph_ta  )
      call out_padbuf (wlnph_m ,l_minx,l_maxx,l_miny,l_maxy,G_nk+1)
      call out_padbuf (wlnph_ta,l_minx,l_maxx,l_miny,l_maxy,G_nk+1)

!     Determine number of output levels
      nk_src= l_nk
      if (Out3_sfcdiag_L) nk_src= l_nk+1

!     Obtain HU from HUT1 and physics diag level
      hu(:,:,1:l_nk) = gmm_hut1(:,:,1:l_nk)
      hu(:,:,l_nk+1) = qdiag
      call out_padbuf(hu,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)
!
!     Store temperature TT (in tt)
!
      tt(1:l_ni,1:l_nj,1:G_nk) = pw_tt_plus(1:l_ni,1:l_nj,1:G_nk)
      tt(:,:,G_nk+1)= tdiag
      call out_padbuf(tt,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)

!     Obtain Virtual temperature from TT1 and physics diag level
!     On diag level there is no hydrometeor.

      vt(:,:,1:l_nk)= tt1(:,:,1:l_nk)
      vt(:,:,l_nk+1)= tt1(:,:,l_nk  )
      if (Out3_sfcdiag_L) then
         call mfotvt (vt(l_minx,l_miny,nk_src),tt(l_minx,l_miny,nk_src),&
                      hu(l_minx,l_miny,nk_src),l_ninj,1,l_ninj)
      end if
      call out_padbuf(vt,l_minx,l_maxx,l_miny,l_maxy,l_nk+1)

      p0= pw_p0_plus
      call out_padbuf(p0,l_minx,l_maxx,l_miny,l_maxy,1)

!_________________________________________________________________
!
!     2.0    Output 2D variables
!_________________________________________________________________
!     output 2D fields on 0mb (pressure)
      knd=2
      Out_stag_S= 'MS '

      if (pnme /= 0)then
         call out_fstecr(fis0,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
              'ME  ',Outd_convmult(pnme,set),Outd_convadd(pnme,set),&
              knd,-1,1,ind0, 1, Outd_nbit(pnme,set),.false. )
            if(Schm_sleve_L)then
               call out_fstecr(orols,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
                  'MELS',Outd_convmult(pnme,set),Outd_convadd(pnme,set),&
                  knd,-1,1,ind0, 1, Outd_nbit(pnme,set),.false. )
            endif
         end if
      if (pnmx /= 0)then
            call out_fstecr(fis0,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
              'MX  ',Outd_convmult(pnmx,set),Outd_convadd(pnmx,set),&
              knd,-1,1,ind0, 1, Outd_nbit(pnmx,set),.false. )
         end if
      if (pnla /= 0) &
          call out_fstecr(geomh_latrx,1,l_ni,1,l_nj,hyb0, &
              'LA  ',Outd_convmult(pnla,set),Outd_convadd(pnla,set),&
              knd,-1,1,ind0, 1, Outd_nbit(pnla,set),.false. )
      if (pnlo /= 0) &
          call out_fstecr(geomh_lonrx,1,l_ni,1,l_nj,hyb0, &
              'LO  ',Outd_convmult(pnlo,set),Outd_convadd(pnlo,set),&
              knd,-1,1, ind0, 1, Outd_nbit(pnlo,set),.false. )
!_______________________________________________________________________
!
!     3.0    Precomputations for output over pressure levels or PN or
!            GZ on thermo levels
!
!        The underground extrapolation can use precalculated
!        temperatures over fictitious underground geopotential levels.
!        The levels in meters are stored in Out3_lieb_levels(Out3_lieb_nk).
!        Out3_lieb_levels is a user's given parameters.
!_______________________________________________________________________
!
      nk_under = Out3_lieb_nk

      If (.not.done_L) then
          lastdt= Lctl_step-1
          done_L=.true.
          allocate ( ttx (l_minx:l_maxx,l_miny:l_maxy,nk_under),&
                     htx (l_minx:l_maxx,l_miny:l_maxy,nk_under),&
                     gzt (l_minx:l_maxx,l_miny:l_maxy,G_nk+1  ),&
                     gzm (l_minx:l_maxx,l_miny:l_maxy,G_nk+1  ),&
                     wlao(l_minx:l_maxx,l_miny:l_maxy))
!         Store WLAO (latitude in rad)
          deg2rad= pi_8/180.d0
          do j=1,l_nj
          do i=1,l_ni
             wlao (i,j) = geomh_latrx(i,j) * deg2rad
          end do
          end do
      end if

!     Compute GZ on thermo levels in gzt
!     Compute ttx and htx (underground)

      if ( lastdt /= Lctl_step ) then

         if (.not.Dynamics_sw_L) then
            gzm(1:l_ni,1:l_nj,1:G_nk+1) = grav_8 * GVM%zmom_8(1:l_ni,1:l_nj,1:G_nk+1)
            gzt(1:l_ni,1:l_nj,1:G_nk+1) = grav_8 * GVM%ztht_8(1:l_ni,1:l_nj,1:G_nk+1)
         else
            gzm(1:l_ni,1:l_nj,1:G_nk+1) = grav_8 *(qt1(1:l_ni,1:l_nj,1:G_nk+1) + Ver_z_8%m(1))
            gzt(1:l_ni,1:l_nj,1:G_nk+1) = grav_8 *(qt1(1:l_ni,1:l_nj,1:G_nk+1) + Ver_z_8%m(1))
         endif

         call out_liebman (ttx, htx, vt, gzt, fis0, wlao, &
                        l_minx,l_maxx,l_miny,l_maxy,Out3_lieb_nk,nk_src)

      end if
                           
      lastdt = Lctl_step

!     Compute PN
      if (pnpn /= 0) then
         w2= log(p0)
         call pnm (w1, vt(l_minx,l_miny,nk_src), fis0, w2, wlao, &
                   ttx, htx, nk_under, l_minx, l_maxx, l_miny, l_maxy, 1)
         if (Outd_filtpass(pnpn,set) > 0) &
             call filter( w1,Outd_filtpass(pnpn,set),Outd_filtcoef(pnpn,set),&
                           l_minx,l_maxx,l_miny,l_maxy,1)
         call out_fstecr( w1,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
              'PN  ',Outd_convmult(pnpn,set),Outd_convadd(pnpn,set), &
              knd,-1,1, ind0, 1, Outd_nbit(pnpn,set),.false. )
      end if

!     Compute P0
      if (pnp0 /= 0) then
         do j=l_miny,l_maxy
         do i=l_minx,l_maxx
            w1(i,j) = p0(i,j)
         end do
         end do
         if (Outd_filtpass(pnp0,set) > 0)&
         call filter( w1,Outd_filtpass(pnp0,set),Outd_filtcoef(pnp0,set), &
                       l_minx,l_maxx,l_miny,l_maxy,1)
         call out_fstecr (w1,l_minx,l_maxx,l_miny,l_maxy,hyb0,&
              'P0  ',Outd_convmult(pnp0,set),Outd_convadd(pnp0,set), &
              knd,-1,1,ind0, 1, Outd_nbit(pnp0,set),.false.)
      end if

      if (pnww /= 0) then
         call calomeg_w(myomega,orols,wt1(l_minx,l_miny,1),tt1(l_minx,l_miny,1),wlnph_ta,l_minx,l_maxx,l_miny,l_maxy,G_nk)
      end if

      if (pnth /= 0) then
         do k= 1, G_nk+1
            do j= 1,l_nj
            do i= 1,l_ni
               th(i,j,k)= tt(i,j,k)*(theta_p0/exp(wlnph_ta(i,j,k)))**cappa_8
            end do
            end do
         end do
      end if

      if (Level_typ_S(levset) == 'M') then  ! Output on model levels

!       Setup the indexing for output
         Out_stag_S= 'MT '
         knd= Level_kind_ip1
         allocate ( indo(G_nk+1) )
         call out_slev ( Level(1,levset), Level_max(levset), &
                          Level_momentum,indo,nko,near_sfc_L)
         write_diag_lev= near_sfc_L .and. out3_sfcdiag_L
         indo(G_nk+1) = G_nk+1
         
!        Retrieve vertical coordinate description
         if ( .not. associated (hybm) ) then
            nullify(ip1m,hybm,hybt,hybt_w)
            istat = vgrid_wb_get('ref-m',vcoord,ip1m)
            deallocate(ip1m); nullify(ip1m)
            if (vgd_get(vcoord,'VCDM - vertical coordinate (m)',hybm) /= VGD_OK) istat = VGD_ERROR
            if (vgd_get(vcoord,'VCDT - vertical coordinate (t)',hybt) /= VGD_OK) istat = VGD_ERROR
            istat = vgd_free(vcoord)
            allocate(hybt_w(G_nk))
            ! For vertical motion quantities, we place level NK at the surface
            hybt_w(1:G_nk)= hybt(1:G_nk)
         end if
         hybt_gnk1(1)=hybt(G_nk+1)
         hybt_gnk2(1)=hybt(G_nk+2)

         if (psl1 /= 0)then
            call out_fstecr(sw_f1,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S1  ',Outd_convmult(psl1,set),Outd_convadd(psl1,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl1,set),.false. )
         end if

         if (psl2 /= 0)then
            call out_fstecr(sw_f2,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S2  ',Outd_convmult(psl2,set),Outd_convadd(psl2,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl2,set),.false. )
         end if

         if (psl3 /= 0)then
            call out_fstecr(sw_f3,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S3  ',Outd_convmult(psl3,set),Outd_convadd(psl3,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl3,set),.false. )
         end if

         if (psl4 /= 0)then
            call out_fstecr(sw_f4,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S4  ',Outd_convmult(psl4,set),Outd_convadd(psl4,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl4,set),.false. )
         end if

         if (psl5 /= 0)then
            call out_fstecr(sw_f5,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S5  ',Outd_convmult(psl5,set),Outd_convadd(psl5,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl5,set),.false. )
         end if

         if (psl6 /= 0)then
            call out_fstecr(sw_f6,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S6  ',Outd_convmult(psl6,set),Outd_convadd(psl6,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl6,set),.false. )
         end if

         if (psl7 /= 0)then
            call out_fstecr(sw_f7,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S7  ',Outd_convmult(psl7,set),Outd_convadd(psl7,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl7,set),.false. )
         end if

         if (psl8 /= 0)then
            call out_fstecr(sw_f8,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S8  ',Outd_convmult(psl8,set),Outd_convadd(psl8,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl8,set),.false. )
         end if

         if (psl9 /= 0)then
            call out_fstecr(sw_f9,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S9  ',Outd_convmult(psl9,set),Outd_convadd(psl9,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl9,set),.false. )
         end if

         if (psl10 /= 0)then
            call out_fstecr(sw_f10,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S10 ',Outd_convmult(psl10,set),Outd_convadd(psl10,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl10,set),.false. )
         end if

         if (psl11 /= 0)then
            call out_fstecr(sw_f11,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S11 ',Outd_convmult(psl11,set),Outd_convadd(psl11,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl11,set),.false. )
         end if

         if (psl12 /= 0)then
            call out_fstecr(sw_f12,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'S12 ',Outd_convmult(psl12,set),Outd_convadd(psl12,set),&
               knd,-1,G_nk,indo,nko,Outd_nbit(psl12,set),.false. )
         end if

         if (pngz /= 0)then
            Out_stag_S= 'MM '
            call out_fstecr(gzm,l_minx,l_maxx,l_miny,l_maxy,hybm, &
               'GZ  ',Outd_convmult(pngz,set),Outd_convadd(pngz,set),&
               knd,-1,G_nk+1,indo,nko,Outd_nbit(pngz,set),.false. )
            Out_stag_S= 'MT ' ; nkindx= nko
            if (near_sfc_L) nkindx= nko+1
            call out_fstecr(gzt,l_minx,l_maxx,l_miny,l_maxy,hybt, &
               'GZ  ',Outd_convmult(pngz,set),Outd_convadd(pngz,set),&
               knd,-1,G_nk+1,indo,nkindx,Outd_nbit(pngz,set),.false. )
         end if

         if (pnvt /= 0)then
            call out_fstecr(vt,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                 'VT  ',Outd_convmult(pnvt,set),Outd_convadd(pnvt,set),&
                 knd,-1,G_nk+1,indo,nko,Outd_nbit(pnvt,set),.false. )
            if (write_diag_lev) then
               Out_stag_S(3:3)= 'D'
               call out_fstecr(vt(l_minx,l_miny,G_nk+1),&
                               l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                    'VT  ',Outd_convmult(pnvt,set),Outd_convadd(pnvt,set),&
                    Level_kind_diag,-1,1,ind0,1,Outd_nbit(pnvt,set),.false. )
               Out_stag_S(3:3)= ' '
            end if
         end if
         if (pnth /= 0) then
            call out_fstecr(th,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                    'TH  ',Outd_convmult(pnth,set),Outd_convadd(pnth,set),&
                    knd,-1,G_nk+1,indo,nko,Outd_nbit(pnth,set),.false. )
            if (write_diag_lev) then
               Out_stag_S(3:3)= 'D'
               call out_fstecr(th(l_minx,l_miny,G_nk+1),&
                                  l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                       'TH  ',Outd_convmult(pnth,set),Outd_convadd(pnth,set),&
                       Level_kind_diag,-1,1,ind0,1,Outd_nbit(pnth,set),.false. )
               Out_stag_S(3:3)= ' '
            end if
         end if

         if (pntt /= 0)then
            call out_fstecr(tt,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                 'TT  ' ,Outd_convmult(pntt,set),Outd_convadd(pntt,set), &
                 knd,-1, G_nk+1,indo,nko,Outd_nbit(pntt,set),.false. )
            if (write_diag_lev) then
               Out_stag_S(3:3)= 'D'
               call out_fstecr(tt(l_minx,l_miny,G_nk+1),&
                               l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                 'TT  ',Outd_convmult(pntt,set),Outd_convadd(pntt,set),&
                 Level_kind_diag,-1,1,ind0,1,Outd_nbit(pntt,set),.false. )
               Out_stag_S(3:3)= ' '
            end if
         end if

         if (pnes /= 0.or.pnpx /= 0.or.pntw /= 0.or.pntd /= 0.or.pnhr /= 0) then

            allocate ( px_ta(l_minx:l_maxx,l_miny:l_maxy,G_nk+1),&
                       px_m (l_minx:l_maxx,l_miny:l_maxy,G_nk+1) )

!            Compute PX on thermo levels
             px_ta(:,:,1:G_nk)= exp(wlnph_ta(:,:,1:G_nk))
             px_ta(:,:,G_nk+1)= p0

!            Compute PX on momentum levels
             px_m(:,:,1:G_nk)= exp(wlnph_m(:,:,1:G_nk))
             px_m(:,:,G_nk+1)= p0

         end if

         if (pnpx /= 0)then
            Out_stag_S= 'MM '
            call out_fstecr(px_m,l_minx,l_maxx,l_miny,l_maxy,hybm, &
                 'PX  ',Outd_convmult(pnpx,set),Outd_convadd(pnpx,set), &
                  knd,-1,G_nk+1,indo,nko,Outd_nbit(pnpx,set),.false. )
            Out_stag_S= 'MT ' ; nkindx= nko
            if (near_sfc_L) nkindx= nko+1
            call out_fstecr(px_ta,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                  'PX  ',Outd_convmult(pnpx,set),Outd_convadd(pnpx,set),&
                  knd,-1,G_nk+1,indo,nkindx,Outd_nbit(pnpx,set),.false. )
         end if

         if (pnes /= 0.or.pntw /= 0.or.pntd /= 0.or.pnhr /= 0) &
               allocate (t8 (l_minx:l_maxx,l_miny:l_maxy,G_nk+1) )

         if (pntw /= 0) then
!        Calculate THETAW TW (t8=TW) (px=PX)
             call mthtaw4 (t8,hu,tt, px_ta,satues_l, &
                           .true.,trpl_8,l_ninj,nk_src,l_ninj)
             call out_fstecr(t8,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                  'TW  ',Outd_convmult(pntw,set),Outd_convadd(pntw,set), &
                  knd,-1,G_nk+1, indo, nko, Outd_nbit(pntw,set),.false. )
             if (write_diag_lev) then
                Out_stag_S(3:3)= 'D'
                call out_fstecr(t8(l_minx,l_miny,G_nk+1),&
                     l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                     'TW  ',Outd_convmult(pntw,set),Outd_convadd(pntw,set),&
                     Level_kind_diag,-1,1,ind0,1, Outd_nbit(pntw,set),.false. )
                Out_stag_S(3:3)= ' '
             end if
         end if

         if (pnes /= 0 .or. pntd /= 0) then
!        Calculate ES (t8=ES) (px=PX)
            call mhuaes3 (t8,hu,tt,px_ta,satues_l, &
                                  l_ninj,nk_src,l_ninj)

            t8(1:l_ni,1:l_nj,1:nk_src) = min(t8(1:l_ni,1:l_nj,1:nk_src), ES_MAX)
            if (Out3_cliph_L) then
               t8(1:l_ni,1:l_nj,1:nk_src) = max(t8(1:l_ni,1:l_nj,1:nk_src), 0.)
            end if

            if (pnes /= 0) then
               call out_fstecr(t8,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                    'ES  ',Outd_convmult(pnes,set),Outd_convadd(pnes,set),&
                    knd,-1,G_nk+1,indo,nko,Outd_nbit(pnes,set),.false. )
               if (write_diag_lev) then
                  Out_stag_S(3:3)= 'D'
                  call out_fstecr(t8(l_minx,l_miny,G_nk+1), &
                                  l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                       'ES  ',Outd_convmult(pnes,set),Outd_convadd(pnes,set),&
                       Level_kind_diag,-1,1,ind0,1,Outd_nbit(pnes,set),.false. )
                       Out_stag_S(3:3)= ' '
               end if
            end if

            if (pntd /= 0) then
!            Calculate TD (tt=TT,t8=old ES, t8=TD=TT-ES)
               do k= 1,nk_src
                  do j= 1,l_nj
                  do i= 1,l_ni
                     t8(i,j,k) = tt(i,j,k) - t8(i,j,k)
                  end do
                  end do
               end do
               call out_fstecr(t8,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                    'TD  ',Outd_convmult(pntd,set),Outd_convadd(pntd,set),&
                    knd,-1,G_nk+1,indo,nko,Outd_nbit(pntd,set),.false. )
               if (write_diag_lev) then
                  Out_stag_S(3:3)= 'D'
                  call out_fstecr(t8(l_minx,l_miny,G_nk+1), &
                       l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                       'TD  ',Outd_convmult(pntd,set),Outd_convadd(pntd,set),&
                       Level_kind_diag,-1,1,ind0,1,Outd_nbit(pntd,set),.false. )
                   Out_stag_S(3:3)= ' '
               end if
            end if
         end if

         if (pnhr /= 0) then
!            Calculate HR (t8=HR,tt=TT,px=PX)
            call mfohr4 (t8,hu,tt,px_ta,l_ninj,nk_src,l_ninj,satues_l)
            if ( Out3_cliph_L ) then
               do k= 1,nk_src
                  do j= 1,l_nj
                  do i= 1,l_ni
                     t8(i,j,k)= max ( min( t8(i,j,k), 1.0 ), 0. )
                  end do
                  end do
               end do
            end if
            call out_fstecr(t8,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                 'HR  ',Outd_convmult(pnhr,set),Outd_convadd(pnhr,set),&
                 knd,-1,G_nk+1,indo,nko,Outd_nbit(pnhr,set),.false. )
            if (write_diag_lev) then
               Out_stag_S(3:3)= 'D'
               call out_fstecr(t8(l_minx,l_miny,G_nk+1), &
                    l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2, &
                    'HR  ',Outd_convmult(pnhr,set),Outd_convadd(pnhr,set),&
                    Level_kind_diag,-1,1,ind0,1,Outd_nbit(pnhr,set),.false. )
               Out_stag_S(3:3)= ' '
            end if
         end if

         if (pnww /= 0) then
            call out_fstecr(myomega,l_minx,l_maxx,l_miny,l_maxy,hybt_w, &
                 'WW  ',Outd_convmult(pnww,set),Outd_convadd(pnww,set),&
                 knd,-1,G_nk,indo,nko,Outd_nbit(pnww,set),.false. )
         end if

         if (pnzz /= 0) then
            call out_fstecr(wt1(l_minx,l_miny,1),l_minx,l_maxx,l_miny,l_maxy,hybt_w, &
                 'ZZ  ',Outd_convmult(pnzz,set),Outd_convadd(pnzz,set),&
                 knd,-1,G_nk,indo,nko,Outd_nbit(pnzz,set),.false. )
         end if

         deallocate (indo)

         if (pnes /= 0.or.pnpx /= 0.or.pntw /= 0.or.pntd /= 0.or.pnhr /= 0) &
            deallocate (px_ta,px_m)
         if (pnes /= 0.or.pntw /= 0.or.pntd /= 0.or.pnhr /= 0) &
            deallocate (t8)

      else   ! Output on pressure levels

         Out_stag_S= 'MP '
         nko= Level_max(levset)
         allocate ( hu_pres(l_minx:l_maxx,l_miny:l_maxy,nko), &
                    vt_pres(l_minx:l_maxx,l_miny:l_maxy,nko), &
                    tt_pres(l_minx:l_maxx,l_miny:l_maxy,nko), &
                    td_pres(l_minx:l_maxx,l_miny:l_maxy,nko), &
                    px_pres(l_minx:l_maxx,l_miny:l_maxy,nko), &
                    w5     (l_minx:l_maxx,l_miny:l_maxy,nko), &
                    w6     (l_minx:l_maxx,l_miny:l_maxy,nko), &
                    cible  (l_minx:l_maxx,l_miny:l_maxy,nko), &
                    indo(nko), rf(nko) , prprlvl(nko) )

         knd=2 !for pressure output

         do i = 1, nko !Setup the indexing for output
            indo     (i)= i
            rf       (i)= Level(i,levset)
            prprlvl  (i)= rf(i) * 100.0
            cible(:,:,i)= log(prprlvl(i))
         end do

! Compute HU (hu_pres=HU,px_ta=vert.der)

         call vertint2 ( hu_pres,cible,nko, hu,wlnph_ta,nk_src,&
                         l_minx,l_maxx,l_miny,l_maxy          ,&
                    1,l_ni,1,l_nj, inttype=Out3_vinterp_type_S )

         if ( Out3_cliph_L ) then
            do k= 1, nko
               do j= 1,l_nj
                  do i= 1,l_ni
                     hu_pres(i,j,k) = amax1( hu_pres(i,j,k), 0. )
                  end do
               end do
            end do
         end if

! Compute GZ,VT (w5=GZ_pres, vt_pres=VT_pres)

         call prgzvta( w5, vt_pres, prprlvl, nko , &
                       gzt, vt, wlnph_ta, wlao   , &
                       ttx, htx, nk_under,.false., &
                       Out3_linbot, l_minx,l_maxx,l_miny,l_maxy,nk_src)

         call out_padbuf(vt_pres,l_minx,l_maxx,l_miny,l_maxy,nko)

        if (pngz /= 0) then
           if (Outd_filtpass(pngz,set) > 0)then
              call filter( w5,Outd_filtpass(pngz,set),Outd_filtcoef(pngz,set), &
                            l_minx,l_maxx,l_miny,l_maxy,nko)
           end if
           call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
              'GZ  ',Outd_convmult(pngz,set),Outd_convadd(pngz,set), &
              knd,-1,nko,indo,nko,Outd_nbit(pngz,set),.false. )
        end if

        if (pntt /= 0.or.pntd /= 0.or.pnhr /= 0) then

! Compute TT (tt_pres=TT,vt_pres=VT,hu_pres=HU)
           call mfottv2 (tt_pres,vt_pres,hu_pres,l_minx,l_maxx,l_miny,l_maxy, &
                         nko,1,l_ni,1,l_nj,.false.)
        end if

        if ( pnes /= 0.or.pntw /= 0.or.pntd /= 0.or.pnhr /= 0) then
! Compute PX for ES,TD,HR
            do k=1,nko
               do j= 1, l_nj
               do i= 1, l_ni
                  px_pres(i,j,k) = prprlvl(k)
               end do
               end do
            end do
            call out_padbuf(px_pres,l_minx,l_maxx,l_miny,l_maxy,nko)
            call out_padbuf(tt_pres,l_minx,l_maxx,l_miny,l_maxy,nko)
            call out_padbuf(hu_pres,l_minx,l_maxx,l_miny,l_maxy,nko)
        end if

        if (pntw /= 0) then
! Compute THETAW TW (w5=TW_pres) (px_pres=PX)
            call mfottv2 (w6,vt_pres,hu_pres,l_minx,l_maxx, &
                        l_miny,l_maxy,nko,1,l_ni,1,l_nj,.false.)
            call out_padbuf(w6,l_minx,l_maxx,l_miny,l_maxy,nko)
            call mthtaw4 (w5,hu_pres,w6, &
                           px_pres,satues_l, &
                           .true.,trpl_8,l_ninj,nko,l_ninj)
            if (Outd_filtpass(pntw,set) > 0) &
                call filter( w5,Outd_filtpass(pntw,set),Outd_filtcoef(pntw,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
            call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                'TW  ',Outd_convmult(pntw,set),Outd_convadd(pntw,set), &
                knd,-1,nko, indo, nko, Outd_nbit(pntw,set),.false. )
        end if

        if (pnes /= 0.or.pntd /= 0) then
! Compute ES (w5=ES_pres,hu_pres=HU,w2=VT,px_pres=PX)
            call mfottv2 (w6,vt_pres,hu_pres,l_minx,l_maxx, &
                        l_miny,l_maxy,nko,1,l_ni,1,l_nj,.false.)
            call out_padbuf(w6,l_minx,l_maxx,l_miny,l_maxy,nko)
            call mhuaes3 (w5, hu_pres,w6, px_pres,satues_l, &
                          l_ninj, nko, l_ninj)
            w5(1:l_ni,1:l_nj,1:nko) = min(w5(1:l_ni,1:l_nj,1:nko), ES_MAX)
            if ( Out3_cliph_L ) then
               w5(1:l_ni,1:l_nj,1:nko) = max(w5(1:l_ni,1:l_nj,1:nko),0.)
            end if

            if (pntd /= 0) then
! Compute TD (tt_pres=TT,w5=ES, TD=TT-ES)
              do k=1,nko
                 do j= 1, l_nj
                 do i= 1, l_ni
                    td_pres(i,j,k) = tt_pres(i,j,k) - w5(i,j,k)
                 end do
                 end do
              end do
              call filter( td_pres,Outd_filtpass(pntd,set),Outd_filtcoef(pntd,set), &
                            l_minx,l_maxx,l_miny,l_maxy, nko )
              call out_fstecr(td_pres,l_minx,l_maxx,l_miny,l_maxy,rf, &
                'TD  ',Outd_convmult(pntd,set),Outd_convadd(pntd,set),&
                knd,-1,nko,indo,nko,Outd_nbit(pntd,set),.false. )
            end if

            if (pnes /= 0) then
                if (Outd_filtpass(pnes,set) > 0) &
                    call filter( w5,Outd_filtpass(pnes,set),Outd_filtcoef(pnes,set), &
                                  l_minx,l_maxx,l_miny,l_maxy,nko )
                call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                   'ES  ',Outd_convmult(pnes,set),Outd_convadd(pnes,set),&
                   knd,-1,nko,indo,nko,Outd_nbit(pnes,set),.false.)
            end if
        end if

        if (pnhr /= 0) then
! Compute HR (w5=HR_pres:hu_pres=HU,tt_pres=TT,px_pres=PX)
           call mfohr4 (w5,hu_pres,tt_pres,px_pres,l_ninj,nko,l_ninj,satues_l)
           if ( Out3_cliph_L ) then
              do k=1,nko
                 do j= 1, l_nj
                    do i= 1, l_ni
                       w5(i,j,k) = min( w5(i,j,k), 1.0 )
                       w5(i,j,k) = max( w5(i,j,k), 0.  )
                    end do
                 end do
              end do
           end if
           if (Outd_filtpass(pnhr,set) > 0) &
                call filter( w5,Outd_filtpass(pnhr,set),Outd_filtcoef(pnhr,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
           call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                'HR  ',Outd_convmult(pnhr,set),Outd_convadd(pnhr,set), &
                knd,-1,nko, indo, nko, Outd_nbit(pnhr,set),.false. )
        end if

        if (pnvt /= 0) then
            if (Outd_filtpass(pnvt,set) > 0) &
                call filter( vt_pres,Outd_filtpass(pnvt,set),Outd_filtcoef(pnvt,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
            call out_fstecr(vt_pres,l_minx,l_maxx,l_miny,l_maxy,rf, &
                 'VT  ',Outd_convmult(pnvt,set),Outd_convadd(pnvt,set), &
                 knd,-1,nko,indo, nko, Outd_nbit(pnvt,set),.false. )
        end if

         if (pnth /= 0) then
            call vertint2 ( w5,cible,nko, th,wlnph_ta,G_nk+1          ,&
                            l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                           inttype=Out3_vinterp_type_S )
            call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                 'TH  ',Outd_convmult(pnth,set),Outd_convadd(pnth,set), &
                 knd,-1,nko, indo, nko, Outd_nbit(pnth,set),.false. )
         end if

        if (pntt /= 0) then
            if (Outd_filtpass(pntt,set) > 0) &
                call filter( tt_pres,Outd_filtpass(pntt,set),Outd_filtcoef(pntt,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
            call out_fstecr(tt_pres,l_minx,l_maxx,l_miny,l_maxy,rf,  &
                 'TT  ',Outd_convmult(pntt,set),Outd_convadd(pntt,set), &
                 knd,-1,nko, indo, nko, Outd_nbit(pntt,set),.false. )
        end if

        if (pnww /= 0) then
            call vertint2 ( w5,cible,nko, myomega,wlnph_ta,G_nk         ,&
                            l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                            inttype=Out3_vinterp_type_S )
            if (Outd_filtpass(pnww,set) > 0) &
                call filter( w5,Outd_filtpass(pnww,set),Outd_filtcoef(pnww,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
             call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                  'WW  ',Outd_convmult(pnww,set),Outd_convadd(pnww,set), &
                  knd,-1,nko, indo, nko, Outd_nbit(pnww,set),.false. )
        end if

         if (pnzz /= 0) then
           www(l_minx:,l_miny:,1:) => wt1(l_minx:,l_miny:,1:)
           call vertint2 ( w5,cible,nko, www,wlnph_ta,G_nk         ,&
                            l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj,&
                            inttype=Out3_vinterp_type_S )
            if (Outd_filtpass(pnzz,set) > 0) &
                call filter( w5,Outd_filtpass(pnzz,set),Outd_filtcoef(pnzz,set), &
                              l_minx,l_maxx,l_miny,l_maxy,nko )
             call out_fstecr(w5,l_minx,l_maxx,l_miny,l_maxy,rf, &
                  'ZZ  ',Outd_convmult(pnzz,set),Outd_convadd(pnzz,set),&
                  knd,-1,nko, indo, nko, Outd_nbit(pnzz,set),.false. )
        end if

        deallocate(indo,rf,prprlvl,cible)
        deallocate(w5,w6,px_pres,hu_pres,td_pres,tt_pres,vt_pres)
      end if

      if (pnww /= 0) deallocate (myomega)
      if (pnth /= 0) deallocate (th)
!
!-------------------------------------------------------------------
!
      return
      end subroutine out_thm
