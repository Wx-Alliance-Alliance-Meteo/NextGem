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

!**s/r out_dyn_casc - model output for cascade

      subroutine out_dyn_casc()
      use dyn_fisl_options
      use dynkernel_options
      use glb_ld
      use cstv
      use HORgrid_options
      use gmm_geof
      use rmn_gmm
      use gmm_pw
      use gmm_vt1
      use grdc_options
      use levels
      use metric
      use out3
      use svro_mod
      use out_meta
      use out_mod
      use out_options
      use out_vref
      use outp
      use tdpack
      use phy_itf, only: phy_get
      use vertical_interpolation
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_free,VGD_OK,VGD_ERROR
      use vgrid_wb, only: vgrid_wb_get
      implicit none

      character(len=512) varname
      logical, save :: done=.false.
      integer :: i,j,k, nbits, istat, indo(G_nk+2)
      integer, dimension(:), pointer  :: ip1m
      real, dimension(:    ), pointer :: hybm,hybt
      real, dimension(:,:,:), pointer :: ptr3d,tr2,wlnph_ta,wlnph_m

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,1:G_nk+1),target :: tr1

      type(vgrid_descriptor) :: vcoord
      real, dimension(1) :: hybm_gnk2, hybt_gnk2, hyb0
      integer, dimension(1) :: ind0
      real :: conv
!
!------------------------------------------------------------------
!     
      do k=1,G_nk+2
         indo(k) = k
      end do
      ind0(1)=1

      nullify (ip1m,hybm,hybt)
      istat = vgrid_wb_get('ref-m',vcoord,ip1m)
      deallocate(ip1m); nullify(ip1m)
      istat = vgd_get (vcoord,'VCDM - vertical coordinate (m)',hybm)
      istat = vgd_get (vcoord,'VCDT - vertical coordinate (t)',hybt)
      istat = vgd_free(vcoord)

      hyb0(1)=0.0
      hybm_gnk2(1)=hybm(G_nk+2)
      hybt_gnk2(1)=hybt(G_nk+2)

      Out_nfstecr = 0
      Out_reduc_l = .true.
      Out_stride  = 1            ! can only be one for now
      Out_gridi0  = max( 1   , Grdc_gid)
      Out_gridin  = min( G_ni, Grdc_gif)
      Out_gridj0  = max( 1   , Grdc_gjd)
      Out_gridjn  = min( G_nj, Grdc_gjf)
      Out_prefix_S= 'casc'
      
      if ( .not. OUTs_server_L) then
         call out_href ( 'Mass_point', Grdc_gid, Grdc_gif, 1,&
                                     Grdc_gjd, Grdc_gjf, 1 )

         call out_vref_itf ( etiket=Out3_etik_S )
      endif
      
      nbits= Grdc_nbits(1)
      if (done) nbits= Grdc_nbits(2)
      done= .true.
      
      call OUTs_metaS ()

      Out_stag_S= 'MT '
      conv = -tcdk_8
      call out_fstecr ( pw_tt_plus,l_minx,l_maxx,l_miny,l_maxy,hybt,&
                         'TT  ',1., conv,Level_kind_ip1,-1,G_nk,indo,&
                         G_nk,nbits,.false. )
      if (Out3_sfcdiag_L) then
         Out_stag_S(3:3)= 'D'
         call out_fstecr ( tdiag,l_minx,l_maxx,l_miny,l_maxy ,&
                            hybt_gnk2, 'TT  ', 1., conv,4,-1,1,&
                            ind0,1,nbits, .false. )
      end if
                         
      Out_stag_S= 'MM '
      conv = 1.d0 / knams_8
      call out_fstecr ( pw_uu_plus, l_minx,l_maxx,l_miny,l_maxy, hybm,&
                         'UU  ',conv, 0., Level_kind_ip1,-1,G_nk,indo ,&
                         G_nk,nbits,.false. )
      call out_fstecr ( pw_vv_plus, l_minx,l_maxx,l_miny,l_maxy, hybm,&
                         'VV  ',conv, 0., Level_kind_ip1,-1,G_nk,indo ,&
                         G_nk,nbits,.false. )
      if (Out3_sfcdiag_L) then
         Out_stag_S(3:3)= 'D'
         call out_fstecr ( udiag,l_minx,l_maxx,l_miny,l_maxy,hybm_gnk2,&
                    'UU  ' , conv, 0., 4,-1,1,ind0,1,nbits,.false. )
         call out_fstecr (vdiag, l_minx,l_maxx,l_miny,l_maxy,hybm_gnk2,&
                    'VV  ' , conv, 0., 4,-1,1,ind0,1,nbits,.false. )
      end if

      Out_stag_S= 'MS '
      call out_fstecr ( pw_p0_plus,l_minx,l_maxx,l_miny,l_maxy,hyb0,&
         'P0  ',.01, 0., 2,-1,1, ind0, 1, nbits, .false. )
         conv = 1.d0 / grav_8
      call out_fstecr(fis0,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
                  'ME  ', conv, 0.,2,-1,1, ind0, 1, nbits, .false. )
      if (Schm_sleve_L) then
         call out_fstecr(orols,l_minx,l_maxx,l_miny,l_maxy,hyb0, &
                  'MELS', conv ,0.,2,-1,1, ind0, 1, nbits, .false. )
      end if

      Out_stag_S= 'MT '
      call out_fstecr ( wt1 ,l_minx,l_maxx,l_miny,l_maxy, hybt,&
                         'WT1 ',1., 0.,Level_kind_ip1,-1,G_nk  ,&
                          indo,G_nk,nbits,.false. )
      call out_fstecr ( zdt1,l_minx,l_maxx,l_miny,l_maxy, hybt,&
                         'ZDT1',1., 0.,Level_kind_ip1,-1,G_nk  ,&
                         indo,G_nk,nbits,.false. )

      Out_stag_S= 'MM '
      call out_fstecr ( qt1,l_minx,l_maxx,l_miny,l_maxy, hybm,&
                       'QT1 ',1., 0.,Level_kind_ip1,-1,G_nk+1 ,&
                         indo,G_nk+1,nbits,.false. )

      do k=1,Grdc_ntr
         nullify (tr2)
         varname = 'TR/'//trim(Grdc_trnm_S(k))//':P'
         istat= gmm_get (varname,tr2)
         Out_stag_S= 'MT '
         call out_fstecr ( tr2 ,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                            Grdc_trnm_S(k),1.,0.,Level_kind_ip1,-1,&
                            G_nk, indo, G_nk, nbits,.false. )
         if ( Out3_sfcdiag_L ) then
            Out_stag_S(3:3)= 'D'
            if (trim(varname)=='TR/HU:P') then
               if (istat == 0) &
               call out_fstecr ( qdiag ,l_minx,l_maxx,l_miny,l_maxy, &
                                  hybt_gnk2,Grdc_trnm_S(k),1.,0.,4, &
                                  -1,1,ind0,1,nbits,.false. )
            else
               tr1(:,:,G_nk+1) = tr2(:,:,G_nk)
               ptr3d => tr1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:G_nk+1)
               istat = phy_get(ptr3d, varname, F_npath='VO', F_bpath='D')
               if (istat == 0) &
               call out_fstecr ( tr1(l_minx,l_miny,G_nk+1), &
                                l_minx,l_maxx,l_miny,l_maxy, &
                                hybt_gnk2,Grdc_trnm_S(k),1.,0.,4, &
                                -1,1,ind0,1,nbits,.false. )
            endif
         end if
      end do

      deallocate (hybm,hybt)
      call OUTs_metaF (Out_nfstecr, OUTs_nvar_indx)
!
!------------------------------------------------------------------
!
      return
      end

