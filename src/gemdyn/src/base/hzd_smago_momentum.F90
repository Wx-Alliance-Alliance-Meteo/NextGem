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

!**s/r hzd_smago_momentum   - applies horizontal diffusion based on the Smagorinsky approach

      subroutine hzd_smago_momentum()
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_vt0
      use hvdif_options
      use HORgrid_options
      use lun
      implicit none
!
!Author:  Syed Husain
!
      logical :: switch_on_wzd
      integer :: HLT_start, HLT_end, local_np
!
!-------------------------------------------------------------------
!
      if( (hzd_smago_param <= 0.) .and. (hzd_smago_lnr(1) <=0.) ) return

      switch_on_wzd   = (Hzd_lnr <= 0.)

      if ( (hzd_smago_param > 0.) .or. (hzd_smago_lnr(1) > 0.) ) then
         call HLT_split (-2, 5*(G_nk+6)-3, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( wt0(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call hzd_smago_in_split(ut0,vt0,wt0,tt0,zdt0, &
               l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)
      endif
            
      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut0,vt0,l_minx,l_maxx,l_miny,l_maxy,G_nk)

         if (switch_on_wzd) then
            call yyg_xchng_hlt (zdt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
            call yyg_xchng_hlt (wt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
         end if
      end if

 1000 format(3X,'SMAGO MOMENTUM DIFFUSION : (S/R HZD_SMAGO_MOMENTUM)')
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_smago_momentum
