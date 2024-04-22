!COMP_ARCH=intel13sp1u2 ; -suppress=-C
!COMP_ARCH=intel-2016.1.156; -suppress=-C

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

      subroutine itf_ens_hzd () 
      use gem_options
      use phy_itf
      use ens_gmm_var
      use ens_options
      use ens_param
      use mem_tstp
      use HORgrid_options
      use glb_ld
      use gmm_vt1
      use cstv
      use ens_spp
      use omp_lib
      use ptopo
      implicit none

      integer k,istat,iend(3),dim
      integer :: HLT_start, HLT_end, local_np
      real, dimension(:,:,:), pointer :: ptr_3d,ug, vg, ug_s, vg_s

!-------------------------------------------------------------------

      if(.not. ens_conf) return

      call ens_markov_spp ()

      if ( .not.  Ens_skeb_conf)  return
      call ens_markov_skeb ()

      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
      ug   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(      1:)
      vg   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(  dim+1:)
      ug_s (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(2*dim+1:)
      vg_s (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(3*dim+1:)

      if(Ens_skeb_gwd)then

         iend = (/-1,-1,l_nk/)
         ptr_3d => ugwdt1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
         istat = phy_get(ptr_3d,'phytd_udis',F_npath='V',F_bpath='P',F_end=iend)
         ptr_3d => vgwdt1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,:)
         istat = phy_get(ptr_3d,'phytd_vdis',F_npath='V',F_bpath='P',F_end=iend)
!!$omp do
         do k= 1, G_nk
            ug(l_minx:l_maxx, l_miny:Grd_lphy_j0-1, k) = 0. ; ug(l_minx:l_maxx,Grd_lphy_jn+1:l_maxy,k) = 0.
            vg(l_minx:l_maxx, l_miny:Grd_lphy_j0-1, k) = 0. ; vg(l_minx:l_maxx,Grd_lphy_jn+1:l_maxy,k) = 0.
            ug(l_minx:Grd_lphy_i0-1,  l_miny:l_maxy,k) = 0. ; ug(Grd_lphy_in+1:l_maxx,l_miny:l_maxy,k) = 0.
            vg(l_minx:Grd_lphy_i0-1,  l_miny:l_maxy,k) = 0. ; vg(Grd_lphy_in+1:l_maxx,l_miny:l_maxy,k) = 0.
            ug(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,k) = ugwdt1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,k)
            vg(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,k) = vgwdt1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,k)
            ug_s(:,:,k)= 0. ; vg_s(:,:,k)= 0.
         end do
!!$omp end do

!Stager ug, vg
         call HLT_split (1, l_nk, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( ug(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
         call gem_xch_halo ( vg(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
         call hwnd_stag_hlt ( ug_s, vg_s, ug, vg ,l_minx,l_maxx,l_miny,l_maxy,l_nk,.true. )

!!$omp do
         do k= 1, G_nk
            ugwdt1(:,:,k) = ug_s(:,:,k) * Cstv_dt_8
            vgwdt1(:,:,k) = vg_s(:,:,k) * Cstv_dt_8
         end do
!!$omp enddo
         end if

!!$omp do
         do k= 1, G_nk
            difut1(:,:,k) = ut1(:,:,k) - difut1(:,:,k)
            difvt1(:,:,k) = vt1(:,:,k) - difvt1(:,:,k)
         end do
!!$omp enddo

         call ens_skeb_apply (ugwdt1,vgwdt1,difut1,difvt1)

           
!
!
!-------------------------------------------------------------------
!
      return
      end subroutine itf_ens_hzd 
