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
!**s/r hzd_lid_sponge

      subroutine hzd_lid_sponge ()
      use hzd_mod
      use hvdif_options
      use gem_options
      use HORgrid_options
      use gmm_vt1
      use glb_ld
      use glb_pil
      use, intrinsic :: iso_fortran_env
      implicit none

      integer :: iter, i0,in,inu,j0,jn,inv,jnv, available
      integer :: HLT_start2, lcl2, HLT_start5, lcl5, HLT_end
!
!     ---------------------------------------------------------------
!
      i0  = 2         - G_halox * (1 - west )
      j0  = 2         - G_haloy * (1 - south)
      in  = l_ni  - 1 + G_halox * (1 - east )
      inu = l_niu - 1 + G_halox * (1 - east )
      jn  = l_nj  - 1 + G_haloy * (1 - north)
      inv = l_ni  - 1 + G_halox * (1 - east )
      jnv = l_njv - 1 + G_haloy * (1 - north)
      available= 0
      
      call HLT_split (-2, 2*(G_nk+6)-3, lcl2, HLT_start2, HLT_end)
      call HLT_split (-2, 5*(G_nk+6)-3, lcl5, HLT_start5, HLT_end)

      do iter= 1, Vspng_niter
         if (available<1)   then
            if (Grd_yinyang_L) then
            call yyg_xchng_vec_uv2uv (ut1(l_minx,l_miny,1), vt1(l_minx,l_miny,1),l_minx,l_maxx,l_miny,l_maxy,G_nk)
            call gem_xch_halo (ut1(l_minx,l_miny,HLT_start2),&
                               l_minx,l_maxx,l_miny,l_maxy,lcl2,-1)
            call yyg_xchng_hlt (tt1(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,&
                            Vspng_nk, .false., 'CUBIC', .true. )            
            call yyg_xchng_hlt (wt1(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,&
                            Vspng_nk, .false., 'CUBIC', .true. )            
            call yyg_xchng_hlt (zdt1(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,&
                            Vspng_nk, .false., 'CUBIC', .true. )            
            else
               call gem_xch_halo (wt1(l_minx,l_miny,HLT_start5),&
                             l_minx,l_maxx,l_miny,l_maxy,lcl5,-1)
            endif
            available=G_halox
         endif            
         call hzd_flt5pt ( ut1(l_minx,l_miny,1), Hzd_geom_u, l_minx,l_maxx,l_miny,l_maxy,&
                           Vspng_nk, Vspng_coef_8, i0,inu,j0,jn )
         call hzd_flt5pt ( vt1(l_minx,l_miny,1), Hzd_geom_v, l_minx,l_maxx,l_miny,l_maxy,&
                           Vspng_nk, Vspng_coef_8, i0,inv,j0,jnv )
         call hzd_flt5pt ( tt1(l_minx,l_miny,1), Hzd_geom_q, l_minx,l_maxx,l_miny,l_maxy,&
                           Vspng_nk, Vspng_coef_8, i0,in,j0,jn )
         call hzd_flt5pt ( wt1(l_minx,l_miny,1), Hzd_geom_q, l_minx,l_maxx,l_miny,l_maxy,&
                           Vspng_nk, Vspng_coef_8, i0,in,j0,jn )
         call hzd_flt5pt ( zdt1(l_minx,l_miny,1), Hzd_geom_q, l_minx,l_maxx,l_miny,l_maxy,&
                           Vspng_nk, Vspng_coef_8, i0,in,j0,jn )
         available=available-1
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine hzd_lid_sponge
