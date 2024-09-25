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

!**s/r heights - Compute heights

      subroutine heights ()
      use, intrinsic :: iso_fortran_env
      use HORgrid_options
      use metric
      use gmm_geof
      use yyg_param
      implicit none

      integer :: i,j,k
      integer :: HLT_start, HLT_end, local_np
      real, dimension (:,:,:), pointer :: oro
!
!     ---------------------------------------------------------------
!
      if (Grd_yinyang_L) then
         call yyg_xchng_hlt (fis0, l_minx,l_maxx,l_miny,l_maxy, &
                             l_ni,l_nj, 6, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, 6, local_np, HLT_start, HLT_end)
         oro(l_minx:l_maxx,l_miny:l_maxy,1:6) => orography(1:)
         call gem_xch_halo (oro(l_minx,l_minx,HLT_start), l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      endif
      
      call lvl_heights ( GVM%zmom_8, GVM%ztht_8, &
                         fis0, orols, l_minx,l_maxx,l_miny,l_maxy)
      call lvl_heights ( zmomu_8, zthtu_8, &
                         fis0u, orolsu, l_minx,l_maxx,l_miny,l_maxy)
      call lvl_heights ( zmomv_8, zthtv_8, &
                         fis0v, orolsv, l_minx,l_maxx,l_miny,l_maxy)
      do k=1,G_nk
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               dgzm(i,j,k)=GVM%zmom_8(i,j,k)-GVM%zmom_8(i,j,k+1)
               dgzt(i,j,k)=GVM%ztht_8(i,j,k)-GVM%ztht_8(i,j,k+1)
            end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine heights
