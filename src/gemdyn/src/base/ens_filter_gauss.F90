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

!**s/r ens_filter_ggauss - Gaussian grid space filter
!
      subroutine ens_filter_gauss(dsp_lcl)
      use gem_options
      use ens_options
      use ens_param
      use glb_ld
      use gmm_vt1
      use, intrinsic :: iso_fortran_env
      implicit none
!
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: dsp_lcl
!
!author
!     Lubos Spacek - rpn - oct 2012

      integer i, j, k
      integer local_np, HLT_start, HLT_end

!!


      call HLT_split (1, G_nk, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( dsp_lcl(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )

!!$omp do
      do k=1,l_nk
    	 do j=l_miny,l_maxy
            do i=l_minx,l_maxx
               psi_local(i,j,k)=dsp_lcl(i,j,k)
            end do
         end do
      end do
!!$omp enddo

      ! x-axis
!!$omp do
      do k=1,l_nk
    	 do j=1,l_nj
            do i=1,l_ni
               dsp_lcl(i,j,k)= sum(fg1(j,:)*psi_local(i-3:i+3,j,k))
            end do
         end do
      end do
!!$omp enddo

      call gem_xch_halo ( dsp_lcl(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )

!!$omp do
      do k=1,l_nk
    	 do j=l_miny,l_maxy
            do i=l_minx,l_maxx
               psi_local(i,j,k)=dsp_lcl(i,j,k)
            end do
         end do
      end do
!!$omp enddo

      !y-axis
!!$omp do
      do k=1,l_nk
         do j=1,l_nj
            do i=1,l_ni
               dsp_lcl(i,j,k)= sum(fg2*psi_local(i,j-3:j+3,k))
            end do
         end do
      end do
!!$omp enddo

      end subroutine ens_filter_gauss
