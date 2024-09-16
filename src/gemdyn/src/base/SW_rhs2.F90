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
!	Compute bdf rhs terms
!
!**********************************************************************
!
      subroutine SW_rhs2 ( F_dt_8, i0, j0, k0, in, jn, k0t )

      use HORgrid_options
      use gem_options
      use step_options
      use dyn_fisl_options
      use coriolis
      use geomh
      use tdpack
      use mem_tstp
      use gmm_geof
      use gmm_vt1
      use gmm_vt2
      use glb_ld
      use cstv
      use dcst
      use ver
      use metric
      use mem_tstp
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
          integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
            HLT_np, HLT_start, HLT_end
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, km
      real(kind=REAL64) :: invT_n_8, w3
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!     __________________________________________________________________
!
      invT_n_8 = 3.0/(2.0*F_dt_8)
      w3 = half/(cpd_8*Cstv_Tstr_8)

!***********************************************************
! Compute the rhs terms that come from the bdf method      *
!***********************************************************


!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0-1, jn    !following the indices of nli terms
            do i= i0-1, in !following the indices of nli terms
              
              !---rhsu---
              rhsu(i,j,k) = (4.0/3.0)*invT_n_8*rhsu_mid(i,j,k) - (one/3.0)*invT_n_8*rhsu_dep(i,j,k) 

              !---rhsv---
              rhsv(i,j,k) = (4.0/3.0)*invT_n_8*rhsv_mid(i,j,k) - (one/3.0)*invT_n_8*rhsv_dep(i,j,k) 

              !---rhsc---
              rhsc(i,j,k) = (4.0/3.0)*invT_n_8*rhsc_mid(i,j,k) - (one/3.0)*invT_n_8*rhsc_dep(i,j,k) + invT_n_8*fis0(i,j)

            end do
         end do
      end do
!!$omp enddo

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsv(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsc(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

!     ---------------------------------------------------------------
!
      return
      end subroutine SW_rhs2
