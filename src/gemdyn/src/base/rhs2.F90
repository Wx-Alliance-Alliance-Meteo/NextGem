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
!**s/r rhs2 - Compute bdf rhs terms

      subroutine rhs2 ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use gem_options
      use adz_mem
      use mem_tstp
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: invT_8
      real(kind=REAL64), parameter :: one=1.d0
!     __________________________________________________________________
!
      invT_8 = 3.0/(2.0*F_dt_8)

!***********************************************************
! Compute the rhs terms that come from the bdf method      *
!***********************************************************

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= Adz_j0, Adz_jn
            do i= Adz_i0u, Adz_inu
              rhsu(i,j,k) = (4.0/3.0)*invT_8*rhsu_mid(i,j,k) - (one/3.0)*invT_8*rhsu_dep(i,j,k) 
            end do
         end do
      end do
!!$omp enddo
!!$omp do collapse(2)
      do k=k0, l_nk
         do j= Adz_j0v, Adz_jnv
            do i= Adz_i0, Adz_in
              rhsv(i,j,k) = (4.0/3.0)*invT_8*rhsv_mid(i,j,k) - (one/3.0)*invT_8*rhsv_dep(i,j,k) 
            end do
         end do
      end do
!!$omp enddo
!!$omp do collapse(2)
     do k=k0, l_nk
         do j= Adz_j0, Adz_jn
            do i= Adz_i0, Adz_in
              rhsc(i,j,k) = (4.0/3.0)*invT_8*rhsc_mid(i,j,k) - (one/3.0)*invT_8*rhsc_dep(i,j,k) 
            end do
         end do
      end do
!!$omp enddo
!!$omp do collapse(2)
      do k=k0t,l_nk
         do j= Adz_j0, Adz_jn
            do i= Adz_i0, Adz_in
              rhsw(i,j,k) = (4.0/3.0)*invT_8*rhsw_mid(i,j,k) - (one/3.0)*invT_8*rhsw_dep(i,j,k) 
              rhst(i,j,k) = (4.0/3.0)*invT_8*rhst_mid(i,j,k) - (one/3.0)*invT_8*rhst_dep(i,j,k)
              rhsf(i,j,k) = (4.0/3.0)*invT_8*rhsf_mid(i,j,k) - (one/3.0)*invT_8*rhsf_dep(i,j,k) 
             ! true_rhst(i,j,k) = rhst(i,j,k)
            end do
         end do
      end do
!!$omp enddo
      call HLT_split (1, 2*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine rhs2
