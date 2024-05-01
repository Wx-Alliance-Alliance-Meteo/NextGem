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

!**s/r adz_wnds_ext - extended arrays for advection winds

      subroutine adz_wnds_ext ()
      use gmm_vt1
      use gmm_vt2
      use gmm_pw
      use adz_mem
      implicit none

      integer i,j,k
      integer :: HLT_start, HLT_end, local_np
!
!     ---------------------------------------------------------------
!
!!$omp do collapse(2)
      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               Adz_uu_ext(i,j,k) = pw_uu_moins(i,j,k)
               Adz_vv_ext(i,j,k) = pw_vv_moins(i,j,k)
               Adz_ww_ext(i,j,k) =        zdt1(i,j,k)
            end do
         end do
      end do
!!$omp enddo

      call HLT_split (1, 3*l_nk, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( Adz_uu_ext(Adz_lminx,Adz_lminy,HLT_start),&
                Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, local_np,-1) 
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_wnds_ext
