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

!**s/r VS3_vinterp - Compute 3rd order weights for vertical interpolation

      subroutine VS3_vinterp ( )
      use, intrinsic :: iso_fortran_env
      use dynkernel_options
      use metric
      use glb_ld
      use ver
      implicit none

      integer :: k
      real(kind=REAL64) :: z1,z2,z3,z4,pt
!     
!     ---------------------------------------------------------------
!
      allocate ( VS3m2t(4,1:G_nk), VS3t2m(4,1:G_nk) )
      VS3m2t= 0.d0 ; VS3t2m= 0.d0
      
      if(Dynamics_sw_L) return
      do k=1, l_nk
         z1 = Ver_ext%m(k-1)
         z2 = Ver_ext%m(k  )
         z3 = Ver_ext%m(k+1)
         z4 = Ver_ext%m(k+2)
         pt = Ver_ext%t(k  )
         VS3m2t(1,k) = ((pt - z2) * (pt - z3) * (pt - z4) ) / ((z1 - z2) * (z1 - z3) * (z1 - z4) )
         VS3m2t(2,k) = ((pt - z1) * (pt - z3) * (pt - z4) ) / ((z2 - z1) * (z2 - z3) * (z2 - z4) )
         VS3m2t(3,k) = ((pt - z1) * (pt - z2) * (pt - z4) ) / ((z3 - z1) * (z3 - z2) * (z3 - z4) )
         VS3m2t(4,k) = ((pt - z1) * (pt - z2) * (pt - z3) ) / ((z4 - z1) * (z4 - z2) * (z4 - z3) )
      end do
      
      do k= 1, G_nk
         z1 = Ver_ext%t(k-2)
         z2 = Ver_ext%t(k-1)
         z3 = Ver_ext%t(k  )
         z4 = Ver_ext%t(k+1)
         pt = Ver_ext%m(k  )
         VS3t2m(1,k) = ((pt - z2) * (pt - z3) * (pt - z4) ) / ((z1 - z2) * (z1 - z3) * (z1 - z4) )
         VS3t2m(2,k) = ((pt - z1) * (pt - z3) * (pt - z4) ) / ((z2 - z1) * (z2 - z3) * (z2 - z4) )
         VS3t2m(3,k) = ((pt - z1) * (pt - z2) * (pt - z4) ) / ((z3 - z1) * (z3 - z2) * (z3 - z4) )
         VS3t2m(4,k) = ((pt - z1) * (pt - z2) * (pt - z3) ) / ((z4 - z1) * (z4 - z2) * (z4 - z3) )
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine VS3_vinterp
