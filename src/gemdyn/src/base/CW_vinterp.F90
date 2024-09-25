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

!**s/r CW_vinterp - Compute cubic weights for vertical interpolation

      subroutine CW_vinterp ( )
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
      allocate ( CWm2t(4,0:G_nk+1), CWt2m(4,0:G_nk+1) )
      CWm2t= 0.d0 ; CWt2m= 0.d0
      
      if(Dynamics_sw_L) return
      do k=2, l_nk-1 ! this assumes the data is valid at G_nk+1
         z1 = Ver_a_8%m(k-1)
         z2 = Ver_a_8%m(k  )
         z3 = Ver_a_8%m(k+1)
         z4 = Ver_a_8%m(k+2)
         pt = Ver_a_8%t(k  )
         CWm2t(1,k) = ((pt - z2) * (pt - z3) * (pt - z4) ) / ((z1 - z2) * (z1 - z3) * (z1 - z4) )
         CWm2t(2,k) = ((pt - z1) * (pt - z3) * (pt - z4) ) / ((z2 - z1) * (z2 - z3) * (z2 - z4) )
         CWm2t(3,k) = ((pt - z1) * (pt - z2) * (pt - z4) ) / ((z3 - z1) * (z3 - z2) * (z3 - z4) )
         CWm2t(4,k) = ((pt - z1) * (pt - z2) * (pt - z3) ) / ((z4 - z1) * (z4 - z2) * (z4 - z3) )
      end do
      CWm2t(2,   1) = 0.5 ; CWm2t(3,1   ) = 0.5
      CWm2t(2,G_nk) = 0.5 ; CWm2t(3,G_nk) = 0.5
      
      do k= 3, G_nk-1
         z1 = Ver_a_8%t(k-2)
         z2 = Ver_a_8%t(k-1)
         z3 = Ver_a_8%t(k  )
         z4 = Ver_a_8%t(k+1)
         pt = Ver_a_8%m(k  )
         CWt2m(1,k) = ((pt - z2) * (pt - z3) * (pt - z4) ) / ((z1 - z2) * (z1 - z3) * (z1 - z4) )
         CWt2m(2,k) = ((pt - z1) * (pt - z3) * (pt - z4) ) / ((z2 - z1) * (z2 - z3) * (z2 - z4) )
         CWt2m(3,k) = ((pt - z1) * (pt - z2) * (pt - z4) ) / ((z3 - z1) * (z3 - z2) * (z3 - z4) )
         CWt2m(4,k) = ((pt - z1) * (pt - z2) * (pt - z3) ) / ((z4 - z1) * (z4 - z2) * (z4 - z3) )
      end do
      CWt2m(2,1) = 1.d0 ! constant extrapolation temporary fudge
      CWt2m(2,2) = (Ver_z_8%t(2)-Ver_z_8%m(2))/(Ver_z_8%t(2)-Ver_z_8%t(1))
      CWt2m(3,2) = 1.d0 - CWt2m(2,2)
      CWt2m(2,G_nk) = (Ver_z_8%t(G_nk)-Ver_z_8%m(G_nk))/(Ver_z_8%t(G_nk)-Ver_z_8%t(G_nk-1))
      CWt2m(3,G_nk) = 1.d0 - CWt2m(2,G_nk)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine CW_vinterp
