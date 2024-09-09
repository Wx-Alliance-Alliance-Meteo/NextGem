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

!**s/r QW_vinterp - Compute quintic weights for vertical interpolation

      subroutine QW_vinterp ( )
      use, intrinsic :: iso_fortran_env
      use dynkernel_options
      use metric
      use glb_ld
      use ver
      implicit none

      integer :: k
      real(kind=REAL64) :: z1,z2,z3,z4,z5,z6,pt
!     
!     ---------------------------------------------------------------
!
      allocate ( QWm2t(6,0:G_nk+1), QWt2m(6,0:G_nk+1) )
      QWm2t= 0.d0 ; QWt2m= 0.d0
      
      if(Dynamics_sw_L) return

!      do k=3,G_nk-3
      do k=3,G_nk-2 ! this assumes the data is valid at G_nk+1
         z1 = Ver_a_8%m(k-2)
         z2 = Ver_a_8%m(k-1)
         z3 = Ver_a_8%m(k  )
         z4 = Ver_a_8%m(k+1)
         z5 = Ver_a_8%m(k+2)
         z6 = Ver_a_8%m(k+3)
         pt = Ver_a_8%t(k  )
         QWm2t(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6))
         QWm2t(2,k) = ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z2 - z1)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6))
         QWm2t(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z3 - z1)*(z3 - z2)*(z3 - z4)*(z3 - z5)*(z3 - z6))
         QWm2t(4,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / ((z4 - z1)*(z4 - z2)*(z4 - z3)*(z4 - z5)*(z4 - z6))
         QWm2t(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / ((z5 - z1)*(z5 - z2)*(z5 - z3)*(z5 - z4)*(z5 - z6))
         QWm2t(6,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / ((z6 - z1)*(z6 - z2)*(z6 - z3)*(z6 - z4)*(z6 - z5))
      end do
      z1 = Ver_a_8%m(1)
      z2 = Ver_a_8%m(2)
      z3 = Ver_a_8%m(3)
      z4 = Ver_a_8%m(4)      
      pt = Ver_a_8%t(2)
      QWm2t(2,2) = ((pt - z2)*(pt - z3)*(pt - z4)) / ((z1 - z2)*(z1 - z3)*(z1 - z4))
      QWm2t(3,2) = ((pt - z1)*(pt - z3)*(pt - z4)) / ((z2 - z1)*(z2 - z3)*(z2 - z4))
      QWm2t(4,2) = ((pt - z1)*(pt - z2)*(pt - z4)) / ((z3 - z1)*(z3 - z2)*(z3 - z4))
      QWm2t(5,2) = ((pt - z1)*(pt - z2)*(pt - z3)) / ((z4 - z1)*(z4 - z2)*(z4 - z3))
      QWm2t(3,1) = 0.5 ; QWm2t(4,1) = 0.5
!      z1 = Ver_a_8%m(G_nk-3)
!      z2 = Ver_a_8%m(G_nk-2)
!      z3 = Ver_a_8%m(G_nk-1)
!      z4 = Ver_a_8%m(G_nk  )      
!      pt = Ver_a_8%t(G_nk-2)
!      QWm2t(2,G_nk-2) = ((pt - z2)*(pt - z3)*(pt - z4)) / ((z1 - z2)*(z1 - z3)*(z1 - z4))
!      QWm2t(3,G_nk-2) = ((pt - z1)*(pt - z3)*(pt - z4)) / ((z2 - z1)*(z2 - z3)*(z2 - z4))
!      QWm2t(4,G_nk-2) = ((pt - z1)*(pt - z2)*(pt - z4)) / ((z3 - z1)*(z3 - z2)*(z3 - z4))
!      QWm2t(5,G_nk-2) = ((pt - z1)*(pt - z2)*(pt - z3)) / ((z4 - z1)*(z4 - z2)*(z4 - z3))
!      QWm2t(3,G_nk-1) = 0.5 ; QWm2t(4,G_nk-1) = 0.5
!      QWm2t(3,G_nk  ) = 0.5 ; QWm2t(4,G_nk  ) = 0.5 ! temporary
      z1 = Ver_a_8%m(G_nk-2)
      z2 = Ver_a_8%m(G_nk-1)
      z3 = Ver_a_8%m(G_nk  )
      z4 = Ver_a_8%m(G_nk+1)      
      pt = Ver_a_8%t(G_nk-1)
      QWm2t(2,G_nk-1) = ((pt - z2)*(pt - z3)*(pt - z4)) / ((z1 - z2)*(z1 - z3)*(z1 - z4))
      QWm2t(3,G_nk-1) = ((pt - z1)*(pt - z3)*(pt - z4)) / ((z2 - z1)*(z2 - z3)*(z2 - z4))
      QWm2t(4,G_nk-1) = ((pt - z1)*(pt - z2)*(pt - z4)) / ((z3 - z1)*(z3 - z2)*(z3 - z4))
      QWm2t(5,G_nk-1) = ((pt - z1)*(pt - z2)*(pt - z3)) / ((z4 - z1)*(z4 - z2)*(z4 - z3))
      QWm2t(3,G_nk  ) = 0.5 ; QWm2t(4,G_nk  ) = 0.5
      
      do k=4,G_nk-2
         z1 = Ver_a_8%t(k-3)
         z2 = Ver_a_8%t(k-2)
         z3 = Ver_a_8%t(k-1)
         z4 = Ver_a_8%t(k  )
         z5 = Ver_a_8%t(k+1)
         z6 = Ver_a_8%t(k+2)
         pt = Ver_a_8%m(k  )
         QWt2m(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6))
         QWt2m(2,k) = ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z2 - z1)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6))
         QWt2m(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / ((z3 - z1)*(z3 - z2)*(z3 - z4)*(z3 - z5)*(z3 - z6))
         QWt2m(4,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / ((z4 - z1)*(z4 - z2)*(z4 - z3)*(z4 - z5)*(z4 - z6))
         QWt2m(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / ((z5 - z1)*(z5 - z2)*(z5 - z3)*(z5 - z4)*(z5 - z6))
         QWt2m(6,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / ((z6 - z1)*(z6 - z2)*(z6 - z3)*(z6 - z4)*(z6 - z5))
      end do
      z1 = Ver_a_8%t(1)
      z2 = Ver_a_8%t(2)
      z3 = Ver_a_8%t(3)
      z4 = Ver_a_8%t(4)      
      pt = Ver_a_8%m(3)
      QWt2m(4,1) = 1.d0 ! constant extrapolation temporary fudge
      QWt2m(2,3) = ((pt - z2)*(pt - z3)*(pt - z4)) / ((z1 - z2)*(z1 - z3)*(z1 - z4))
      QWt2m(3,3) = ((pt - z1)*(pt - z3)*(pt - z4)) / ((z2 - z1)*(z2 - z3)*(z2 - z4))
      QWt2m(4,3) = ((pt - z1)*(pt - z2)*(pt - z4)) / ((z3 - z1)*(z3 - z2)*(z3 - z4))
      QWt2m(5,3) = ((pt - z1)*(pt - z2)*(pt - z3)) / ((z4 - z1)*(z4 - z2)*(z4 - z3))
      QWt2m(3,2) = (Ver_z_8%t(2)-Ver_z_8%m(2))/(Ver_z_8%t(2)-Ver_z_8%t(1))
      QWt2m(4,2) = 1.d0 - QWt2m(3,2)
      z1 = Ver_a_8%t(G_nk-3)
      z2 = Ver_a_8%t(G_nk-2)
      z3 = Ver_a_8%t(G_nk-1)
      z4 = Ver_a_8%t(G_nk  )      
      pt = Ver_a_8%m(G_nk-1)
      QWt2m(2,G_nk-1) = ((pt - z2)*(pt - z3)*(pt - z4)) / ((z1 - z2)*(z1 - z3)*(z1 - z4))
      QWt2m(3,G_nk-1) = ((pt - z1)*(pt - z3)*(pt - z4)) / ((z2 - z1)*(z2 - z3)*(z2 - z4))
      QWt2m(4,G_nk-1) = ((pt - z1)*(pt - z2)*(pt - z4)) / ((z3 - z1)*(z3 - z2)*(z3 - z4))
      QWt2m(5,G_nk-1) = ((pt - z1)*(pt - z2)*(pt - z3)) / ((z4 - z1)*(z4 - z2)*(z4 - z3))
      QWt2m(3,G_nk) = (Ver_z_8%t(G_nk)-Ver_z_8%m(G_nk))/(Ver_z_8%t(G_nk)-Ver_z_8%t(G_nk-1))
      QWt2m (4,G_nk) = 1.d0 - QWt2m(3,G_nk)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine QW_vinterp
