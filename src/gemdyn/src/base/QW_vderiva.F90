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

!**s/r QW_vinterp - Compute quintic weights for vertical derivatives

      subroutine QW_vderiva ( )
      use, intrinsic :: iso_fortran_env
      use metric
      use glb_ld
      use ver
      implicit none

      integer :: k
      real(kind=REAL64) :: p1, p2, p3, p4, p5, p6, l1, l2, l3, l4, l5, l6
      real(kind=REAL64) :: z1,z2,z3,z4,z5,z6,pt
!     
!     ---------------------------------------------------------------
!
      allocate ( QDm2t(6,0:G_nk+1), QDt2m(6,0:G_nk+1) )
      QDm2t= 0.d0 ; QDt2m= 0.d0

!      do k=3,G_nk-3
      do k=3,G_nk-2 ! this assumes the data is valid at G_nk+1
         z1 = Ver_a_8%m(k-2)
         z2 = Ver_a_8%m(k-1)
         z3 = Ver_a_8%m(k  )
         z4 = Ver_a_8%m(k+1)
         z5 = Ver_a_8%m(k+2)
         z6 = Ver_a_8%m(k+3)
         pt = Ver_a_8%t(k  )

         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6)
         l2 = (z1 - z2)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6)
         l3 = (z1 - z3)*(z2 - z3)*(z3 - z4)*(z3 - z5)*(z3 - z6)
         l4 = (z1 - z4)*(z2 - z4)*(z3 - z4)*(z4 - z5)*(z4 - z6)
         l5 = (z1 - z5)*(z2 - z5)*(z3 - z5)*(z4 - z5)*(z5 - z6)
         l6 = (z1 - z6)*(z2 - z6)*(z3 - z6)*(z4 - z6)*(z5 - z6)

         QDm2t(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l1

         QDm2t(2,k) = - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l2

         QDm2t(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l3

         QDm2t(4,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l4

         QDm2t(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l5

         QDm2t(6,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l6
      end do

      QDm2t(4,1) = 1.d0/(Ver_z_8%m(2) - Ver_z_8%m(1))
      QDm2t(4,2) = 1.d0/(Ver_z_8%m(3) - Ver_z_8%m(2))
!      QDm2t(4,G_nk-2) = 1.d0/(Ver_z_8%m(G_nk-1) - Ver_z_8%m(G_nk-2))
      QDm2t(4,G_nk-1) = 1.d0/(Ver_z_8%m(G_nk  ) - Ver_z_8%m(G_nk-1))
      QDm2t(4,G_nk  ) = 1.d0/(Ver_z_8%m(G_nk+1) - Ver_z_8%m(G_nk  ))
      QDm2t(3,1) = -QDm2t(4,1)
      QDm2t(3,2) = -QDm2t(4,2)
!      QDm2t(3,G_nk-2) = -QDm2t(4,G_nk-2)
      QDm2t(3,G_nk-1) = -QDm2t(4,G_nk-1)
      QDm2t(3,G_nk  ) = -QDm2t(4,G_nk  )

      do k=4,G_nk-2
         z1 = Ver_a_8%t(k-3)
         z2 = Ver_a_8%t(k-2)
         z3 = Ver_a_8%t(k-1)
         z4 = Ver_a_8%t(k  )
         z5 = Ver_a_8%t(k+1)
         z6 = Ver_a_8%t(k+2)
         pt = Ver_a_8%m(k  )

         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6)
         l2 = (z1 - z2)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6)
         l3 = (z1 - z3)*(z2 - z3)*(z3 - z4)*(z3 - z5)*(z3 - z6)
         l4 = (z1 - z4)*(z2 - z4)*(z3 - z4)*(z4 - z5)*(z4 - z6)
         l5 = (z1 - z5)*(z2 - z5)*(z3 - z5)*(z4 - z5)*(z5 - z6)
         l6 = (z1 - z6)*(z2 - z6)*(z3 - z6)*(z4 - z6)*(z5 - z6)

         QDt2m(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l1

         QDt2m(2,k) = - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l2

         QDt2m(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l3

         QDt2m(4,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l4

         QDt2m(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l5

         QDt2m(6,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l6
      end do
      QDt2m(4,1) = 1.d0/(Ver_z_8%t(1) - Ver_z_8%t(0))
      QDt2m(4,2) = 1.d0/(Ver_z_8%t(2) - Ver_z_8%t(1))
      QDt2m(4,3) = 1.d0/(Ver_z_8%t(3) - Ver_z_8%t(2))
      QDt2m(4,G_nk-1) = 1.d0/(Ver_z_8%t(G_nk-1) - Ver_z_8%t(G_nk-2))
      QDt2m(4,G_nk  ) = 1.d0/(Ver_z_8%t(G_nk  ) - Ver_z_8%t(G_nk-1))
      QDt2m(3,1) = -QDt2m(4,1)
      QDt2m(3,2) = -QDt2m(4,2)
      QDt2m(3,3) = -QDt2m(4,3)
      QDt2m(3,G_nk-1) = -QDt2m(4,G_nk-1)
      QDt2m(3,G_nk  ) = -QDt2m(4,G_nk  )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine QW_vderiva
