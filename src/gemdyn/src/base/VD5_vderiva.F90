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

!**s/r VD5_vinterp - Compute 5th order weights for vertical derivatives

      subroutine VD5_vderiva ( )
      use, intrinsic :: iso_fortran_env
      use dynkernel_options
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
      allocate ( VD5m2t(6,1:G_nk), VD5t2m(6,1:G_nk) )
      VD5m2t= 0.d0 ; VD5t2m= 0.d0

      if(Dynamics_sw_L) return

      do k=1,G_nk
         z1 = Ver_ext%m(k-2)
         z2 = Ver_ext%m(k-1)
         z3 = Ver_ext%m(k  )
         z4 = Ver_ext%m(k+1)
         z5 = Ver_ext%m(k+2)
         z6 = Ver_ext%m(k+3)
         pt = Ver_ext%t(k  )

         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6)
         l2 = (z1 - z2)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6)
         l3 = (z1 - z3)*(z2 - z3)*(z3 - z4)*(z3 - z5)*(z3 - z6)
         l4 = (z1 - z4)*(z2 - z4)*(z3 - z4)*(z4 - z5)*(z4 - z6)
         l5 = (z1 - z5)*(z2 - z5)*(z3 - z5)*(z4 - z5)*(z5 - z6)
         l6 = (z1 - z6)*(z2 - z6)*(z3 - z6)*(z4 - z6)*(z5 - z6)

         VD5m2t(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l1

         VD5m2t(2,k) = - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l2

         VD5m2t(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l3

         VD5m2t(4,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l4

         VD5m2t(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l5

         VD5m2t(6,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l6
      end do

      do k=1,G_nk
         z1 = Ver_ext%t(k-3)
         z2 = Ver_ext%t(k-2)
         z3 = Ver_ext%t(k-1)
         z4 = Ver_ext%t(k  )
         z5 = Ver_ext%t(k+1)
         z6 = Ver_ext%t(k+2)
         pt = Ver_ext%m(k  )

         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)*(z1 - z5)*(z1 - z6)
         l2 = (z1 - z2)*(z2 - z3)*(z2 - z4)*(z2 - z5)*(z2 - z6)
         l3 = (z1 - z3)*(z2 - z3)*(z3 - z4)*(z3 - z5)*(z3 - z6)
         l4 = (z1 - z4)*(z2 - z4)*(z3 - z4)*(z4 - z5)*(z4 - z6)
         l5 = (z1 - z5)*(z2 - z5)*(z3 - z5)*(z4 - z5)*(z5 - z6)
         l6 = (z1 - z6)*(z2 - z6)*(z3 - z6)*(z4 - z6)*(z5 - z6)

         VD5t2m(1,k) = ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l1 &
         + ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l1

         VD5t2m(2,k) = - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l2 &
         - ((pt - z3)*(pt - z4)*(pt - z5)*(pt - z6)) / l2

         VD5t2m(3,k) = ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z1)*(pt - z4)*(pt - z5)*(pt - z6)) / l3 &
         + ((pt - z2)*(pt - z4)*(pt - z5)*(pt - z6)) / l3

         VD5t2m(4,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z2)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z1)*(pt - z3)*(pt - z5)*(pt - z6)) / l4 &
         - ((pt - z2)*(pt - z3)*(pt - z5)*(pt - z6)) / l4

         VD5t2m(5,k) = ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z6)) / l5 &
         + ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z6)) / l5

         VD5t2m(6,k) = - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z4)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z3)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z2)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z1)*(pt - z3)*(pt - z4)*(pt - z5)) / l6 &
         - ((pt - z2)*(pt - z3)*(pt - z4)*(pt - z5)) / l6
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine VD5_vderiva
