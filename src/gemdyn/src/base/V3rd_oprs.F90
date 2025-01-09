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

!**s/r V3rd_oprs - Compute 3rd order weights for vertical operators

      subroutine V3rd_oprs ( )
      use, intrinsic :: iso_fortran_env
      use glb_ld
      use ver
      implicit none

      integer :: k
      real(kind=REAL64) :: p1, p2, p3, p4, l1, l2, l3, l4
      real(kind=REAL64) :: z1,z2,z3,z4,pt
!     
!     ---------------------------------------------------------------
!
      allocate ( VD3m2t(4,lbound(Ver_ext%m,1):ubound(Ver_ext%m,1)-1))
      allocate ( VS3m2t(4,lbound(Ver_ext%m,1):ubound(Ver_ext%m,1)-1))
      allocate ( VD3t2m(4,lbound(Ver_ext%t,1)+3:ubound(Ver_ext%t,1)-2) )
      allocate ( VS3t2m(4,lbound(Ver_ext%t,1)+3:ubound(Ver_ext%t,1)-2) )

      !===> momentum to thermo
      do k=0,G_nk
         z1 = Ver_ext%m(k-1)
         z2 = Ver_ext%m(k  )
         z3 = Ver_ext%m(k+1)
         z4 = Ver_ext%m(k+2)
         pt = Ver_ext%t(k  )
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3m2t(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3m2t(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3m2t(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3m2t(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3m2t(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3m2t(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3m2t(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3m2t(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
      do k=lbound(Ver_ext%m,1),-1
         z1 = Ver_ext%m(k  )
         z2 = Ver_ext%m(k+1)
         z3 = Ver_ext%m(k+2)
         z4 = Ver_ext%m(k+3)
         pt = Ver_ext%t(k  )
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3m2t(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3m2t(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3m2t(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3m2t(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3m2t(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3m2t(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3m2t(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3m2t(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
      do k=G_nk+1,ubound(Ver_ext%m,1)-1
         z1 = Ver_ext%m(k-2)
         z2 = Ver_ext%m(k-1)
         z3 = Ver_ext%m(k  )
         z4 = Ver_ext%m(k+1)
         pt = Ver_ext%t(k  )
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3m2t(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3m2t(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3m2t(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3m2t(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3m2t(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3m2t(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3m2t(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3m2t(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
      
      !===> thermo to momentum
      do k=1,G_nk
         z1 = Ver_ext%t(k-2)
         z2 = Ver_ext%t(k-1)
         z3 = Ver_ext%t(k  )
         z4 = Ver_ext%t(k+1)
         pt = Ver_ext%m(k  )
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3t2m(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3t2m(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3t2m(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3t2m(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3t2m(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3t2m(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3t2m(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3t2m(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
      do k=lbound(Ver_ext%t,1)+3,0
         z1 = Ver_ext%t(k-3)
         z2 = Ver_ext%t(k-2)
         z3 = Ver_ext%t(k-1)
         z4 = Ver_ext%t( k)
         pt = Ver_ext%m( k)
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3t2m(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3t2m(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3t2m(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3t2m(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3t2m(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3t2m(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3t2m(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3t2m(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
      do k=G_nk+1,ubound(Ver_ext%t,1)-2
         z1 = Ver_ext%t(k-1)
         z2 = Ver_ext%t(k  )
         z3 = Ver_ext%t(k+1)
         z4 = Ver_ext%t(k+2)
         pt = Ver_ext%m(k  )
         l1 = (z1 - z2)*(z1 - z3)*(z1 - z4)
         l2 = (z2 - z1)*(z2 - z3)*(z2 - z4)
         l3 = (z3 - z1)*(z3 - z2)*(z3 - z4)
         l4 = (z4 - z1)*(z4 - z2)*(z4 - z3)
         VD3t2m(1,k) = ((pt-z2)*(pt-z3) + (pt-z2)*(pt-z4) + (pt-z3)*(pt-z4)) / l1
         VD3t2m(2,k) = ((pt-z1)*(pt-z3) + (pt-z1)*(pt-z4) + (pt-z3)*(pt-z4)) / l2
         VD3t2m(3,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z4) + (pt-z2)*(pt-z4)) / l3
         VD3t2m(4,k) = ((pt-z1)*(pt-z2) + (pt-z1)*(pt-z3) + (pt-z2)*(pt-z3)) / l4
         VS3t2m(1,k) = ((pt-z2) * (pt-z3) * (pt-z4) ) / l1
         VS3t2m(2,k) = ((pt-z1) * (pt-z3) * (pt-z4) ) / l2
         VS3t2m(3,k) = ((pt-z1) * (pt-z2) * (pt-z4) ) / l3
         VS3t2m(4,k) = ((pt-z1) * (pt-z2) * (pt-z3) ) / l4
      end do
!!$      do k=lbound(VD3m2t,2), ubound(VD3m2t,2)
!!$         write(6,'(i3,4f20.12)') k,VD3m2t(:,k)
!!$      end do
!!$      do k=lbound(VD3m2t,2), ubound(VD3m2t,2)
!!$         write(6,'(i3,4f20.12)') k,VS3m2t(:,k)
!!$      end do
!!$      do k=lbound(VD3t2m,2), ubound(VD3t2m,2)
!!$         write(6,'(i3,4f20.12)') k,VD3t2m(:,k)
!!$      end do
!!$      do k=lbound(VD3t2m,2), ubound(VD3t2m,2)
!!$         write(6,'(i3,4f20.12)') k,VS3t2m(:,k)
!!$      end do
!!$      call gem_stop
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine V3rd_oprs
