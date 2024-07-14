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

!**s/r nest_HOR_gwa

      subroutine nest_HOR_gwa ()
      use lam_options
      use theo_options
      use gmm_vt0
      use mem_nest
      use mem_tracers
      use glb_ld
      use tr3d
      implicit none

      integer i,j,k
      integer :: n,deb
!
!----------------------------------------------------------------------
!
      if ( Theo_periodicX_L .or. &
          ((Lam_blend_Hx <= 0).and.(Lam_blend_Hy <= 0)) ) return

!!$omp do
      do k=1,G_nk
      do j=1,l_nj
      do i=1,l_niu
         ut0(i,j,k) = ut0(i,j,k)*(1.-nest_weightu(i,j,k)) + nest_u(i,j,k)*nest_weightu(i,j,k)
      enddo
      enddo

      do j=1,l_njv
      do i=1,l_ni
         vt0(i,j,k) = vt0(i,j,k)*(1.-nest_weightv(i,j,k)) + nest_v(i,j,k)*nest_weightv(i,j,k)
      enddo
      enddo
      enddo
!!$omp end do nowait
!!$omp do collapse(2)
      do k=1,G_nk
      do j=1,l_nj
      do i=1,l_ni
         tt0(i,j,k) =  tt0(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_t (i,j,k)*nest_weightm(i,j,k)
         wt0(i,j,k) =  wt0(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_w (i,j,k)*nest_weightm(i,j,k)
         zdt0(i,j,k)= zdt0(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_zd(i,j,k)*nest_weightm(i,j,k)
      enddo
      enddo
      enddo
!!$omp end do nowait

!!$omp do collapse(2)
      do k=1,G_nk+1
      do j=1,l_nj
      do i=1,l_ni
         qt0(i,j,k) = qt0(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_q(i,j,k)*nest_weightm(i,j,k)
      enddo
      enddo
      enddo
!!$omp end do nowait

      do n=1,Tr3d_ntr
         deb= (n-1)*G_nk + 1
!!$omp do collapse(2)
         do k=1,G_nk
         do j=1,l_nj
         do i=1,l_ni
            tracers_M(n)%pntr(i,j,k) = &
            tracers_M(n)%pntr(i,j,k)*(1.-nest_weightm(i,j,k)) + &
                   nest_tr(i,j,deb+k-1)*nest_weightm(i,j,k)
         enddo
         enddo
         enddo
!!$omp end do nowait
      enddo

!!$OMP BARRIER

!
!----------------------------------------------------------------------
!
      return
      end subroutine nest_HOR_gwa
