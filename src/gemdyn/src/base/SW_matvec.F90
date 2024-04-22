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

!** matvec_hlt - 3D Matrix-vector product subroutines (H coordinates)
      subroutine SW_matvec ( F_vector, F_minx,F_maxx,F_miny,F_maxy,&
                             F_prod  , F_i0,F_in,F_j0,F_jn, F_nk )
      use dyn_fisl_options
      use dynkernel_options
      use HORgrid_options
      use lam_options
      use glb_ld
      use metric
      use omp_timing
      use sol_mem
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_i0,F_in,F_j0,F_jn,F_nk
      real(kind=REAL64), dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_vector
      real(kind=REAL64), dimension(F_i0:F_in,F_j0:F_jn,F_nk), intent(out) :: F_prod

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, k0, k0t, km, kp, n
      integer :: i0,in,j0,jn
!
!     ---------------------------------------------------------------
!
      call gtmg_start (72, 'MATVEC1', 25 )
      i0 = 1+pil_w ; in = l_ni-pil_e
      j0 = 1+pil_s ; jn = l_nj-pil_n
      k0=1+Lam_gbpil_T
      k0t=k0 ; if (Schm_opentop_L) k0t=k0-1

!!$omp do collapse(2)
      do k= k0, l_nk
         do j=j0, jn
            do i=i0, in
               fdg2(i,j,k)= F_vector(i,j,k)
            end do
         end do
      end do
!!$omp enddo

!!$omp do
      do j= j0, jn
         do i= i0, in
            fdg2(i,j,l_nk+1)=F_vector(i,j,l_nk)
         end do
      end do
!!$omp enddo

      if ( Grd_yinyang_L) then
         call yyg_xchng_hlt (fdg2, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, l_nk+1, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, (l_nk+1), HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( fdg2(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy, HLT_np, 1)
      endif
      call gtmg_stop (72)
      call gtmg_start (73, 'MATVEC2', 25 )

!!$      HLT_j0 = j0
!!$      HLT_jn = jn
!!$      HLT_nk = l_nk
!!$      HLT_nj = HLT_jn - HLT_j0 + 1
!!$      call HLT_split (1, HLT_nj*HLT_nk, HLT_np, HLT_start, HLT_end)
!!$      do n= HLT_start, HLT_end
!!$         k= (n-1)/HLT_nj
!!$         j= n - k*HLT_nj + HLT_j0 - 1
!!$         k= k+1
!!$            km=max(k-1,1)
!!$            kp=k+1
!!$
!!$            do i= i0, in
!!$               F_prod(i,j,k)= &
!!$                  matv2_8 (i,j,k)*fdg2(i-1,j,  k  ) &
!!$               +  matv1_8 (i,j,k)*fdg2(i,  j,  k  ) &
!!$               +  matv3_8 (i,j,k)*fdg2(i+1,j,  k  ) &
!!$               +  matv6_8 (i,j,k)*fdg2(i-1,j,  km ) &
!!$               +  matv4_8 (i,j,k)*fdg2(i,  j,  km ) &
!!$               +  matv8_8 (i,j,k)*fdg2(i+1,j,  km ) &
!!$               +  matv7_8 (i,j,k)*fdg2(i-1,j,  kp ) &
!!$               +  matv5_8 (i,j,k)*fdg2(i,  j,  kp ) &
!!$               +  matv9_8 (i,j,k)*fdg2(i+1,j,  kp ) &
!!$               + matv10_8 (i,j,k)*fdg2(i,  j-1,k  ) &
!!$               + matv11_8 (i,j,k)*fdg2(i,  j+1,k  ) &
!!$               + matv12_8 (i,j,k)*fdg2(i,  j-1,km ) &
!!$               + matv13_8 (i,j,k)*fdg2(i,  j-1,kp ) &
!!$               + matv14_8 (i,j,k)*fdg2(i,  j+1,km ) &
!!$               + matv15_8 (i,j,k)*fdg2(i,  j+1,kp )
!!$            end do            
!!$      end do

!!$omp do collapse(2)
      do k=k0,l_nk
         do j= j0, jn
            km=max(k-1,1)
            kp=k+1
            do i= i0, in
               F_prod(i,j,k)= &
                  Sol_stencilh_8 (i,j,k, 1)*fdg2(i,  j,  k  ) &
               +  Sol_stencilh_8 (i,j,k, 2)*fdg2(i-1,j,  k  ) &
               +  Sol_stencilh_8 (i,j,k, 3)*fdg2(i+1,j,  k  ) &
               +  Sol_stencilh_8 (i,j,k,10)*fdg2(i,  j-1,k  ) &
               +  Sol_stencilh_8 (i,j,k,11)*fdg2(i,  j+1,k  ) 
            end do
         end do
      end do
!!$omp enddo
      call gtmg_stop (73)
!     
!     ---------------------------------------------------------------
!     
      return
      end subroutine SW_matvec
