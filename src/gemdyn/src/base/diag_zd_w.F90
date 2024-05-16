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

!** s/r - Computes diagnostic model vertical velocities zd and w

      subroutine diag_zd_w (F_zd, F_w, F_u, F_v, F_t, F_q, &
              F_metric,Minx, Maxx, Miny, Maxy, Nk, F_zd_L, F_w_L )
      use dyn_fisl_options
      use gem_options
      use geomh
      use glb_ld
      use lun
      use mem_tstp
      use metric
      use tdpack
      use ver
      implicit none

      integer, intent(in) ::  Minx, Maxx, Miny, Maxy, Nk
      logical, intent(in) ::  F_zd_L, F_w_L
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_zd, F_w
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q
      type(Vmetric) , intent(in ) :: F_metric

      integer :: i, j, k, kp, km, i0, in, j0, jn, dim
      real, dimension(:,:), pointer :: rJzX, rJzY, rJzZ
!     ________________________________________________________________
!
      dim=  (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      rJzX  (l_minx:l_maxx,l_miny:l_maxy) => WS1(      1:)
      rJzY  (l_minx:l_maxx,l_miny:l_maxy) => WS1(  dim+1:)
      rJzZ  (l_minx:l_maxx,l_miny:l_maxy) => WS1(2*dim+1:)

!     local grid setup for final results
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      if (l_west)  i0 = 2    - G_halox
      if (l_east)  in = in-1 + G_halox
      if (l_south) j0 = 2    - G_haloy
      if (l_north) jn = jn-1 + G_haloy

      if (F_zd_L) then

! Compute Zdot
!!$omp do collapse(2)
         do k=1, Nk
            do j= Miny, Maxy
               do i= Minx, Maxx
                  F_zd(i,j,k) = 0.
               end do
            end do
         end do
!!$omp enddo
!!$omp do
         do j= Miny, Maxy
            do i= Minx, Maxx
               rJzZ(i,j) = 0.
            end do
         end do
!!$omp enddo
!!$omp single
         do k=Nk,2,-1
            km=max(k-1,1)
            do j=j0,jn
               do i=i0-1,in
                  rJzX(i,j) = 0.5d0*((F_metric%ztht_8(i+1,j,k)-F_metric%ztht_8(i+1,j,km)) &
                            + (F_metric%ztht_8(i  ,j,k)-F_metric%ztht_8(i  ,j,km)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0-1,jn
               do i=i0,in
                  rJzY(i,j) = 0.5d0*((F_metric%ztht_8(i,j+1,k)-F_metric%ztht_8(i,j+1,km)) &
                            + (F_metric%ztht_8(i,j  ,k)-F_metric%ztht_8(i,j  ,km)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0,jn
               do i=i0,in
                  F_zd(i,j,k-1) = rJzZ(i,j)*F_zd(i,j,k) + Ver_dz_8%m(k)*( &
                                (rJzX(i  ,j)*F_u(i  ,j,k)                   &
                               -rJzX(i-1,j)*F_u(i-1,j,k))*geomh_invDX_8(j) &
                              +(rJzY(i,j  )*F_v(i,j  ,k)*geomh_cyV_8(j  )  &
                               -rJzY(i,j-1)*F_v(i,j-1,k)*geomh_cyV_8(j-1))*geomh_invDYM_8(j) )
                  rJzZ(i,j) = (F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k-1))*Ver_idz_8%t(k-1)
                  F_zd(i,j,k-1) = F_zd(i,j,k-1)/rJzZ(i,j)
               end do
            end do
         end do
!!$omp end single

      end if

      if(F_w_L) then

! Compute W (which depends on previous computation of Zdot)

         i0 = 1
         in = l_ni
         j0 = 1
         jn = l_nj
         if (l_west)  i0 = 3    - G_halox
         if (l_east)  in = in-1 + G_halox
         if (l_south) j0 = 3    - G_haloy
         if (l_north) jn = jn-1 + G_haloy
!!$omp do
         do k=1,Nk
            do j= Miny, Maxy
               do i= Minx, Maxx
                  F_w(i,j,k) = 0.
               end do
            end do
            km=max(k-1,1)
            kp=min(k+1,Nk)
            do j=j0,jn
               do i=i0,in
                  F_w(i,j,k) = 0.25d0* ( &
                     (F_u(i,j,kp)*F_metric%mc_Jx_8(i,j,kp)+F_u(i-1,j,kp)*F_metric%mc_Jx_8(i-1,j,kp))   &
                    +(F_u(i,j,k )*F_metric%mc_Jx_8(i,j,k )+F_u(i-1,j,k )*F_metric%mc_Jx_8(i-1,j,k ))   &
                    +(F_v(i,j,kp)*F_metric%mc_Jy_8(i,j,kp)+F_v(i,j-1,kp)*F_metric%mc_Jy_8(i,j-1,kp))   &
                    +(F_v(i,j,k )*F_metric%mc_Jy_8(i,j,k )+F_v(i,j-1,k )*F_metric%mc_Jy_8(i,j-1,k )) ) &
                    +(Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km))  &
                    *(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))*Ver_idz_8%t(k)
               end do
            end do
         end do
!!$omp end do
      end if
!
!     ________________________________________________________________
!
      return
      end subroutine diag_zd_w
