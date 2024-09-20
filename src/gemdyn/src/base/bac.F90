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
!---------------------------------- LICENCE END --------------------------------

!**s/r bac - backsubstitution: obtain new values for variables:u,v,w,t,q,zd
!                                using Sol_lhs

      subroutine bac ( F_dt_8 )
      use dyn_fisl_options
      use geomh
      use sol_mem
      use HORgrid_options
      use tdpack
      use gmm_vt0
      use mem_tstp
      use glb_ld
      use cstv
      use ver
      use metric
      use yyg_param
      use ctrl
      use theo_options
      use glb_pil
      use ldnh
      use stat_mpi
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k,ni,nj
      integer :: HLT_start, HLT_end, local_np
      real(kind=REAL64) :: w5, tau_8, invT_8, Buoy
      real(kind=REAL64), parameter :: one=1.d0
!
!     ---------------------------------------------------------------
!
      if (Ctrl_testcases_adv_L) then
!!$omp single
         call canonical_cases ("BAC")
!!$omp end single
         return
      end if

      tau_8 = (2.d0*F_dt_8) / 3.d0 
      invT_8= one/tau_8

!--extrapolation at the surface and at the lid (k=0 and k=l_nk+1)   
!!$omp do
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            ! pour le moment cette extrapolation ne fonctionne pas
            ! w5=(GVM%zmom_8(i,j,l_nk+1)-GVM%zmom_8(i,j,l_nk))/(GVM%zmom_8(i,j,l_nk)-GVM%zmom_8(i,j,l_nk-1))
            ! Sol_lhs(i,j,l_nk+1)= Sol_lhs(i,j,l_nk)*(1+w5) - w5*Sol_lhs(i,j,l_nk-1)
            ! so we keep the following (from FISL)
            Sol_lhs(i,j,l_nk+1) = GVM%mc_alfas_H_8(i,j)*Sol_lhs(i,j,l_nk  ) &
                                - GVM%mc_betas_H_8(i,j)*Sol_lhs(i,j,l_nk-1) &
                              + GVM%mc_css_H_8  (i,j)* (Rtt(i,j,l_nk)-Ver_wmstar_8(G_nk)*Rtt(i,j,l_nk-1) &
                              + invT_8*(Rzz(i,j,l_nk)-Ver_wmstar_8(G_nk)*Rzz(i,j,l_nk-1)))
            w5= (ver_z_8%m(0)-GVM%zmom_8(i,j,2))/(GVM%zmom_8(i,j,2)-GVM%zmom_8(i,j,1))
            Sol_lhs(i,j,0)= Sol_lhs(i,j,2)*(1+w5) -w5*Sol_lhs(i,j,1)
         enddo
      enddo
!!$omp end do

      call delQ (Sol_lhs, l_minx,l_maxx,l_miny,l_maxy, Qu,Qv,Qw,Qq,0,l_nk+1)

!!$omp do collapse(2)
      do k=1, l_nk+1
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               qt0(i,j,k) = sngl(Sol_lhs(i,j,k))
            end do
         end do
      end do
!!$omp enddo

!!$omp do collapse(2)
      do k=ds_k0, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, l_niu-pil_e
               ut0(i,j,k) = tau_8*(Ruu(i,j,k) - Qu(i,j,k))
            end do
         end do
         do j= ds_j0, l_njv-pil_n
            do i= ds_i0, ds_in
               vt0(i,j,k) = tau_8*(Rvv(i,j,k) - Qv(i,j,k))
            end do
         end do
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               wt0 (i,j,k) = tau_8*(Rtt(i,j,k) - gama_bdf_8*Qw(i,j,k))
               zdt0(i,j,k) = (Rzz(i,j,k) + wt0(i,j,k))
               Buoy = Qq(i,j,k) + wt0(i,j,k)*invT_8 - Rww(i,j,k)
               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy / grav_8 )
! or alternatively
!!$               a = 4.d0*invT_8/3.d0
!!$               b =      invT_8/3.d0
!!$               w5= a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k)&
!!$                   - invT_8* ( log(tt0(i,j,k)/Cstv_Tstr_8) - (one-Cstv_Tstr_8/tt0(i,j,k) ))
!!$               Buoy = half*(Sol_lhs(i,j,k+1)+Sol_lhs(i,j,k))/(cpd_8*Cstv_Tstr_8) + tau_8*(w5-mu_8*wt0(i,j,k))
!!$               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy)
            end do
         end do
      end do
!!$omp enddo nowait
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine bac
