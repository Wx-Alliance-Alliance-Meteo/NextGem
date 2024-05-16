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

!**s/r  bac_H- backsubstitution: obtain new values of the variables:u,v,w,t,q,zd
!                from new q , the right-hand sides (Ru,Rv,Rw,Rt,Rf)
!                             and non-linear terms (Nu,Nv,Nw,Nt)

      subroutine bac ( F_dt_8, i0, j0, k0, in, jn ,k0t )
      use gem_options
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use sol_mem
      use HORgrid_options
      use tdpack
      use gmm_vt0
      use gmm_vt1
      use mem_tstp
      use glb_pil
      use glb_ld
      use cstv
      use ver
      use metric
      use ctrl
      use step_options 
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0, j0, k0, in, jn, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k
      integer :: HLT_start, HLT_end, local_np
      real :: w3,w5,w6,q_ext(l_minx:l_maxx,l_miny:l_maxy,0:l_nk+1)
      real(kind=REAL64) :: tau_m_8, tau_nh_8, invT_m_8, invT_nh_8,&
                           Buoy, wp_8(G_nk  ), wm_8(G_nk  )
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      if (Ctrl_testcases_adv_L) then
!!$omp single
         call canonical_cases ("BAC")
!!$omp end single
         return
      end if
      
      do k=1,G_nk
         wm_8(k) = (Ver_dqdz_8(k)-Ver_z_8%m(k))/&
                   (Ver_dqdz_8(k)-Ver_dqdz_8(k-1))
         wp_8(k) = one-wm_8(k)
      end do

      w3       = 1.d0/(cpd_8*Cstv_Tstr_8)
      tau_m_8  = (2.d0*F_dt_8) / 3.d0 
      tau_nh_8 = (2.d0*F_dt_8) / 3.d0 
      invT_m_8 = one/tau_m_8
      invT_nh_8= one/tau_nh_8

!!$omp do collapse(2)
      do k= k0, l_nk
         do j= j0, jn
            do i= i0, in
               q_ext(i,j,k) = sngl(Sol_lhs(i,j,k))
            end do
         end do
      end do
!!$omp enddo

!--extrapolation at the surface    
!!$omp do
      do j= j0, jn
         do i= i0, in
            w5=(GVM%zmom_8(i,j,l_nk+1)-GVM%zmom_8(i,j,l_nk))/(GVM%zmom_8(i,j,l_nk)-GVM%zmom_8(i,j,l_nk-1))
            ! pour le moment cette extrapolation ne fonctionne pas
!         q_ext(i,j,l_nk+1)= Sol_lhs(i,j,l_nk)*(1+w5) - w5*Sol_lhs(i,j,l_nk-1)
            q_ext(i,j,l_nk+1) = GVM%mc_alfas_H_8(i,j) * q_ext(i,j,l_nk)   &
                         - GVM%mc_betas_H_8(i,j) * q_ext(i,j,l_nk-1) &
                         + GVM%mc_css_H_8(i,j)   * (rhst(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhst(i,j,l_nk-1) &
                            + invT_m_8*(rhsf(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhsf(i,j,l_nk-1))) &
                         - GVM%mc_css_H_8(i,j)   * (nl_t(i,j,l_nk ) - Ver_wmstar_8(G_nk)*nl_t(i,j,l_nk-1))

         enddo
      enddo
!!$omp end do
      if ( Schm_opentop_L ) then
!!$omp do
         do k=0, k0-2
            q_ext(:,:,k) = 0.0
         end do
!!$omp enddo
!!$omp do
         do j= j0, jn
            do i= i0, in
               w5= GVM%mc_alfat_8(i,j) * Sol_lhs(i,j,k0)
               w6= GVM%mc_cst_8  (i,j) * (rhsb(i,j)-nl_b(i,j))
               q_ext(i,j,k0-1) = w5 - w6
            end do
         end do
!!$omp enddo
      else
!--extrapolation at the lid   
!!$omp do
         do j= j0, jn
         do i= i0, in
            w5= (ver_z_8%m(0)-GVM%zmom_8(i,j,2))/(GVM%zmom_8(i,j,2)-GVM%zmom_8(i,j,1))
            q_ext(i,j,0)= Sol_lhs(i,j,2)*(1+w5) -w5*Sol_lhs(i,j,1)
         enddo
         enddo
!!$omp end do
      end if

!!$omp do collapse(2)
      do k=1, l_nk+1
         do j= j0, jn
            do i= i0, in
               qt0(i,j,k) = q_ext(i,j,k)
            end do
         end do
      end do
!!$omp enddo

      call HLT_split (1, G_nk+1, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( qt0(l_minx,l_miny,HLT_start),l_minx,l_maxx,&
                          l_miny,l_maxy,local_np,-1 )
      call HLT_split (0, G_nk+1, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( q_ext(l_minx,l_miny,HLT_start),l_minx,l_maxx,&
                          l_miny,l_maxy,local_np,-1 )

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, l_niu-pil_e

   !           Compute U
   !           ~~~~~~~~~
               ut0(i,j,k) = tau_m_8*(rhsu(i,j,k)-nl_u(i,j,k)  &
                          - (q_ext(i+1,j,k)-q_ext(i,j,k))*geomh_invDX_8(j) + &
                             GVM%mc_Jx_8(i,j,k) * ( &
                             wp_8(k)*half*( (q_ext(i+1,j,k+1)-q_ext(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )  &
                                          + (q_ext(i  ,j,k+1)-q_ext(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) )&
                            +wm_8(k)*half*( (q_ext(i+1,j,k  )-q_ext(i+1,j,k-1))*GVM%mc_iJz_8(i+1,j,k-1)&
                                          + (q_ext(i  ,j,k  )-q_ext(i  ,j,k-1))*GVM%mc_iJz_8(i  ,j,k-1) ) ))
            end do
         end do
      end do
!!$omp enddo nowait
      
!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0, l_njv-pil_n
!DIR$ SIMD
            do i= i0, in

   !           Compute V
   !           ~~~~~~~~~
               vt0(i,j,k) = tau_m_8*(rhsv(i,j,k)-nl_v(i,j,k) &
                          - (q_ext(i,j+1,k)-q_ext(i,j,k))*geomh_invDYMv_8(j) + &
                            GVM%mc_Jy_8(i,j,k) * ( &
                             wp_8(k)*half*( (q_ext(i,j+1,k+1)-q_ext(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )  &
                                          + (q_ext(i,j  ,k+1)-q_ext(i,j  ,k ))*GVM%mc_iJz_8(i,j  ,k ) )&
                            +wm_8(k)*half*( (q_ext(i,j+1,k  )-q_ext(i,j+1,k-1))*GVM%mc_iJz_8(i,j+1,k-1)&
                                          + (q_ext(i,j  ,k  )-q_ext(i,j  ,k-1))*GVM%mc_iJz_8(i,j  ,k-1) ) ))
            end do
         end do
      end do
!!$omp enddo nowait

!!$omp do collapse(2)
      do k=k0t, l_nk
         do j= j0, jn
!DIR$ SIMD
            do i= i0, in

   !           Compute w
   !           ~~~~~~~~~
               wt0(i,j,k) = tau_m_8 * ( (rhst(i,j,k) - nl_t(i,j,k)) &
                          - gama_bdf_8*(GVM%mc_iJz_8(i,j,k)*(q_ext(i,j,k+1)-q_ext(i,j,k)) &
                          - half*mu_8*(q_ext(i,j,k+1)+q_ext(i,j,k))))

   !           Compute zdot
   !           ~~~~~~~~~~~~
               zdt0(i,j,k) = (rhsf(i,j,k) + wt0(i,j,k))


   !           Compute T
   !           ~~~~~~~~~
               Buoy = (q_ext(i,j,k+1)-q_ext(i,j,k))*GVM%mc_iJz_8(i,j,k)  &
                    + wt0(i,j,k)*invT_nh_8 - (rhsw(i,j,k) - nl_w(i,j,k)) 
               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy / grav_8 )


            end do
         end do
      end do
!!$omp enddo
!
!     ---------------------------------------------------------------
!
      return
      end subroutine bac
