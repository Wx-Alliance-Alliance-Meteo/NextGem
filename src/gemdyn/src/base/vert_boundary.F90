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

!**s/r vert_boundary - 3D_elliptic lower boundary RHS computation for GEM_H

      subroutine vert_boundary ( i0, j0, in, jn )
      use geomh
      use gem_options
      use lam_options
      use HORgrid_options
      use dyn_fisl_options
      use glb_ld
      use mem_tstp
      use sol_mem
      use gmm_vt0
      use cstv
      use ver
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: i0,j0,in,jn

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i,j,k0,dim
      real(kind=REAL64) :: add_v8,bdd_v8,cdd_v8,w1
      real(kind=REAL64), parameter :: half=0.5d0
      real             , dimension(:,:,:), pointer :: wkf
      real(kind=REAL64), dimension(:,:,:), pointer :: wka, wkb
!
!     ---------------------------------------------------------------
!
      k0=1+Lam_gbpil_T
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      wkf(l_minx:l_maxx,l_miny:l_maxy,1:2) => WS1(1:)
      wka(l_minx:l_maxx,l_miny:l_maxy,1:2) => WS1_8(1:)
      wkb(l_minx:l_maxx,l_miny:l_maxy,1:2) => WS1_8(2*dim+1:)
!!$omp do
      do j=j0,jn
         wkf(:,j,:)= 0.0
         do i=i0,in
            wkf(i,j,1)= GVM%mc_css_H_8(i,j) * (rhst(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhst(i,j,l_nk-1)&
                     + Cstv_invT_m_8*(rhsf(i,j,l_nk)-Ver_wmstar_8(G_nk)*rhsf(i,j,l_nk-1)))   &
                     - GVM%mc_css_H_8(i,j) * (nl_t(i,j,l_nk )-Ver_wmstar_8(G_nk)*nl_t(i,j,l_nk-1))
            wkf(i,j,2)= -GVM%mc_cst_8(i,j) * (rhsb(i,j)-nl_b(i,j)) 
         end do
      end do
!!$omp enddo
!!$omp do
      do i= l_minx, l_maxx
         wkf(i,l_miny:j0-1,:) = 0.0
         wkf(i,jn+1:l_maxy,:) = 0.0
      end do
!!$omp enddo

      if ((Grd_yinyang_L)) then
         call yyg_xchng_hlt (wkf, l_minx,l_maxx,l_miny,l_maxy, &
                         l_ni,l_nj, 2, .false., 'CUBIC', .true.)
      else
         call HLT_split (1, 2, HLT_np, HLT_start, HLT_end)
         call gem_xch_halo ( wkf(l_minx,l_miny,HLT_start), &
                              l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      endif
                          
!!$omp do
      do j= pil_s, l_nj-pil_n
         do i= pil_w, l_ni-pil_e
            wka(i,j,1) = - GVM%mc_Jx_8(i,j,l_nk) * &
                        Ver_wp_8%m(l_nk)*half*( wkf(i+1,j,1)*GVM%mc_iJz_8(i+1,j,l_nk) &
                                            + wkf(i  ,j,1)*GVM%mc_iJz_8(i  ,j,l_nk) )
            wkb(i,j,1) = - GVM%mc_Jy_8(i,j,l_nk) * &
                        Ver_wp_8%m(l_nk)*half*( wkf(i,j+1,1)*GVM%mc_iJz_8(i,j+1,l_nk) &
                                            + wkf(i,j  ,1)*GVM%mc_iJz_8(i,j  ,l_nk) )
         end do
      end do
!!$omp enddo
      if (Schm_opentop_L) then
!!$omp do
      do j= pil_s, l_nj-pil_n
         do i= pil_w, l_ni-pil_e
            wka(i,j,2) = - GVM%mc_Jx_8(i,j,k0) * &
                        Ver_wm_8%m(k0)*half*( (-wkf(i+1,j,2))*GVM%mc_iJz_8(i+1,j,k0-1) &
                                            + (-wkf(i  ,j,2))*GVM%mc_iJz_8(i  ,j,k0-1) )
            wkb(i,j,2) = - GVM%mc_Jy_8(i,j,k0) * &
                        Ver_wm_8%m(k0)*half*( (-wkf(i,j+1,2))*GVM%mc_iJz_8(i,j+1,k0-1) &
                                            + (-wkf(i,j  ,2))*GVM%mc_iJz_8(i,j  ,k0-1) )
         end do
      end do
!!$omp enddo
      endif

      if (.not.Grd_yinyang_L) then
!!$omp do
      do i= l_minx,l_maxx
         if (l_south ) then
            wka(i,pil_s,:) = 0.
            wkb(i,pil_s,:) = 0.
         endif
         if (l_north ) then
            wkb(i,l_nj-pil_n,:) = 0.
         endif
      end do
!!$omp end do
!!$omp do
      do j= l_miny,l_maxy
         if (l_west ) then
            wka(pil_w,j,:) = 0.
            wkb(pil_w,j,:) = 0.
         endif
         if (l_east ) then
            wka(l_ni-pil_e,j,:) = 0.
         endif
      end do
!!$omp end do
      endif
      
!!$omp do
      do j=j0,jn
         do i=i0,in
            w1 = gama_8*(wkf(i,j,1)*GVM%mc_iJz_8(i,j,l_nk) - mu_8*half*wkf(i,j,1))
            add_v8 = (wka (i,j,1)-wka (i-1,j,1))*geomh_invDXM_8(j) &
                     + half * (GVM%mc_Ix_8(i,j,l_nk)*(wka(i,j,1)+wka(i-1,j,1)))
            bdd_v8 = (wkb (i,j,1)*geomh_cyM_8(j)-wkb (i,j-1,1)*&
                     geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                     + half * (GVM%mc_Iy_8(i,j,l_nk)*(wkb(i,j,1)+wkb(i,j-1,1)))
            cdd_v8 = w1*Ver_idz_8%m(l_nk) +(GVM%mc_Iz_8(i,j,l_nk)-epsi_8)*(Ver_wp_8%m(l_nk)*w1)

            Sol_rhs(i,j,l_nk)=Sol_rhs(i,j,l_nk)-Cstv_hco0_8*(add_v8+bdd_v8+cdd_v8)

         end do
      end do
!!$omp enddo
      if (Schm_opentop_L) then
!!$omp do
      do j=j0,jn
         do i=i0,in
            w1 = gama_8*(-wkf(i,j,2)*GVM%mc_iJz_8(i,j,k0-1) - mu_8*half*wkf(i,j,2))
            add_v8 = (wka (i,j,2)-wka (i-1,j,2))*geomh_invDXM_8(j) &
                     + half * (GVM%mc_Ix_8(i,j,k0)*(wka(i,j,2)+wka(i-1,j,2)))
            bdd_v8 = (wkb (i,j,2)*geomh_cyM_8(j)-wkb (i,j-1,2)*&
                     geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                     + half * (GVM%mc_Iy_8(i,j,k0)*(wkb(i,j,2)+wkb(i,j-1,2)))
            cdd_v8=(-w1)*Ver_idz_8%m(k0) +(GVM%mc_Iz_8(i,j,k0)-epsi_8)*(Ver_wm_8%m(k0)*w1) 
            Sol_rhs(i,j,k0)=Sol_rhs(i,j,k0)-Cstv_hco0_8*(add_v8+bdd_v8+cdd_v8)
         end do
      end do
!!$omp enddo
      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine vert_boundary

