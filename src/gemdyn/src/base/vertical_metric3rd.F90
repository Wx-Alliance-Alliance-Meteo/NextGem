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

!**s/r vertical_metric - Compute vertical metric coefficients

      subroutine vertical_metric3rd ()
      use, intrinsic :: iso_fortran_env
      use gem_options
      use dyn_fisl_options
      use HORgrid_options
      use geomh
      use dcst
      use gmm_vt1
      use tdpack
      use glb_ld
      use metric
      use gmm_geof
      use ver
      implicit none

      integer :: i,j,k,kstart
      integer :: HLT_start, HLT_end, local_np, k0,kn
      real(kind=REAL64) :: Jzu, Jzv, Jzq, Jx, Jy
      real(kind=REAL64), dimension(:,:,:), pointer :: M_Jzu , M_Jzv , M_Jzq
      real, parameter :: one=1.d0, half=.5d0
!
!     ---------------------------------------------------------------
!
      kstart= ubound(VM3%zmom,3)
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            VM3%pstar(i,j,kstart)=log(1.d5)-grav_8*VM3%zmom(i,j,kstart)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do
      do j=1-G_haloy,l_nj+G_haloy
         do k=kstart-1,lbound(VM3%zmom,3),-1
            do i=1-G_halox,l_ni+G_halox
               VM3%pstar(i,j,k)=VM3%pstar(i,j,k+1)+grav_8*(VM3%zmom(i,j,k+1)-VM3%zmom(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
      end do

      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            GVM%lg_pstar_8(i,j,G_nk+1)=log(1.d5)-grav_8*GVM%zmom_8(i,j,G_nk+1)/(rgasd_8*Cstv_Tstr_8)
         end do
      end do

      do j=1-G_haloy,l_nj+G_haloy
         do k=G_nk,1,-1
            do i=1-G_halox,l_ni+G_halox
               GVM%lg_pstar_8(i,j,k)=GVM%lg_pstar_8(i,j,k+1)+grav_8*(GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))/(rgasd_8*Cstv_Tstr_8)
            end do
         end do
         do i=1-G_halox,l_ni+G_halox
            GVM%ztht_8(i,j,G_nk)= GVM%zmom_8(i,j,G_nk+1) !temporary for mc_Ix_8 and mc_Iy_8 below
         end do
      end do

      k0 = lbound(VM3%zmom,3)
      kn = ubound(VM3%zmom,3)-1
      allocate ( M_Jzu(l_minx:l_maxx,l_miny:l_maxy,k0:kn),&
                 M_Jzv(l_minx:l_maxx,l_miny:l_maxy,k0:kn),&
                 M_Jzq(l_ni,l_nj,k0:kn) )

      do k=1,ubound(VM3%zmom,3)-2
         do j=1-G_haloy+1,l_nj+G_haloy-2
            do i=1-G_halox+1,l_ni+G_halox-2
               Jzu = VM3%ztht_u(i,j,k-2) * VD3t2m(1,k)&
                    +VM3%ztht_u(i,j,k-1) * VD3t2m(2,k)&
                    +VM3%ztht_u(i,j,k  ) * VD3t2m(3,k)&
                    +VM3%ztht_u(i,j,k+1) * VD3t2m(4,k)
               Jzv = VM3%ztht_v(i,j,k-2) * VD3t2m(1,k)&
                    +VM3%ztht_v(i,j,k-1) * VD3t2m(2,k)&
                    +VM3%ztht_v(i,j,k  ) * VD3t2m(3,k)&
                    +VM3%ztht_v(i,j,k+1) * VD3t2m(4,k)
               Jx= Hderiv8(VM3%zmom(i-1,j,k),VM3%zmom(i  ,j,k),&
                           VM3%zmom(i+1,j,k),VM3%zmom(i+2,j,k),geomh_invDX_8(j))
               Jy= Hderiv8(VM3%zmom(i,j-1,k),VM3%zmom(i  ,j,k),&
                           VM3%zmom(i,j+1,k),VM3%zmom(i,j+2,k),geomh_invDY_8)
               M_Jxozu(i,j,k)= Jx  / Jzu ! First  term in matvec portion of eqn 58
               M_Jyozv(i,j,k)= Jy  / Jzv ! Second term in matvec portion of eqn 58
               M_Jzu(i,j,k) = log(Jzu)
               M_Jzv(i,j,k) = log(Jzv)
            end do
         end do
      end do
      do k=lbound(VM3%zmom,3)+2,0
      do j=1-G_haloy+1,l_nj+G_haloy-2
         do i=1-G_halox+1,l_ni+G_halox-2
            Jzu = VM3%ztht_u(i,j,k-3) * VD3t2m(1,k)&
                 +VM3%ztht_u(i,j,k-2) * VD3t2m(2,k)&
                 +VM3%ztht_u(i,j,k-1) * VD3t2m(3,k)&
                 +VM3%ztht_u(i,j,k  ) * VD3t2m(4,k)
            Jzv = VM3%ztht_v(i,j,k-3) * VD3t2m(1,k)&
                 +VM3%ztht_v(i,j,k-2) * VD3t2m(2,k)&
                 +VM3%ztht_v(i,j,k-1) * VD3t2m(3,k)&
                 +VM3%ztht_v(i,j,k  ) * VD3t2m(4,k)
            Jx= Hderiv8(VM3%zmom(i-1,j,k),VM3%zmom(i  ,j,k),&
                        VM3%zmom(i+1,j,k),VM3%zmom(i+2,j,k),geomh_invDX_8(j))
            Jy= Hderiv8(VM3%zmom(i,j-1,k),VM3%zmom(i  ,j,k),&
                        VM3%zmom(i,j+1,k),VM3%zmom(i,j+2,k),geomh_invDY_8)
            M_Jxozu(i,j,k)= Jx  / Jzu ! First  term in matvec portion of eqn 58
            M_Jyozv(i,j,k)= Jy  / Jzv ! Second term in matvec portion of eqn 58
            M_Jzu(i,j,k) = log(Jzu)
            M_Jzv(i,j,k) = log(Jzv)
         end do
      end do
      end do

      do k=0,G_nk
         do j=1,l_nj
            do i=1,l_ni
               Jzq = VM3%zmom  (i,j,k-1) * VD3m2t(1,k)&
                    +VM3%zmom  (i,j,k  ) * VD3m2t(2,k)&
                    +VM3%zmom  (i,j,k+1) * VD3m2t(3,k)&
                    +VM3%zmom  (i,j,k+2) * VD3m2t(4,k)
               M_iJzq (i,j,k)= one / Jzq ! Third  term in matvec portion of eqn 58
               M_Jzq(i,j,k) = log(Jzq)
            end do
         end do
      end do
      do k=lbound(VM3%zmom,3),-1
         do j=1,l_nj
            do i=1,l_ni
               Jzq = VM3%zmom  (i,j,k  ) * VD3m2t(1,k)&
                    +VM3%zmom  (i,j,k+1) * VD3m2t(2,k)&
                    +VM3%zmom  (i,j,k+2) * VD3m2t(3,k)&
                    +VM3%zmom  (i,j,k+3) * VD3m2t(4,k)
               M_iJzq (i,j,k)= one / Jzq ! Third  term in matvec portion of eqn 58
               M_Jzq(i,j,k) = log(Jzq)
            end do
         end do
      end do
      do k=G_nk+1,ubound(VM3%zmom,3)-1
         do j=1,l_nj
            do i=1,l_ni
               Jzq = VM3%zmom  (i,j,k-2) * VD3m2t(1,k)&
                    +VM3%zmom  (i,j,k-1) * VD3m2t(2,k)&
                    +VM3%zmom  (i,j,k  ) * VD3m2t(3,k)&
                    +VM3%zmom  (i,j,k+1) * VD3m2t(4,k)
               M_iJzq (i,j,k)= one / Jzq ! Third  term in matvec portion of eqn 58
               M_Jzq(i,j,k) = log(Jzq)
            end do
         end do
      end do
         
      do k=1,G_nk
         do j=1,l_nj
            do i=1,l_ni
               ! 5th  term in elliptic RHS portion of eqn 58
               ! 5th  term in matvec portion of eqn 58
               M_logJzu(i,j,k)= Hderiv8(M_Jzu(i-1,j,k),M_Jzu(i  ,j,k),&
                                        M_Jzu(i+1,j,k),M_Jzu(i+2,j,k),geomh_invDX_8(j))
               ! 6th  term in elliptic RHS portion of eqn 58
               ! 6th  term in matvec portion of eqn 58
               M_logJzv(i,j,k)= Hderiv8(M_Jzv(i,j-1,k),M_Jzv(i  ,j,k),M_Jzv(i,j+1,k),M_Jzv(i,j+2,k),geomh_invDY_8)
               ! 7th and 9th  term in elliptic RHS portion of eqn 58
               ! 7th  term in matvec portion of eqn 58
               M_logJzq(i,j,k)= M_Jzq(i,j,k-2) * VD3t2m(1,k)&
                               +M_Jzq(i,j,k-1) * VD3t2m(2,k)&
                               +M_Jzq(i,j,k  ) * VD3t2m(3,k)&
                               +M_Jzq(i,j,k+1) * VD3t2m(4,k)
            end do
         end do
      end do
!!$omp do collapse(2)
      do k=1,G_nk
         do j=1-G_haloy+1,l_nj+G_haloy-2
            do i=1-G_halox+1,l_ni+G_halox-2
               GVM%mc_Jx_8 (i,j,k)= Hderiv8(GVM%zmom_8(i-1,j,k),GVM%zmom_8(i  ,j,k),&
                                            GVM%zmom_8(i+1,j,k),GVM%zmom_8(i+2,j,k),geomh_invDX_8(j))
               GVM%mc_Jy_8 (i,j,k)= Hderiv8(GVM%zmom_8(i,j-1,k),GVM%zmom_8(i  ,j,k),&
                                            GVM%zmom_8(i,j+1,k),GVM%zmom_8(i,j+2,k),geomh_invDY_8)
               GVM%mc_iJz_8(i,j,k)=one/(GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))
               GVM%mc_Ix_8(i,j,k)=log( (zthtu_8(i,j,k)-zthtu_8(i,j,k-1))/(zthtu_8(i-1,j,k)-zthtu_8(i-1,j,k-1)) )*geomh_invDX_8(j)
               GVM%mc_Iy_8(i,j,k)=log( (zthtv_8(i,j,k)-zthtv_8(i,j,k-1))/(zthtv_8(i,j-1,k)-zthtv_8(i,j-1,k-1)) )*geomh_invDY_8
               GVM%mc_Iz_8(i,j,k)=log( (GVM%zmom_8(i,j,k+1)-GVM%zmom_8(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
                                  /(GVM%zmom_8(i,j,k)-GVM%zmom_8(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
               GVM%mc_logJz_8(i,j,k)= 0.0
            end do
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-2
!DIR$ SIMD
         do i=1-G_halox+1,l_ni+G_halox-2
            GVM%mc_iJz_8(i,j,0)=one/(GVM%zmom_8(i,j,1)-ver_z_8%m(0))
            GVM%mc_css_H_8(i,j) = one/(gama_8*(GVM%mc_iJz_8(i,j,G_nk)-half*mu_8))
         end do
      end do
!!$omp enddo
!!$omp do
      do j=1-G_haloy,l_nj+G_haloy
!DIR$ SIMD
         do i=1-G_halox,l_ni+G_halox
            GVM%ztht_8(i,j,G_nk)= ver_z_8%t(G_nk)+(Ver_b_8%t(G_nk)*fis0(i,j)+Ver_c_8%t(G_nk)*orols(i,j))/grav_8
         end do
      end do
!!$omp enddo

!!$omp do
      do j=1-G_haloy+1,l_nj+G_haloy-2
         do i=1-G_halox+1,l_ni+G_halox-2
            GVM%mc_alfas_H_8(i,j) = ( GVM%mc_iJz_8(i,j,G_nk) + half*mu_8 + Ver_wmstar_8(G_nk)*(GVM%mc_iJz_8(i,j,G_nk-1) -half*mu_8) ) / (GVM%mc_iJz_8(i,j,G_nk)-half*mu_8)
            GVM%mc_betas_H_8(i,j) =                                             Ver_wmstar_8(G_nk)*(GVM%mc_iJz_8(i,j,G_nk-1) +half*mu_8)   / (GVM%mc_iJz_8(i,j,G_nk)-half*mu_8)
         enddo
      enddo
!     !$omp enddo
      deallocate (M_Jzu , M_Jzv , M_Jzq)
!
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine vertical_metric3rd
