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
!**s/r elliptic_rhs_schar - Compute right hand side of the elliptic problem

      subroutine SW_elliptic_rhs ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use geomh
      use coriolis
      use adz_mem
      use gmm_vt0
      use gmm_geof
      use mem_tstp
      use sol_mem
      use metric
      use dcst
      use cstv
      use tdpack
      use ver
      use glb_pil
      use stat_mpi
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      integer :: km,i00,inn,j00,jnn,dim,ub
      real, dimension(:,:,:), pointer :: wkf
      real, dimension(:,:,:), pointer :: tots, logT, logQ, Rt
      real(kind=REAL64) :: Rqq,Nqq,tau_8,invT_8,a,b,c,barz,barzp,div
      real(kind=REAL64) :: w0,w1,w2,w3,w4,w5,w6,w7, dudx,dvdy,ubx,vby
      real(kind=REAL64) :: Nwww,Nttt,t_interp,u_interp,v_interp
      real(kind=REAL64) :: wka(l_minx:l_maxx,l_miny:l_maxy),&
                           wkb(l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL128) :: d1
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      Sol_rhs=0.
!     if (Schm_POSO == 3) then
!        call SW_elliptic_rhs3rd( F_dt_8, k0, k0t )
!        return
!     endif

      i00= ds_i0-1 ; inn= ds_in
      j00= ds_j0-1 ; jnn= ds_jn
      if (.not.Grd_yinyang_L) then
         i00= ds_i0-1+min(pil_w,1)
         inn= ds_in  -min(pil_e,1)
         j00= ds_j0-1+min(pil_s,1)
         jnn= ds_jn  -min(pil_n,1)
      endif
      
      tau_8  = (2.d0 * F_dt_8 ) / 3.d0
      invT_8 = 1.d0/tau_8
      a      = 4.d0*invT_8/3.d0
      b      =      invT_8/3.d0
      c      = grav_8 * tau_8

      do k=1, l_nk
         do j= Adz_j0 , Adz_jn
         do i= Adz_i0u, Adz_inu
            Ruu(i,j,k) = a*rhsu_mid(i,j,k) - b*rhsu_dep(i,j,k) 
         end do
         end do
         do j= Adz_j0v, Adz_jnv
         do i= Adz_i0 , Adz_in
            Rvv(i,j,k) = a*rhsv_mid(i,j,k) - b*rhsv_dep(i,j,k) 
         end do
         end do
         do j= ds_j0, ds_jn
         do i= i00, inn
            v_interp = 0.25d0*(vt0(i  ,j,k)+vt0(i  ,j-1,k)+&
                               vt0(i+1,j,k)+vt0(i+1,j-1,k))
            Nuu(i,j,k)=  - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v_interp 
         end do
         end do         
         do j= j00, jnn
         do i= ds_i0, ds_in
            u_interp = 0.25d0*(ut0(i,j,k)+ut0(i-1,j,k)+ut0(i,j+1,k)+ut0(i-1,j+1,k))
            Nvv(i,j,k) = ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp 
         end do
         end do
      end do
      
      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( Ruu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( Rvv(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      do k=1, l_nk
         Ruu(:,:,k)= Ruu(:,:,k) - Nuu(:,:,k)
         Rvv(:,:,k)= Rvv(:,:,k) - Nvv(:,:,k)
      end do

      do k=1, l_nk
         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            Rqq = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k ) + (1.d0-Cstv_swln_8)*invT_8*fis0(i,j)

            div = (ut0 (i,j,k)- ut0 (i-1,j,k))*geomh_invDXM_8(j)     &
                + (vt0 (i,j,k)*geomh_cyM_8(j)-vt0 (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) 

            Nqq = (1.d0-Cstv_swln_8)*(qt0(i,j,k)-fis0(i,j)) * div &                                                                                         
                            -Cstv_swln_8*invT_8*(Cstv_h0inv_8*qt0(i,j,k)-log(Cstv_h0inv_8*(qt0(i,j,k)-fis0(i,j))+1.d0))

            Rqq = Rqq - Nqq

            dudx= (Ruu(i,j,k)-Ruu(i-1,j,k))*geomh_invDXM_8(j)
            dvdy= (Rvv(i,j  ,k)*geomh_cyM_8(j  ) - &
                   Rvv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j)

            Sol_rhs(i,j,k) = (dudx + dvdy)/grav_8 - invT_8*((1.d0-Cstv_swln_8)*Cstv_h0inv_8+Cstv_swln_8)/grav_8*Rqq
         end do
         end do
      end do
      
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_elliptic_rhs
