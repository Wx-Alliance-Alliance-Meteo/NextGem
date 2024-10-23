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

      subroutine SW_elliptic_rhs3rd ( F_dt_8, k0, k0t )
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
      use yyg_param
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      integer :: km,i00,inn,j00,jnn,dim,ub
      real, dimension(:,:,:), pointer :: wkf
      real, dimension(:,:,:), pointer :: tots, logT, logQ, Rt
      real(kind=REAL64), dimension(:,:,:), pointer :: v2u, u2v
      real(kind=REAL64) :: Rqq,Nqq,tau_8,invT_8,a,b,c,barz,barzp,div
      real(kind=REAL64) :: dqdx, dqdy
      real(kind=REAL64) :: w0,w1,w2,w3,w4,w5,w6,w7, dudx,dvdy,ubx,vby
      real(kind=REAL64) :: Nwww,Nttt,t_interp,u_interp,v_interp
      real(kind=REAL64) :: wka(l_minx:l_maxx,l_miny:l_maxy),&
                           wkb(l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL128) :: d1
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      i00= ds_i0-3 ; inn= ds_in+2
      j00= ds_j0-3 ; jnn= ds_jn+2
      if (.not.Grd_yinyang_L) then
         i00= ds_i0-(1-min(pil_w,1))*3
         inn= ds_in+2-3*min(pil_e,1)
         j00= ds_j0-(1-min(pil_s,1))*3
         jnn= ds_jn+2-3*min(pil_n,1)
      endif
      
      tau_8  = (2.d0 * F_dt_8 ) / 3.d0
      invT_8 = 1.d0/tau_8
      a      = 4.d0*invT_8/3.d0
      b      =      invT_8/3.d0
      c      = grav_8 * tau_8

      ub=0
      v2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      u2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk

      call SW_prerhs3rd (v2u,u2v,l_minx,l_maxx,l_miny,l_maxy,G_nk)

      do k=1, l_nk
        !do j= 1, l_nj
        !do i= 1, l_ni
         do j= Adz_j0 , Adz_jn
         do i= Adz_i0u, Adz_inu
            Ruu(i,j,k) = a*rhsu_mid(i,j,k) - b*rhsu_dep(i,j,k) 
         end do
         end do
        !do j= 1, l_nj
        !do i= 1, l_ni
         do j= Adz_j0v, Adz_jnv
         do i= Adz_i0 , Adz_in
            Rvv(i,j,k) = a*rhsv_mid(i,j,k) - b*rhsv_dep(i,j,k) 
         end do
         end do
        !do j= 1, l_nj
        !do i= 1, l_ni
         do j= ds_j0, ds_jn
         do i= i00, inn
            Nuu(i,j,k)=  - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j)*ut0(i,j,k) ) * v2u(i,j,k) 
         end do
         end do         
        !do j= 1, l_nj
        !do i= 1, l_ni
         do j= j00, jnn
         do i= ds_i0, ds_in
            Nvv(i,j,k) = ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j)*u2v(i,j,k) ) * u2v(i,j,k)
         end do
         end do
      end do
      
      do k=1, l_nk
         Ruu(:,:,k)= Ruu(:,:,k) - Nuu(:,:,k)
         Rvv(:,:,k)= Rvv(:,:,k) - Nvv(:,:,k)
      end do
!              call yyg_xchng_8 (Ruu(l_minx,l_miny,1), YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
!                          l_ni,l_nj, l_nk, .false., 'CUBIC', .true.)
!              call yyg_xchng_8 (Rvv(l_minx,l_miny,1), YYG_HALO_q2q, l_minx,l_maxx,l_miny,l_maxy, &
!                          l_ni,l_nj, l_nk, .false., 'CUBIC', .true.)

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( Ruu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( Rvv(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      do k=1, l_nk
         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            Rqq = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k ) + (1.d0-Cstv_swln_8)*invT_8*fis0(i,j)

             dudx = Hderiv8(ut0(i-2,j,k)*1.d0, ut0(i-1,j,k)*1.d0, &
                            ut0(i  ,j,k)*1.d0, ut0(i+1,j,k)*1.d0, geomh_invDXM_8(j))
             dvdy = Hderiv8(vt0(i,j-2,k)*geomh_cyM_8(j-2), &
                            vt0(i,j-1,k)*geomh_cyM_8(j-1), &
                            vt0(i,j  ,k)*geomh_cyM_8(j  ), &
                            vt0(i,j+1,k)*geomh_cyM_8(j+1), &
                            geomh_invDYM_8(j))
            !dudx = Hderiv8(ut0(i-1 ,j,k)*1.d0, ut0(i  ,j,k)*1.d0, &
            !               ut0(i+1 ,j,k)*1.d0, ut0(i+2,j,k)*1.d0, geomh_invDXM_8(j))
            !dvdy = Hderiv8(vt0(i,j-1,k)*geomh_cyM_8(j-1), &
            !               vt0(i,j,k)*geomh_cyM_8(j), &
            !               vt0(i,j+1,k)*geomh_cyM_8(j+1  ), &
            !              vt0(i,j+2,k)*geomh_cyM_8(j+2), &
            !              geomh_invDYM_8(j))
           !dudx= (ut0 (i,j,k)- ut0 (i-1,j,k))*geomh_invDXM_8(j)    
           !dvdy=  (vt0 (i,j,k)*geomh_cyM_8(j)-vt0 (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) 

            Nqq = (1.d0-Cstv_swln_8)*(qt0(i,j,k)-fis0(i,j)) * ( dudx + dvdy ) &                                                                                         
                            -Cstv_swln_8*invT_8*(Cstv_h0inv_8*qt0(i,j,k)-log(Cstv_h0inv_8*(qt0(i,j,k)-fis0(i,j))+1.d0))

            Rqq = Rqq - Nqq

            dudx = Hderiv8(Ruu(i-2,j,k), Ruu(i-1,j,k), &
                           Ruu(i  ,j,k), Ruu(i+1,j,k), geomh_invDXM_8(j))
            dvdy = Hderiv8(Rvv(i,j-2,k)*geomh_cyM_8(j-2), &
                           Rvv(i,j-1,k)*geomh_cyM_8(j-1), &
                           Rvv(i,j  ,k)*geomh_cyM_8(j  ), &
                           Rvv(i,j+1,k)*geomh_cyM_8(j+1), &
                           geomh_invDYM_8(j))
!                  dudx=0.
!                  dvdy=0.
           !dudx = Hderiv8(Ruu(i-1,j,k), Ruu(i,j,k), &
           !               Ruu(i+1  ,j,k), Ruu(i+2,j,k), geomh_invDXM_8(j))
           !dvdy = Hderiv8(Rvv(i,j-1,k)*geomh_cyM_8(j-1), &
           !               Rvv(i,j  ,k)*geomh_cyM_8(j  ), &
           !               Rvv(i,j+1,k)*geomh_cyM_8(j+1), &
           !               Rvv(i,j+2,k)*geomh_cyM_8(j+2), &
           !               geomh_invDYM_8(j+1))

            Sol_rhs(i,j,k) = (dudx + dvdy)/grav_8 - invT_8*((1.d0-Cstv_swln_8)*Cstv_h0inv_8+Cstv_swln_8)/grav_8*Rqq
         end do
         end do
      end do
      
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine SW_elliptic_rhs3rd
