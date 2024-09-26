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
!**s/r elliptic_rhs - Compute right hand side of the elliptic problem - 5th order

      subroutine elliptic_rhs5 ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use geomh
      use coriolis
      use adz_mem
      use gmm_vt0
      use mem_tstp
      use sol_mem
      use metric
      use dcst
      use tdpack
      use ver
      use vgh
      use glb_pil
      use stat_mpi
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      integer :: i00,inn,j00,jnn,dim,ub
      real, dimension(:,:,:), pointer :: tots, logT, logQ
      real(kind=REAL64), dimension(:,:,:), pointer :: t2u, v2u, t2v, u2v
      real(kind=REAL64), dimension(:,:,:), pointer :: dqz2u, dqz2v, dqz2w
      real(kind=REAL64) :: Rqq,tau_8,invT_8,a,b,c,barz,barzp
      real(kind=REAL64) :: w0,w1,w2,w3,w4, dudx,dvdy,ubx,vby
      real(kind=REAL64) :: dqdx, dqdy, ttbz, zzbz, dzrtt
      real(kind=REAL64) :: Ntttdz, Ntttbz, Nzzzdz, Nzzzbz
      real(kind=REAL128) :: dzrzz
      real(kind=REAL64), dimension(1:6) :: Nttt, Nwww
      real(kind=REAL64), parameter :: one=1.d0
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
      dim= l_ni*l_nj
      tots (1:l_ni,1:l_nj,1:l_nk) => WS1(ub+1:) ; ub=ub+dim*l_nk
      logT (1:l_ni,1:l_nj,1:l_nk) => WS1(ub+1:) ; ub=ub+dim*l_nk
      
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      ub=0
      t2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      v2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      t2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      u2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk

      !--- for vertical derivative of q to points u,v,w in rhs---
      dqz2u (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2v (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2w (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk

      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               tots(i,j,k)= tt0(i,j,k)/Cstv_Tstr_8
               logT(i,j,k)= log(tots(i,j,k))
            end do
         end do
      end do
      
      call prerhs5 (t2u, v2u, t2v, u2v, dqz2u, dqz2v, dqz2w,&
                    l_minx,l_maxx,l_miny,l_maxy,G_nk)

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
            dqdx = Hderiv(qt0(i-2,j,k), qt0(i-1,j,k), qt0(i,j,k), &
                          qt0(i+1,j,k), qt0(i+2,j,k), qt0(i+3,j,k), geomh_invDX_8(j))

            Nuu(i,j,k) =  t2u(i,j,k)*dqdx &
                        - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v2u(i,j,k) &
                        - t2u(i,j,k)*(dqdx - dqz2u(i,j,k)) 
         end do
         end do         
        
         do j= j00, jnn
         do i= ds_i0, ds_in
            dqdy = Hderiv(qt0(i,j-2,k),qt0(i,j-1,k),qt0(i  ,j,k),&
                          qt0(i,j+1,k),qt0(i,j+2,k),qt0(i,j+3,k),geomh_invDY_8)

            Nvv(i,j,k)= t2v(i,j,k)*dqdy &
                        + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u2v(i,j,k)) * u2v(i,j,k) &
                        - t2v(i,j,k)*(dqdy - dqz2v(i,j,k))
         end do
         end do

         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            w0= tots(i,j,k)-one
            w1= a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k)
            w2= a*rhsw_mid(i,j,k) - b*rhsw_dep(i,j,k)
            w3= invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
            !w4= w0*(dqz2w(i,j,k) - grav_8*(one-one/tots(i,j,k)))
            w4= w0*(GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) - grav_8*(one-one/tots(i,j,k)))

            Rtt(i,j,k)= w1 - w3
            Rww(i,j,k)= w2 - w4

            Rtt(i,j,k)= gama_bdf_8 * ( c*Rtt(i,j,k) + Rww(i,j,k) )
            Rzz(i,j,k)= a*rhsf_mid(i,j,k) - b*rhsf_dep(i,j,k ) - invT_8*(GVM%ztht_8(i,j,k)-Ver_z_8%t(k))
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

      call fill_Vhalo (Rtt,l_minx,l_maxx,l_miny,l_maxy,1.d0)
      call fill_Vhalo (Rzz,l_minx,l_maxx,l_miny,l_maxy,1.d0)
      
      do k=1, l_nk-1
         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            Rqq = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k )
            dudx = Hderiv8(Ruu(i-3,j,k), Ruu(i-2,j,k), Ruu(i-1,j,k), &
                          Ruu(i  ,j,k), Ruu(i+1,j,k), Ruu(i+2,j,k), geomh_invDXM_8(j))
            dvdy = Hderiv8(Rvv(i,j-3,k)*geomh_cyM_8(j-3), &
                           Rvv(i,j-2,k)*geomh_cyM_8(j-2), &
                           Rvv(i,j-1,k)*geomh_cyM_8(j-1), &
                           Rvv(i,j  ,k)*geomh_cyM_8(j  ), &
                           Rvv(i,j+1,k)*geomh_cyM_8(j+1), &
                           Rvv(i,j+2,k)*geomh_cyM_8(j+2), &
                           geomh_invDYM_8(j))
            ubx = Hstag8(Ruu(i-3,j,k), Ruu(i-2,j,k), Ruu(i-1,j,k), &
                         Ruu(i  ,j,k), Ruu(i+1,j,k), Ruu(i+2,j,k))
            vby = Hstag8(Rvv(i,j-3,k), Rvv(i,j-2,k), Rvv(i,j-1,k), &
                         Rvv(i,j  ,k), Rvv(i,j+1,k), Rvv(i,j+2,k))
            !compute vertical staggering of Rzz: thermo lvl -> mom lvl center point
            zzbz = Rzz(i,j,k-3) * VS5t2m(1,k) &
                 + Rzz(i,j,k-2) * VS5t2m(2,k) & 
                 + Rzz(i,j,k-1) * VS5t2m(3,k) & 
                 + Rzz(i,j,k  ) * VS5t2m(4,k) & 
                 + Rzz(i,j,k+1) * VS5t2m(5,k) & 
                 + Rzz(i,j,k+2) * VS5t2m(6,k)  

            !compute vertical staggering of Rtt: thermo lvl -> mom lvl center point
            ttbz = Rtt(i,j,k-3) * VS5t2m(1,k) & 
                 + Rtt(i,j,k-2) * VS5t2m(2,k) & 
                 + Rtt(i,j,k-1) * VS5t2m(3,k) & 
                 + Rtt(i,j,k  ) * VS5t2m(4,k) & 
                 + Rtt(i,j,k+1) * VS5t2m(5,k) & 
                 + Rtt(i,j,k+2) * VS5t2m(6,k)  

            !compute vertical derivative of Rzz: thermo lvl -> mom lvl center point
            dzrzz = Rzz(i,j,k-3) * VD5t2m(1,k) &
                  + Rzz(i,j,k-2) * VD5t2m(2,k) & 
                  + Rzz(i,j,k-1) * VD5t2m(3,k) & 
                  + Rzz(i,j,k  ) * VD5t2m(4,k) & 
                  + Rzz(i,j,k+1) * VD5t2m(5,k) & 
                  + Rzz(i,j,k+2) * VD5t2m(6,k)  

            !compute vertical derivative of Rtt: thermo lvl -> mom lvl cntr point
            dzrtt = Rtt(i,j,k-3) * VD5t2m(1,k) & 
                  + Rtt(i,j,k-2) * VD5t2m(2,k) & 
                  + Rtt(i,j,k-1) * VD5t2m(3,k) & 
                  + Rtt(i,j,k  ) * VD5t2m(4,k) & 
                  + Rtt(i,j,k+1) * VD5t2m(5,k) & 
                  + Rtt(i,j,k+2) * VD5t2m(6,k)

            ! exact form of eqn 58 in SG notes
            Sol_rhs(i,j,k) = -invT_8*Rqq + dudx + dvdy + invT_8*dzrzz   &
                            + ubx*M_logJzu(i,j,k) + vby*M_logJzv(i,j,k) &
                            + invT_8*M_logJzq(i,j,k)* zzbz + dzrtt      &
                            + ttbz*(M_logJzq(i,j,k) - epsi_8)
         end do
         end do
!         call statf_dm (Sol_rhs(1:,1:,k:k), 'ERHS', k, 'ELLI', 1,ubound(Sol_rhs,1),1,ubound(Sol_rhs,2),1,1,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,1,8)
      end do
!
!     ---------------------------------------------------------------
!
      return
      include 'H5th_ope.inc'
      end subroutine elliptic_rhs5
