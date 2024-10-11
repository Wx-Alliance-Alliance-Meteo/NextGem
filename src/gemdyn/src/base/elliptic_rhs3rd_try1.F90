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

      subroutine elliptic_rhs3rd_try1 ( F_dt_8, k0, k0t )
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
      use glb_pil
      use stat_mpi
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, n, HLT_np, HLT_start, HLT_end
      integer :: km,i00,inn,j00,jnn,dim,ub
      integer :: km1,km2,km3,kp1,kp2,kp3
      real, dimension(:,:,:), pointer :: wkf
      real, dimension(:,:,:), pointer :: tots, logT, logQ, Rt
      real(kind=REAL64), dimension(:,:,:), pointer :: t2u, v2u, t2v, u2v
      real(kind=REAL64), dimension(:,:,:), pointer :: dqz2u, dqz2v, dqz2w
      real(kind=REAL64) :: Rqq,tau_8,invT_8,a,b,c,barz,barzp
      real(kind=REAL64) :: dqdx, dqdy, ttbz, zzbz, dzrtt
      real(kind=REAL64) :: Ntttdz, Ntttbz, Nzzzdz, Nzzzbz
      real(kind=REAL64) :: t_interp,u_interp,v_interp
      real(kind=REAL64) :: wka(l_minx:l_maxx,l_miny:l_maxy),&
                           wkb(l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL128) :: d1
      real(kind=REAL64) :: w0,w1,w2,w3,w4,w5,w7, dudx,dvdy,ubx,vby
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
      real(kind=REAL64), dimension(-1:2) :: t2qv, t2qu, v2q, u2q, dq2u, dq2v
      real(kind=REAL64) :: duu, dvv 
      real(kind=REAL128) :: dzrzz,rtm(l_nk)
      real(kind=REAL64), dimension(1:6) :: Nttt, Nwww, w6
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) :: dqzv, dqzu
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
      Rt   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(ub+1:) ; ub=ub+dim*l_nk
      wkf  (l_minx:l_maxx,l_miny:l_maxy,1:1)    => WS1(ub+1:) ; ub=ub+dim
      ub=0

      t2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      v2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      t2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      u2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk

      !---vertical derivative of q to points u,v,w in rhs---
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


      !--- prepare RHS ---
      do k=1, l_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         kp3=min(k+3,G_nk)
         do j= 1, l_nj
           do i= 1, l_ni
             do n=-1,2
                t2qv(n) = tt0(i+n,j,km2)*CWt2m(1,k) + tt0(i+n,j,km1)*CWt2m(2,k) &
                        + tt0(i+n,j,k  )*CWt2m(3,k) + tt0(i+n,j,kp1)*CWt2m(4,k)

               v2q(n) = Hstag(vt0(i,j-2,k),vt0(i,j-1,k),&
                              vt0(i,j,k  ),vt0(i,j+1,k))
               t2qu(n) = tt0(i,j+n,km2)*CWt2m(1,k) + tt0(i,j+n,km1)*CWt2m(2,k) &
                       + tt0(i,j+n,k  )*CWt2m(3,k) + tt0(i,j+n,kp1)*CWt2m(4,k)

               u2q(n) = Hstag(ut0(i-2,j,k),ut0(i-1,j,k),&
                              ut0(i,j,k  ),ut0(i+1,j,k))
               dq2u(n) = qt0(i+n,j,km1) * CDm2t(1,k)&
                        +qt0(i+n,j,k  ) * CDm2t(2,k)&
                        +qt0(i+n,j,k+1) * CDm2t(3,k)&
                        +qt0(i+n,j,kp2) * CDm2t(4,k)

               dq2v(n) = qt0(i,j+n,km1) * CDm2t(1,k)&
                        +qt0(i,j+n,k  ) * CDm2t(2,k)&
                        +qt0(i,j+n,k+1) * CDm2t(3,k)&
                        +qt0(i,j+n,kp2) * CDm2t(4,k)
             end do
             t2u(i,j,k)= Hstag8(t2qv(-1),t2qv(0),t2qv(1),t2qv(2))/Cstv_Tstr_8-1.d0
             v2u(i,j,k)= Hstag8( v2q(-1), v2q(0), v2q(1), v2q(2))
             t2v(i,j,k)= Hstag8(t2qu(-1),t2qu(0),t2qu(1),t2qu(2))/Cstv_Tstr_8-1.d0
             u2v(i,j,k)= Hstag8( u2q(-1), u2q(0), u2q(1), u2q(2))
             dqz2w(i,j,k) = M_iJzq(i,j,k)*dq2u(0) !for Rw
             dqzu(i,j,k)= Hstag8(dq2u(-1),dq2u(0),dq2u(1),dq2u(2)) !values in thermo, u-grid
             dqzv(i,j,k)= Hstag8(dq2v(-1),dq2v(0),dq2v(1),dq2v(2)) !values in thermo, v-grid
         end do
         end do
      end do

      !---now interpolate back to momentum level---
      !Note: dqdzu and dqdzv are already staggered appropriatly from previous loop
      do k=1,l_nk
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)
         do j= 1, l_nj 
            do i= 1, l_ni
               duu = dqzu(i,j,km2) * CWt2m(1,k)&
                   + dqzu(i,j,km1) * CWt2m(2,k)&
                   + dqzu(i,j,k  ) * CWt2m(3,k)&
                   + dqzu(i,j,kp1) * CWt2m(4,k)

               dvv = dqzv(i,j,km2) * CWt2m(1,k)&
                   + dqzv(i,j,km1) * CWt2m(2,k)&
                   + dqzv(i,j,k  ) * CWt2m(3,k)&
                   + dqzv(i,j,kp1) * CWt2m(4,k)

               dqz2u(i,j,k)= M_Jxozu(i,j,k) * duu 
               dqz2v(i,j,k)= M_Jyozv(i,j,k) * dvv
            end do
         end do
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( t2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( v2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( t2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( u2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( dqz2u(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( dqz2v(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( dqz2w(l_minx,l_miny,HLT_start),&
                            l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!     
!     ---------------------------------------------------------------
!

      do k=1, l_nk
         km=max(k-1,1)
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

!        do j= ds_j0, ds_jn
!        do i= i00, inn
         do j= 1, l_nj
          do i= 1, l_ni
            dqdx = Hderiv(qt0(i-1,j,k), qt0(i,j,k), &
                          qt0(i+1,j,k), qt0(i+2,j,k), geomh_invDX_8(j))
            Nuu(i,j,k) = - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v2u(i,j,k) &
                        - t2u(i,j,k)*(dqdx - dqz2u(i,j,k)) 
          end do
         end do         
!        do j= j00, jnn
!        do i= ds_i0, ds_in
         do j= 1, l_nj
          do i= 1, l_ni
            dqdy = Hderiv(qt0(i,j-1,k),qt0(i  ,j,k),&
                          qt0(i,j+1,k),qt0(i,j+2,k),geomh_invDY_8)
            Nvv(i,j,k)= ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u2v(i,j,k)) * u2v(i,j,k) &
                        - t2v(i,j,k)*(dqdy - dqz2v(i,j,k))
          end do
         end do

         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
            w0= tots(i,j,k)-one
            w1= a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k)
            w2= a*rhsw_mid(i,j,k) - b*rhsw_dep(i,j,k)
            w3= invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
            w4= w0*(GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) - grav_8*(one-one/tots(i,j,k)))
            Rtt(i,j,k)= w1 - w3
            Rww(i,j,k)= w2 - w4
            Rtt(i,j,k)= gama_bdf_8 * ( c*Rtt(i,j,k) + Rww(i,j,k) )
            Rzz(i,j,k)= a*rhsf_mid(i,j,k) - b*rhsf_dep(i,j,k ) - invT_8*(GVM%ztht_8(i,j,k)-Ver_z_8%t(k))
            
            Rt(i,j,k) = gama_bdf_8 * ( c*w1 + w2 ) ! special one for wkf below
         end do
         end do
      end do
      
      do k=1, l_nk
         Ruu(:,:,k)= Ruu(:,:,k) - Nuu(:,:,k)
         Rvv(:,:,k)= Rvv(:,:,k) - Nvv(:,:,k)
      end do

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo_8 ( Ruu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
      call gem_xch_halo_8 ( Rvv(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      do k=1, l_nk
         km=max(k-1,1)
         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)

         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            Rqq = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k )
             dudx = Hderiv8(Ruu(i-2,j,k), Ruu(i-1,j,k), &
                           Ruu(i  ,j,k), Ruu(i+1,j,k), geomh_invDXM_8(j))
             dvdy = Hderiv8(Rvv(i,j-2,k)*geomh_cyM_8(j-2), &
                           Rvv(i,j-1,k)*geomh_cyM_8(j-1), &
                           Rvv(i,j  ,k)*geomh_cyM_8(j  ), &
                           Rvv(i,j+1,k)*geomh_cyM_8(j+1), &
                           geomh_invDYM_8(j))
             ubx = Hstag8(Ruu(i-2,j,k), Ruu(i-1,j,k), &
                         Ruu(i  ,j,k), Ruu(i+1,j,k))
             vby = Hstag8(Rvv(i,j-2,k), Rvv(i,j-1,k), &
                         Rvv(i,j  ,k), Rvv(i,j+1,k))
            !compute vertical staggering of Rzz: thermo lvl -> mom lvl center point
            zzbz = Rzz(i,j,km2) * CWt2m(1,k) & 
                 + Rzz(i,j,km1) * CWt2m(2,k) & 
                 + Rzz(i,j,k  ) * CWt2m(3,k) & 
                 + Rzz(i,j,kp1) * CWt2m(4,k)

            !compute vertical staggering of Rtt: thermo lvl -> mom lvl center point
            ttbz = Rtt(i,j,km2) * CWt2m(1,k) & 
                 + Rtt(i,j,km1) * CWt2m(2,k) & 
                 + Rtt(i,j,k  ) * CWt2m(3,k) & 
                 + Rtt(i,j,kp1) * CWt2m(4,k)

            !compute vertical derivative of Rzz: thermo lvl -> mom lvl center point
            dzrzz = Rzz(i,j,km2) * CDt2m(1,k) & 
                  + Rzz(i,j,km1) * CDt2m(2,k) & 
                  + Rzz(i,j,k  ) * CDt2m(3,k) & 
                  + Rzz(i,j,kp1) * CDt2m(4,k) 

            !compute vertical derivative of Rtt: thermo lvl -> mom lvl cntr point
            dzrtt = Rtt(i,j,km2) * CDt2m(1,k) & 
                  + Rtt(i,j,km1) * CDt2m(2,k) & 
                  + Rtt(i,j,k  ) * CDt2m(3,k) & 
                  + Rtt(i,j,kp1) * CDt2m(4,k)

            ! exact form of eqn 58 in SG notes
            Sol_rhs(i,j,k) = -invT_8*Rqq + dudx + dvdy + invT_8*dzrzz   &
                            + ubx*M_logJzu(i,j,k) + vby*M_logJzv(i,j,k) &
                            + invT_8*M_logJzq(i,j,k)* zzbz + dzrtt      &
                            + ttbz*(M_logJzq(i,j,k) - epsi_8)

! this exception at k=1 will be removed soon and km should be replaced with simply k-1
            if (k==1   ) then
               w3   = (Ver_idz_8%m(k) + (GVM%mc_Iz_8(i,j,k) - epsi_8)*Ver_wp_8%m(k))
               w5   = a*rhst_mid(i,j,k ) - b*rhst_dep(i,j,k )
               w5   = gama_bdf_8 * ( c * w5   + a*rhsw_mid(i,j,k ) - b*rhsw_dep(i,j,k ))
               w6(1)   = invT_8*( logT(i,j,k ) - (one-one/tots(i,j,k )) )
               Nwww = (tots(i,j,k)-one)*GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) &
                     -(tots(i,j,k)-one)*grav_8*(one-one/tots(i,j,k)) 
               Nttt = gama_bdf_8 * ( c * w6(1)   + Nwww )
               w7   = a*rhsf_mid(i,j,k ) - b*rhsf_dep(i,j,k )
               Sol_rhs(i,j,k)= -invT_8*Rqq + dudx + dvdy + ubx*GVM%mc_Ix_8(i,j,k) + vby*GVM%mc_Iy_8(i,j,k)&
                               + invT_8*w7*(Ver_idz_8%m(1)+GVM%mc_Iz_8(i,j,1)*Ver_wp_8%m(1)) + (w3*w5) - w3*Nttt(1)
            endif
         end do
         end do
      end do

! this exception at k=l_nk will also be removed soon
! precomputation of wka,wkb for special lower(l_nk) boundary condition
      wkf= 0. ; wka=0. ; wkb=0. ; k=l_nk
      do j=ds_j0,ds_jn
         do i=ds_i0,ds_in
            w1 = invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
            w4 = (tots(i,j,k)-one)*GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) &
               -(tots(i,j,k)-one)*grav_8*(one-one/tots(i,j,k))
            w2 = gama_bdf_8 * ( c * w1 + w4 )
            w1 = invT_8*( logT(i,j,k-1) - (one-one/tots(i,j,k-1)) )
            w4= (tots(i,j,k-1)-one)*GVM%mc_iJz_8(i,j,k-1)*(qt0(i,j,k)-qt0(i,j,k-1)) &
               -(tots(i,j,k-1)-one)*grav_8*(one-one/tots(i,j,k-1))
            w3 = gama_bdf_8 * ( c * w1 + w4 )
               
            wkf(i,j,1)= GVM%mc_css_H_8(i,j) * (Rt(i,j,l_nk)-Ver_wmstar_8(G_nk)*Rt(i,j,l_nk-1)&
         + invT_8*(Rzz(i,j,l_nk)-Ver_wmstar_8(G_nk)*Rzz(i,j,l_nk-1))) &
         - GVM%mc_css_H_8(i,j) * (w2-Ver_wmstar_8(G_nk)*w3)
         end do
      end do
      call HLT_split (1, 1, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( wkf(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

! Specific scope is here important for wka,wkb
      do j= j00, ds_jn
      do i= i00, inn
         wka(i,j) = - GVM%mc_Jx_8(i,j,l_nk) * &
                     Ver_wp_8%m(l_nk)*half*( wkf(i+1,j,1)*GVM%mc_iJz_8(i+1,j,l_nk) &
                                           + wkf(i  ,j,1)*GVM%mc_iJz_8(i  ,j,l_nk) )
      end do
      end do
      do j= j00, jnn
      do i= i00, ds_in
         wkb(i,j) = - GVM%mc_Jy_8(i,j,l_nk) * &
                     Ver_wp_8%m(l_nk)*half*( wkf(i,j+1,1)*GVM%mc_iJz_8(i,j+1,l_nk) &
                                           + wkf(i,j  ,1)*GVM%mc_iJz_8(i,j  ,l_nk) )
      end do
      end do
      
      k=l_nk
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            w1 = gama_8*(wkf(i,j,1)*GVM%mc_iJz_8(i,j,l_nk) - mu_8*half*wkf(i,j,1))
            w2= (wka (i,j)-wka (i-1,j))*geomh_invDXM_8(j) &
               + half * (GVM%mc_Ix_8(i,j,l_nk)*(wka(i,j)+wka(i-1,j)))
            w3= (wkb(i,j)*geomh_cyM_8(j)-wkb(i,j-1)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
               + half * (GVM%mc_Iy_8(i,j,l_nk)*(wkb(i,j)+wkb(i,j-1)))
            w4= w1*Ver_idz_8%m(l_nk) +(GVM%mc_Iz_8(i,j,l_nk)-epsi_8)*(Ver_wp_8%m(l_nk)*w1)
            Sol_rhs(i,j,k)= Sol_rhs(i,j,k)-w2-w3-w4
      end do
      end do

!!$      do k=1,l_nk
!!$         call statf_dm (Sol_rhs(1:,1:,k:k), 'ERHS', k, 'ELLI', 1,ubound(Sol_rhs,1),1,ubound(Sol_rhs,2),1,1,1+Glb_pil_w,1+Glb_pil_s,1,G_ni-Glb_pil_e,G_nj-Glb_pil_n,1,8)
!!$      end do
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine elliptic_rhs3rd_try1
