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
!**s/r elliptic_rhs - Compute right hand side of the elliptic problem

      subroutine elliptic_rhs ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use gem_options
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
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      integer :: i0,in,j0,jn,km,i00,inn,j00,jnn
      real, dimension(:,:,:), pointer :: tots, logT, logQ
      real(kind=REAL64) :: tau_8,invT_8,Ldiv,Ndiv,a,ray,b,c,barz,barzp
      real(kind=REAL64) :: w1,w2,w2km,w3,w4,w5,w5km,w6,w6km
      real(kind=REAL64) :: Nwww,Nwkm,Nttt,Ntkm,t_interp,u_interp,v_interp
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
      i00= i0-1 ; inn= in
      j00= j0-1 ; jnn= jn
      if (.not.Grd_yinyang_L) then
         i00= i0-1+min(pil_w,1)
         inn= in  -min(pil_e,1)
         j00= j0-1+min(pil_s,1)
         jnn= jn  -min(pil_n,1)
      endif
      
      tau_8  = (2.d0 * F_dt_8 ) / 3.d0
      invT_8 = 1.d0/tau_8
      a      = 4.d0*invT_8/3.d0
      b      =      invT_8/3.d0
      c      = grav_8 * tau_8
      ray    = Dcst_rayt_8**2

      tots (1:l_ni,1:l_nj,1:l_nk) => WS1(               1:)
      logT (1:l_ni,1:l_nj,1:l_nk) => WS1(l_ni*l_nj*l_nk+1:)

      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               tots(i,j,k)= tt0(i,j,k)/Cstv_Tstr_8
               logT(i,j,k)= log(tots(i,j,k))
            end do
         end do
      end do

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
         do j= Adz_j0, Adz_jn
         do i= Adz_i0, Adz_in
            Rww(i,j,k) = a*rhsw_mid(i,j,k) - b*rhsw_dep(i,j,k)
         end do
         end do
         if (l_north) Rvv(1+pil_w:l_ni-pil_e,l_nj-pil_n,k)=rhsv(1+pil_w:l_ni-pil_e,l_nj-pil_n,k)
         if (l_south) Rvv(1+pil_w:l_ni-pil_e,pil_s,k)=rhsv(1+pil_w:l_ni-pil_e,pil_s,k)
         if (l_east)  Ruu(l_ni-pil_e,1+pil_s:l_nj-pil_n,k)=rhsu(l_ni-pil_e,1+pil_s:l_nj-pil_n,k)
         if (l_west)  Ruu(pil_w,1+pil_s:l_nj-pil_n,k)=rhsu(pil_w,1+pil_s:l_nj-pil_n,k)
         do j= j0, jn
         do i= i00, inn
            barz = Ver_wp_8%m(k)*tt0(i  ,j,k)+Ver_wm_8%m(k)*tt0(i  ,j,km)
            barzp= Ver_wp_8%m(k)*tt0(i+1,j,k)+Ver_wm_8%m(k)*tt0(i+1,j,km)
            t_interp = (barz+barzp)*half/Cstv_Tstr_8
            v_interp = 0.25d0*(vt0(i  ,j,k)+vt0(i  ,j-1,k)+&
                               vt0(i+1,j,k)+vt0(i+1,j-1,k))
            Nuu(i,j,k)= (t_interp-one) * ( qt0(i+1,j,k) - qt0(i,j,k) ) * geomh_invDX_8(j) &
                           - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v_interp &
                           - (t_interp-one) * GVM%mc_Jx_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i+1,j,k+1)-qt0(i+1,j,k ))*GVM%mc_iJz_8(i+1,j,k )   &
                                                 +(qt0(i  ,j,k+1)-qt0(i  ,j,k ))*GVM%mc_iJz_8(i  ,j,k ) ) &
                            +Ver_wm_8%m(k)*half*( (qt0(i+1,j,k  )-qt0(i+1,j,km))*GVM%mc_iJz_8(i+1,j,km)   &
                                                 +(qt0(i  ,j,k  )-qt0(i  ,j,km))*GVM%mc_iJz_8(i  ,j,km) ) )
         end do
         end do         
         do j= j00, jnn
         do i= i0, in
            barz  = Ver_wp_8%m(k)*tt0(i,j  ,k)+Ver_wm_8%m(k)*tt0(i,j  ,km)
            barzp = Ver_wp_8%m(k)*tt0(i,j+1,k)+Ver_wm_8%m(k)*tt0(i,j+1,km)
            t_interp = (barz+barzp)*half/Cstv_Tstr_8
            u_interp = 0.25d0*(ut0(i,j,k)+ut0(i-1,j,k)+ut0(i,j+1,k)+ut0(i-1,j+1,k))
               Nvv(i,j,k) = (t_interp-one) * ( qt0(i,j+1,k) - qt0(i,j,k) ) * geomh_invDY_8 &
                           + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp &
                           - (t_interp-one) * GVM%mc_Jy_8(i,j,k) * ( &
                             Ver_wp_8%m(k)*half*( (qt0(i,j+1,k+1)-qt0(i,j+1,k ))*GVM%mc_iJz_8(i,j+1,k )   &
                                                 +(qt0(i,j  ,k+1)-qt0(i,j  ,k ))*GVM%mc_iJz_8(i,j  ,k ) ) &
                            +Ver_wm_8%m(k)*half*( (qt0(i,j+1,k  )-qt0(i,j+1,km))*GVM%mc_iJz_8(i,j+1,km)   &
                                                 +(qt0(i,j  ,k  )-qt0(i,j  ,km))*GVM%mc_iJz_8(i,j  ,km) ) )
         end do
         end do
      end do
      call HLT_split (1, 2*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( Ruu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      do k=2, l_nk-1
         km=max(k-1,1)
         do j= j0, jn
         do i= i0, in
            w1   = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k )
            w2   = a*rhsf_mid(i,j,k ) - b*rhsf_dep(i,j,k ) - invT_8*(GVM%ztht_8(i,j,k )-Ver_z_8%t(k ))
            w2km = a*rhsf_mid(i,j,km) - b*rhsf_dep(i,j,km) - invT_8*(GVM%ztht_8(i,j,km)-Ver_z_8%t(km))
            w3   = (Ver_idz_8%m(k) + (GVM%mc_Iz_8(i,j,k)- epsi_8)*Ver_wp_8%m(k))
            w4   = (Ver_idz_8%m(k) - (GVM%mc_Iz_8(i,j,k)- epsi_8)*Ver_wm_8%m(k))
            w5   = a*rhst_mid(i,j,k ) - b*rhst_dep(i,j,k )
            w5km = a*rhst_mid(i,j,km) - b*rhst_dep(i,j,km)
            w5   = gama_bdf_8 * ( c * w5   + a*rhsw_mid(i,j,k ) - b*rhsw_dep(i,j,k ))
            w5km = gama_bdf_8 * ( c * w5km + a*rhsw_mid(i,j,km) - b*rhsw_dep(i,j,km))

            w6   = invT_8*( logT(i,j,k ) - (one-one/tots(i,j,k )) )
            w6km = invT_8*( logT(i,j,km) - (one-one/tots(i,j,km)) )
            Nwww = (tots(i,j,k)-one)*GVM%mc_iJz_8(i,j,k)*(qt0(i,j,k+1)-qt0(i,j,k)) &
                  -(tots(i,j,k)-one)*grav_8*(one-one/tots(i,j,k)) 
            Nwkm = (tots(i,j,km)-one)*GVM%mc_iJz_8(i,j,km)*(qt0(i,j,k)-qt0(i,j,k-1)) &
                  -(tots(i,j,km)-one)*grav_8*(one-one/tots(i,j,km)) 
            Nttt = gama_bdf_8 * ( c * w6   + Nwww )
            Ntkm = gama_bdf_8 * ( c * w6km + Nwkm )
            Nww(i,j,k) = Nwww
            Ntt(i,j,k) = Nttt
            Rzz(i,j,k) = w2
            Rtt(i,j,k) = w5
            Ldiv= (Ruu(i,j,k)-Ruu(i-1,j,k))*geomh_invDXM_8(j) &
               + (Rvv(i,j  ,k)*geomh_cyM_8(j)- &
                  Rvv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j)  &
               + half*(GVM%mc_Ix_8(i,j,k)*(Ruu(i,j,k)+Ruu(i-1,j,k)) &
               + GVM%mc_Iy_8(i,j,k)*(Rvv(i,j,k)+Rvv(i,j-1,k)) )
            Ndiv=(Nuu(i,j,k)-Nuu(i-1,j,k))*geomh_invDXM_8(j) &
               + (Nvv(i,j,k)*geomh_cyM_8(j)- &
                  Nvv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
               + half*( GVM%mc_Ix_8(i,j,k)*(Nuu(i,j,k)+Nuu(i-1,j,k)) &
               + GVM%mc_Iy_8(i,j,k)*(Nvv(i,j,k)+Nvv(i,j-1,k)) )

            RHS_sol(i,j,k) = ray*(Ldiv - invT_8*w1 + invT_8*( (w2-w2km)*Ver_idz_8%m(k) + GVM%mc_Iz_8(i,j,k)*(Ver_wp_8%m(k)*w2+Ver_wm_8%m(k)*w2km)) +&
                                 w3 * w5 - w4 * w5km - Ndiv - (w3 * Nttt - w4 * Ntkm ))
         end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine elliptic_rhs
