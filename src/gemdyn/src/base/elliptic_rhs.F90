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
      integer :: km, km3, km2, km1, kp1,kp2,i00,inn,j00,jnn,dim,ub
      real, dimension(:,:,:), pointer :: wkf
      real, dimension(:,:,:), pointer :: tots, logT, logQ, Rt!, t2u, v2u, t2v, u2v
      real(kind=REAL64), dimension(:,:,:), pointer  :: t2u, v2u, t2v, u2v
      real(kind=REAL64), dimension(:,:,:), pointer  :: dqz2u, dqz2v, dqz2w
      real(kind=REAL64) :: Rqq,tau_8,invT_8,a,b,c,barz,barzp
      real(kind=REAL64) :: w0,w1,w2,w3,w4,w5,w7, dudx,dvdy,ubx,vby
      real(kind=REAL64) :: dqdx, dqdy, dzrtt, ttbz, zzbz
      real(kind=REAL64) :: Ntttdz, Ntttbz, Nzzzdz, Nzzzbz
      real(kind=REAL64) :: wka(l_minx:l_maxx,l_miny:l_maxy), wkb(l_minx:l_maxx,l_miny:l_maxy)
      real(kind=REAL128) :: dzrzz
      real(kind=REAL64), dimension(1:6) :: Nttt, Nwww, w6
      real(kind=REAL64), parameter :: zero=0.d0, one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
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
      
      !---added my code to Michel's existing to compute vertical derivative of q to u,q,w points---
      call tt2wnd(t2u, v2u, t2v, u2v, dqz2u, dqz2v, dqz2w, l_minx,l_maxx,l_miny,l_maxy,G_nk)
      
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
              w4= w0*(dqz2w(i,j,k) - grav_8*(one-one/tots(i,j,k)))

              Rtt(i,j,k)= w1 - w3
              Rww(i,j,k)= w2 - w4

              Rtt(i,j,k)= gama_bdf_8 * ( c*Rtt(i,j,k) + Rww(i,j,k) )
              Rzz(i,j,k)= a*rhsf_mid(i,j,k) - b*rhsf_dep(i,j,k ) - invT_8*(GVM%ztht_8(i,j,k)-Ver_z_8%t(k))
            
              Rt(i,j,k) = gama_bdf_8 * ( c*w1 + w2 ) ! special one for wkf below
           end do
         end do
      end do
      
      call HLT_split (1, 2*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( Ruu(l_minx,l_miny,HLT_start),&
                 l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      do k=1, l_nk
         Ruu(:,:,k)= Ruu(:,:,k) - Nuu(:,:,k)
         Rvv(:,:,k)= Rvv(:,:,k) - Nvv(:,:,k)
      end do

      do k=1, l_nk

         km1=max(k-1,1)
         km2=max(k-2,1)
         km3=max(k-3,1)
         kp1=min(k+1,G_nk)
         kp2=min(k+2,G_nk)

         do j= ds_j0, ds_jn
           do i= ds_i0, ds_in

            Rqq = a*rhsc_mid(i,j,k ) - b*rhsc_dep(i,j,k )

            !dx(Ru): horizontal derivative of Ru -> use Hder function
            dudx = Hderiv(Ruu(i-2,j,k), Ruu(i-1,j,k), Ruu(i,j,k), &
                          Ruu(i+1,j,k), Ruu(i+2,j,k), Ruu(i+3,j,k), geomh_invDXM_8(j))

            !dy(Rv): horizonal derivative of Rv -> use Hder function
            dvdy = Hderiv8(Rvv(i,j-2,k)*geomh_cyM_8(j-2), &
                           Rvv(i,j-1,k)*geomh_cyM_8(j-1), &
                           Rvv(i,j  ,k)*geomh_cyM_8(j)  , &
                           Rvv(i,j+1,k)*geomh_cyM_8(j+1), &
                           Rvv(i,j+2,k)*geomh_cyM_8(j+2), &
                           Rvv(i,j+3,k)*geomh_cyM_8(j+3), &
                           geomh_invDYM_8(j))

            !\overbar{Ruu}: horizonal staggering -> use Hder func
            ubx = Hstag(Ruu(i-2,j,k), Ruu(i-1,j,k), Ruu(i  ,j,k), &
                        Ruu(i+1,j,k), Ruu(i+2,j,k), Ruu(i+3,j,k))

            !\overbar{Rvv}: horizonal staggering -> use Hder func
            vby = Hstag(Rvv(i-2,j,k), Rvv(i-1,j,k), Rvv(i  ,j,k), &
                        Rvv(i+1,j,k), Rvv(i+2,j,k), Rvv(i+3,j,k))

            !compute vertical staggering of Rzz: thermo lvl -> mom lvl center point
            zzbz = Rzz(i,j,km3) * QWt2m(1,k) &
                 + Rzz(i,j,km2) * QWt2m(2,k) & 
                 + Rzz(i,j,km1) * QWt2m(3,k) & 
                 + Rzz(i,j,k  ) * QWt2m(4,k) & 
                 + Rzz(i,j,kp1) * QWt2m(5,k) & 
                 + Rzz(i,j,kp2) * QWt2m(6,k)  

            !compute vertical staggering of Rtt: thermo lvl -> mom lvl center point
            ttbz = Rtt(i,j,km3) * QWt2m(1,k) & 
                 + Rtt(i,j,km2) * QWt2m(2,k) & 
                 + Rtt(i,j,km1) * QWt2m(3,k) & 
                 + Rtt(i,j,k  ) * QWt2m(4,k) & 
                 + Rtt(i,j,kp1) * QWt2m(5,k) & 
                 + Rtt(i,j,kp2) * QWt2m(6,k)  

            !compute vertical derivative of Rzz: thermo lvl -> mom lvl center point
            !d1  = (Rzz(i,j,k)-Rzz(i,j,km))*Ver_idz_8%m(k) ! very sensitive... hence temporary storage in quad precicion
            dzrzz = Rzz(i,j,km3) * QDt2m(1,k) &
                  + Rzz(i,j,km2) * QDt2m(2,k) & 
                  + Rzz(i,j,km1) * QDt2m(3,k) & 
                  + Rzz(i,j,k  ) * QDt2m(4,k) & 
                  + Rzz(i,j,kp1) * QDt2m(5,k) & 
                  + Rzz(i,j,kp2) * QDt2m(6,k)  

            !compute vertical derivative of Rtt: thermo lvl -> mom lvl cntr point
            dzrtt = Rtt(i,j,km3) * QDt2m(1,k) & 
                  + Rtt(i,j,km2) * QDt2m(2,k) & 
                  + Rtt(i,j,km1) * QDt2m(3,k) & 
                  + Rtt(i,j,k  ) * QDt2m(4,k) & 
                  + Rtt(i,j,kp1) * QDt2m(5,k) & 
                  + Rtt(i,j,kp2) * QDt2m(6,k)  

            ! exact form of eqn 58 in SG notes
            Sol_rhs(i,j,k) = -invT_8*Rqq + dudx + dvdy + invT_8*dzrzz   &
                            + ubx*M_logJzu(i,j,k) + vby*M_logJzv(i,j,k) &
                            + invT_8*M_logJzq(i,j,k)* zzbz              &
                            + dzrtt                                     &
                            + ttbz*(M_logJzq(i,j,k) - epsi_8) 

! this exception at k=1 will be removed soon and km should be replaced with simply k-1

            !still needed update this 
            if (k==1   ) then
               !w3   = (Ver_idz_8%m(k) + (GVM%mc_Iz_8(i,j,k) - epsi_8)*Ver_wp_8%m(k))

               w5   = a*rhst_mid(i,j,k ) - b*rhst_dep(i,j,k )
               w5   = gama_bdf_8 * ( c * w5   + a*rhsw_mid(i,j,k ) - b*rhsw_dep(i,j,k ))

               !need vertical operators of Nt, so compute this 6 times
               w6(1)   = invT_8*( logT(i,j,km3) - (one-one/tots(i,j,km3)) )
               w6(2)   = invT_8*( logT(i,j,km2) - (one-one/tots(i,j,km2)) )
               w6(3)   = invT_8*( logT(i,j,km1) - (one-one/tots(i,j,km1)) )
               w6(4)   = invT_8*( logT(i,j,k  ) - (one-one/tots(i,j,k  )) )
               w6(5)   = invT_8*( logT(i,j,kp1) - (one-one/tots(i,j,kp1)) )
               w6(6)   = invT_8*( logT(i,j,kp2) - (one-one/tots(i,j,kp2)) )
               !w7   = a*rhsf_mid(i,j,k ) - b*rhsf_dep(i,j,k ) !why is it missing -(z-\zeta)?

               !replace with 5th order der
               Nwww(1) = (tots(i,j,km3)-one)*dqz2w(i,j,km3)  &
                       - (tots(i,j,km3)-one)*grav_8*(one-one/tots(i,j,km3)) 
               Nwww(2) = (tots(i,j,km2)-one)*dqz2w(i,j,km2)  &
                       - (tots(i,j,km2)-one)*grav_8*(one-one/tots(i,j,km2)) 
               Nwww(3) = (tots(i,j,km1)-one)*dqz2w(i,j,km1)  &
                       - (tots(i,j,km1)-one)*grav_8*(one-one/tots(i,j,km1)) 
               Nwww(4) = (tots(i,j,k  )-one)*dqz2w(i,j,k )  &
                       - (tots(i,j,k  )-one)*grav_8*(one-one/tots(i,j,k  )) 
               Nwww(5) = (tots(i,j,kp1)-one)*dqz2w(i,j,kp1)  &
                       - (tots(i,j,kp1)-one)*grav_8*(one-one/tots(i,j,kp1)) 
               Nwww(6) = (tots(i,j,kp2)-one)*dqz2w(i,j,kp2)  &
                       - (tots(i,j,kp2)-one)*grav_8*(one-one/tots(i,j,kp2)) 

               Nttt(1) = gama_bdf_8 * ( c * w6(1)  + Nwww(1) )
               Nttt(2) = gama_bdf_8 * ( c * w6(2)  + Nwww(2) )
               Nttt(3) = gama_bdf_8 * ( c * w6(3)  + Nwww(3) )
               Nttt(4) = gama_bdf_8 * ( c * w6(4)  + Nwww(4) )
               Nttt(5) = gama_bdf_8 * ( c * w6(5)  + Nwww(5) )
               Nttt(6) = gama_bdf_8 * ( c * w6(6)  + Nwww(6) )

               !remove grouping of w3 and w7 to 2 seperate operations
               !Note: QW and QD weights default to linear at k=1 

               !vertical interpolation of Nttt
               Ntttbz = Nttt(1) * QWt2m(1,k) & 
                      + Nttt(2) * QWt2m(2,k) & 
                      + Nttt(3) * QWt2m(3,k) & 
                      + Nttt(4) * QWt2m(4,k) & 
                      + Nttt(5) * QWt2m(5,k) & 
                      + Nttt(6) * QWt2m(6,k)

               !vertical derivative  of Nttt
               Ntttdz = Nttt(1) * QDt2m(1,k) & 
                      + Nttt(2) * QDt2m(2,k) & 
                      + Nttt(3) * QDt2m(3,k) & 
                      + Nttt(4) * QDt2m(4,k) & 
                      + Nttt(5) * QDt2m(5,k) & 
                      + Nttt(6) * QDt2m(6,k)

               !vertical derivative of Rz
               Nzzzdz = (a*rhsf_mid(i,j,km3) - b*rhsf_mid(i,j,km3)) * QDt2m(1,k) & 
                      + (a*rhsf_mid(i,j,km2) - b*rhsf_mid(i,j,km2)) * QDt2m(2,k) & 
                      + (a*rhsf_mid(i,j,km1) - b*rhsf_mid(i,j,km1)) * QDt2m(3,k) & 
                      + (a*rhsf_mid(i,j,k  ) - b*rhsf_mid(i,j,k  )) * QDt2m(4,k) & 
                      + (a*rhsf_mid(i,j,kp1) - b*rhsf_mid(i,j,kp1)) * QDt2m(5,k) & 
                      + (a*rhsf_mid(i,j,kp2) - b*rhsf_mid(i,j,kp2)) * QDt2m(6,k) 
         
               !vertical interpolation of Rz 
               Nzzzbz = (a*rhsf_mid(i,j,km3) - b*rhsf_mid(i,j,km3)) * QWt2m(1,k) & 
                      + (a*rhsf_mid(i,j,km2) - b*rhsf_mid(i,j,km2)) * QWt2m(2,k) & 
                      + (a*rhsf_mid(i,j,km1) - b*rhsf_mid(i,j,km1)) * QWt2m(3,k) & 
                      + (a*rhsf_mid(i,j,k  ) - b*rhsf_mid(i,j,k  )) * QWt2m(4,k) & 
                      + (a*rhsf_mid(i,j,kp1) - b*rhsf_mid(i,j,kp1)) * QWt2m(5,k) & 
                      + (a*rhsf_mid(i,j,kp2) - b*rhsf_mid(i,j,kp2)) * QWt2m(6,k) 

               Sol_rhs(i,j,k) = -invT_8*Rqq + dudx + dvdy                   &
                              +  ubx*M_logJzu(i,j,k) + vby*M_logJzu(i,j,k)  &
                              +  Ntttdz + Ntttbz*(M_logJzq(i,j,k) - epsi_8) &
                              +  invT_8*(Nzzzdz + M_logJzq(i,j,1)*Nzzzbz)
                              !+  invT_8*w7*(Ver_idz_8%m(1)+M_logJzq(i,j,1)*Ver_wp_8%m(1)) + (w3*w5) - w3*Nttt
            endif
          end do
         end do
      end do
      
!---DID NOT UPDATE LOWER BOUNDARY CONDITION---

! this exception at k=l_nk will also be removed soon
! precomputation of wka,wkb for special lower(l_nk) boundary condition
      wkf= 0. ; wka=0. ; wkb=0. ; k=l_nk
      do j=ds_j0,ds_jn
         do i=ds_i0,ds_in

            w1 = invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
            w4 = (tots(i,j,k)-one)*dqz2w(i,j,k) & !derivative updated
               - (tots(i,j,k)-one)*grav_8*(one-one/tots(i,j,k))

            w2 = gama_bdf_8 * ( c * w1 + w4 )

            w1 = invT_8*( logT(i,j,k-1) - (one-one/tots(i,j,k-1)) )
            w4= (tots(i,j,k-1)-one)*dqz2w(i,j,k-1) & !derivative updated
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
!
!     ---------------------------------------------------------------
!
      return
      include 'H5th_ope.inc'
      end subroutine elliptic_rhs
