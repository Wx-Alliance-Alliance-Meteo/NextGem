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

      subroutine elliptic_rhs3rdVH ( F_dt_8, k0, k0t )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use HORgrid_options
      use geomh
      use glb_pil
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
      use stat_mpi
      use gmm_geof
      use bubble_options
      implicit none

      integer, intent(in) :: k0, k0t
      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k, HLT_np, HLT_start, HLT_end
      integer :: i00,inn,j00,jnn,dim,ub
      real, dimension(:,:,:), pointer :: tots, logT, logQ
      real(kind=REAL64), dimension(:,:,:), pointer :: t2u, v2u, t2v, u2v
      real(kind=REAL64), dimension(:,:,:), pointer :: dqz2u, dqz2v, dqz2w
      real(kind=REAL64) :: Rqq,tau_8,invT_8,a,b,c,barz,barzp
      real(kind=REAL64) :: w0,w1,w2,w3,w4,dudx,dvdy,ubx,vby
      real(kind=REAL64) :: ttbz, zzbz, dzrtt, difgz
      real(kind=REAL64) :: Ntttdz, Ntttbz, Nzzzdz, Nzzzbz,b1
      real(kind=REAL64), dimension(:,:,:), allocatable :: advf,advw,advt
      real(kind=REAL128) :: dzrzz
      real(kind=REAL64), dimension(1:6) :: Nttt, Nwww
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
      tots (1:l_ni,1:l_nj,-1:l_nk+1) => WS1(ub+1:) ; ub=ub+dim*(l_nk+3)
      logT (1:l_ni,1:l_nj,-1:l_nk+1) => WS1(ub+1:) ; ub=ub+dim*(l_nk+3)
      
      ub=0
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      t2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      v2u  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      t2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk
      u2v  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1_8(ub+1:) ; ub=ub+dim*l_nk

      !---vertical derivative of q to points u,v,w in rhs---
      dqz2u (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2v (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2w (l_minx:l_maxx, l_miny:l_maxy, -1:l_nk+1) => WS1_8(ub+1:); ub=ub+dim*(l_nk+3)
      
      allocate (advf (1:l_ni,1:l_nj,-4:l_nk+4),&
                advw (1:l_ni,1:l_nj,-4:l_nk+4),&
                advt (1:l_ni,1:l_nj,-4:l_nk+4))

      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               advf(i,j,k)= a*rhsf_mid(i,j,k) - b*rhsf_dep(i,j,k )
               advw(i,j,k)= a*rhsw_mid(i,j,k) - b*rhsw_dep(i,j,k )
               advt(i,j,k)= a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k )
            end do
         end do
      end do

      do k= 0,-4,-1
         do j=1, l_nj
            do i= 1, l_ni
               advf(i,j,k)=(VM3%ztht(i,j,k) - Ver_ext%t(k))*invT_8
               advw(i,j,k)=0.
            end do
         end do
      end do
      do k= l_nk+1,l_nk+4
         do j=1, l_nj
            do i= 1, l_ni
               advf(i,j,k)=(VM3%ztht(i,j,k) - Ver_ext%t(k))*invT_8
               advw(i,j,k)=0.
            end do
         end do
      end do

      k=0
      do j=1, l_nj
         do i= 1, l_ni
            b1=  qt0(i,j,k-1) * VS3m2t(1,k) & 
                +qt0(i,j,k  ) * VS3m2t(2,k) & 
                +qt0(i,j,k+1) * VS3m2t(3,k) & 
                +qt0(i,j,k+2) * VS3m2t(4,k)
            advt (i,j,k)= (log (ext_t(i,j,k)/Cstv_Tstr_8) - b1/(cpd_8*Cstv_Tstr_8))*invT_8
         end do
      end do
      do k= -1,-3,-1
      do j=1, l_nj
         do i= 1, l_ni
            b1=  qt0(i,j,k  ) * VS3m2t(1,k) & 
                +qt0(i,j,k+1) * VS3m2t(2,k) & 
                +qt0(i,j,k+2) * VS3m2t(3,k) & 
                +qt0(i,j,k+3) * VS3m2t(4,k)
                advt (i,j,k)= (log (ext_t(i,j,k)/Cstv_Tstr_8) - b1/(cpd_8*Cstv_Tstr_8))*invT_8
         end do
      end do
      end do
      do k= l_nk+1,l_nk+3
      do j=1, l_nj
         do i= 1, l_ni
            b1=  qt0(i,j,k-2) * VS3m2t(1,k) & 
                +qt0(i,j,k-1) * VS3m2t(2,k) & 
                +qt0(i,j,k  ) * VS3m2t(3,k) & 
                +qt0(i,j,k+1) * VS3m2t(4,k)
            advt (i,j,k)= (log (ext_t(i,j,k)/Cstv_Tstr_8) - b1/(cpd_8*Cstv_Tstr_8))*invT_8
         end do
      end do
      end do

      do k=-1, l_nk+1
         do j=1, l_nj
            do i= 1, l_ni
               tots(i,j,k)= ext_t(i,j,k)/Cstv_Tstr_8
               logT(i,j,k)= log(tots(i,j,k))
            end do
         end do
      end do
      
      call prerhs3rdVH (t2u, v2u, t2v, u2v, dqz2u, dqz2v, dqz2w,&
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
            Nuu(i,j,k) =  t2u(i,j,k)*dqz2u(i,j,k) &
                        - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut0(i,j,k) ) * v2u(i,j,k)
         end do
         end do         
        
         do j= j00, jnn
         do i= ds_i0, ds_in
            Nvv(i,j,k)= t2v(i,j,k)*dqz2v(i,j,k) &
                        + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u2v(i,j,k)) * u2v(i,j,k)
         end do
         end do

         do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            w0= tots(i,j,k)-one
            w2= advw(i,j,k)
            w4= w0*(dqz2w(i,j,k) - grav_8*(one-one/tots(i,j,k)))
            Rww(i,j,k)= w2 - w4
         end do
         end do
      end do

      do k=-1, l_nk+1
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
            w0= tots(i,j,k)-one
            w1= advt(i,j,k)
            w2= advw(i,j,k)
            w3= invT_8*( logT(i,j,k) - (one-one/tots(i,j,k)) )
            w4= w0*(dqz2w(i,j,k) - grav_8*(one-one/tots(i,j,k)))
            Rtt(i,j,k)= gama_bdf_8 * ( c*(w1 - w3) + (w2 - w4))
            Rzz(i,j,k)= advf(i,j,k) - invT_8*(VM3%ztht(i,j,k)-Ver_ext%t(k))
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
            zzbz = Rzz(i,j,k-2) * VS3t2m(1,k) & 
                 + Rzz(i,j,k-1) * VS3t2m(2,k) & 
                 + Rzz(i,j,k  ) * VS3t2m(3,k) & 
                 + Rzz(i,j,k+1) * VS3t2m(4,k)

            !compute vertical staggering of Rtt: thermo lvl -> mom lvl center point
            ttbz = Rtt(i,j,k-2) * VS3t2m(1,k) & 
                 + Rtt(i,j,k-1) * VS3t2m(2,k) & 
                 + Rtt(i,j,k  ) * VS3t2m(3,k) & 
                 + Rtt(i,j,k+1) * VS3t2m(4,k)

            !compute vertical derivative of Rzz: thermo lvl -> mom lvl center point
            dzrzz = Rzz(i,j,k-2) * VD3t2m(1,k) & 
                  + Rzz(i,j,k-1) * VD3t2m(2,k) & 
                  + Rzz(i,j,k  ) * VD3t2m(3,k) & 
                  + Rzz(i,j,k+1) * VD3t2m(4,k) 

            !compute vertical derivative of Rtt: thermo lvl -> mom lvl cntr point
            dzrtt = Rtt(i,j,k-2) * VD3t2m(1,k) & 
                  + Rtt(i,j,k-1) * VD3t2m(2,k) & 
                  + Rtt(i,j,k  ) * VD3t2m(3,k) & 
                  + Rtt(i,j,k+1) * VD3t2m(4,k)

            ! exact form of eqn 58 in SG notes
            Sol_rhs(i,j,k) = -invT_8*Rqq + dudx + dvdy + invT_8*dzrzz   &
                            + ubx*M_logJzu(i,j,k) + vby*M_logJzv(i,j,k) &
                            + invT_8*M_logJzq(i,j,k)* zzbz + dzrtt      &
                            + ttbz*(M_logJzq(i,j,k) - epsi_8)
         end do
         end do
      end do

      deallocate (advf,advw,advt)
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine elliptic_rhs3rdVH
      
