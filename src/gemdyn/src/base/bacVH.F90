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

!**s/r bac - backsubstitution: obtain new values for variables:u,v,w,t,q,zd
!                              using Sol_lhs

      subroutine bacVH ( F_dt_8 )
      use dyn_fisl_options
      use HORgrid_options
      use gem_options
      use geomh
      use sol_mem
      use tdpack
      use gmm_vt0
      use mem_tstp
      use glb_ld
      use cstv
      use ver
      use metric
      use ctrl
      use stat_mpi
      use, intrinsic :: iso_fortran_env
      implicit none

      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: i, j, k,ni,nj,ub,dim
      integer :: HLT_start, HLT_end, HLT_np
      real(kind=REAL64) :: w5, tau_8, invT_8, Buoy, a,b,b1,q,Qq2,qbz2,Qw2
      real(kind=REAL64), parameter :: one=1.d0
      real(kind=REAL64), dimension(:,:,:), pointer :: dqz2u, dqz2v, dqz2w
!
!     ---------------------------------------------------------------
!
      if (Ctrl_testcases_adv_L) then
!!$omp single
         call canonical_cases ("BAC")
!!$omp end single
         return
      end if

      print*, 'BAC avant update'
      i=l_ni/2
      j=l_nj/2+1
      do k=1, l_nk
         Qq2= (qt0(i,j,k-1) * VD3m2t(1,k)&
              +qt0(i,j,k  ) * VD3m2t(2,k)&
              +qt0(i,j,k+1) * VD3m2t(3,k)&
              +qt0(i,j,k+2) * VD3m2t(4,k))*M_iJzq(i,j,k)
         qbz2= qt0(i,j,k-1) * VS3m2t(1,k)&
              +qt0(i,j,k  ) * VS3m2t(2,k)&
              +qt0(i,j,k+1) * VS3m2t(3,k)&
              +qt0(i,j,k+2) * VS3m2t(4,k)
         Qw2= Qq2 - mu_8*qbz2
         write(6,'(i3,4(1pe14.6))') k,qt0(i,j,k  ),Rtt(i,j,k) , gama_bdf_8*Qw2,abs(Rtt(i,j,k) -gama_bdf_8*Qw2)/abs(Rtt(i,j,k))
      end do
      
      ub=0
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      dqz2u (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2v (l_minx:l_maxx, l_miny:l_maxy, 1:l_nk) => WS1_8(ub+1:); ub=ub+dim*l_nk
      dqz2w (l_minx:l_maxx, l_miny:l_maxy, -1:l_nk+1) => WS1_8(ub+1:); ub=ub+dim*(l_nk+3)
      
      tau_8 = (2.d0*F_dt_8) / 3.d0 
      invT_8= one/tau_8

      do k= 0,-3,-1
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            q= (VM3%zmom(i,j,k+1)-VM3%zmom(i,j,k))/(VM3%zmom(i,j,k+2)-VM3%zmom(i,j,k+1))
            Sol_lhs(i,j,k)= (1+q)*Sol_lhs(i,j,k+1) - q*Sol_lhs(i,j,k+2)
      end do
      end do
      end do
      do k= l_nk+1, l_nk+4
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            q= (VM3%zmom(i,j,k)-VM3%zmom(i,j,k-1))/(VM3%zmom(i,j,k-1)-VM3%zmom(i,j,k-2))
            Sol_lhs(i,j,k)= (1+q)*Sol_lhs(i,j,k-1) - q*Sol_lhs(i,j,k-2)
      end do
      end do
      end do
      qt0(ds_i0:ds_in,ds_j0:ds_jn,:) = Sol_lhs(ds_i0:ds_in,ds_j0:ds_jn,:)
      call mirror ()
      call HLT_split (-3, l_nk+4, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo (qt0(l_minx,l_miny,HLT_start),l_minx,l_maxx, l_miny,l_maxy, HLT_np,-1)
      Sol_lhs = qt0

      print*, 'BAC apres update'
      i=l_ni/2
      j=l_nj/2+1
      do k=1, l_nk
         Qq2= (qt0(i,j,k-1) * VD3m2t(1,k)&
              +qt0(i,j,k  ) * VD3m2t(2,k)&
              +qt0(i,j,k+1) * VD3m2t(3,k)&
              +qt0(i,j,k+2) * VD3m2t(4,k))*M_iJzq(i,j,k)
         qbz2= qt0(i,j,k-1) * VS3m2t(1,k)&
              +qt0(i,j,k  ) * VS3m2t(2,k)&
              +qt0(i,j,k+1) * VS3m2t(3,k)&
              +qt0(i,j,k+2) * VS3m2t(4,k)
         Qw2= Qq2 - mu_8*qbz2
         write(6,'(i3,4(1pe14.6))') k,qt0(i,j,k  ),Rtt(i,j,k) , gama_bdf_8*Qw2,abs(Rtt(i,j,k) -gama_bdf_8*Qw2)/abs(Rtt(i,j,k))
      end do

      call dqdz3rd ( Sol_lhs, dqz2u , dqz2v, dqz2w, &
                     l_minx,l_maxx,l_miny,l_maxy,G_nk,-3, l_nk+4)

      do k= 0, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               b1=  Sol_lhs(i,j,k-1) * VS3m2t(1,k) & 
                   +Sol_lhs(i,j,k  ) * VS3m2t(2,k) & 
                   +Sol_lhs(i,j,k+1) * VS3m2t(3,k) & 
                   +Sol_lhs(i,j,k+2) * VS3m2t(4,k)
               dqz2w(i,j,k) = dqz2w(i,j,k) - mu_8*b1
            end do
         end do
      end do
      do j= ds_j0, ds_jn
         do i= ds_i0, ds_in
            b1= Sol_lhs(i,j,-1) * VS3m2t(1,-1) & 
               +Sol_lhs(i,j, 0) * VS3m2t(2,-1) & 
               +Sol_lhs(i,j, 1) * VS3m2t(3,-1) & 
               +Sol_lhs(i,j, 2) * VS3m2t(4,-1)                      
            dqz2w(i,j,-1) = dqz2w(i,j,-1) - mu_8*b1
            b1= Sol_lhs(i,j,l_nk-1) * VS3m2t(1,l_nk+1) & 
               +Sol_lhs(i,j,l_nk  ) * VS3m2t(2,l_nk+1) & 
               +Sol_lhs(i,j,l_nk+1) * VS3m2t(3,l_nk+1) & 
               +Sol_lhs(i,j,l_nk+2) * VS3m2t(4,l_nk+1)
            dqz2w(i,j,l_nk+1) = dqz2w(i,j,l_nk+1) - mu_8*b1
         end do
      end do

      a = 4.d0*invT_8/3.d0
      b =      invT_8/3.d0
!!$omp do collapse(2)
      do k=ds_k0, l_nk
         do j= ds_j0, ds_jn
            do i= ds_i0, l_niu-pil_e
               ut0(i,j,k) = tau_8*(Ruu(i,j,k) - dqz2u(i,j,k))
            end do
         end do
         do j= ds_j0, l_njv-pil_n
            do i= ds_i0, ds_in
               vt0(i,j,k) = tau_8*(Rvv(i,j,k) - dqz2v(i,j,k))
            end do
         end do
         do j= ds_j0, ds_jn
            do i= ds_i0, ds_in
               b1= Sol_lhs(i,j,k-1) * VS3m2t(1,k)& 
                  +Sol_lhs(i,j,k  ) * VS3m2t(2,k)& 
                  +Sol_lhs(i,j,k+1) * VS3m2t(3,k)& 
                  +Sol_lhs(i,j,k+2) * VS3m2t(4,k)
               q= (Sol_lhs(i,j,k-1) * VD3m2t(1,k)&
                  +Sol_lhs(i,j,k  ) * VD3m2t(2,k)&
                  +Sol_lhs(i,j,k+1) * VD3m2t(3,k)&
                  +Sol_lhs(i,j,k+2) * VD3m2t(4,k))*M_iJzq(i,j,k)
               wt0 (i,j,k) = tau_8*(Rtt(i,j,k) - gama_bdf_8*dqz2w(i,j,k))
               zdt0(i,j,k) = (Rzz(i,j,k) + wt0(i,j,k))
               Buoy = q + wt0(i,j,k)*invT_8 - Rww(i,j,k)
               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy / grav_8 )
!!$
!!$               w5= a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k)&
!!$                   - invT_8* ( log(tt0(i,j,k)/Cstv_Tstr_8) - (one-Cstv_Tstr_8/tt0(i,j,k) ))
!!$               Buoy = b1/(cpd_8*Cstv_Tstr_8) + tau_8*(w5-mu_8*wt0(i,j,k))
!!$                  if ((i==l_ni/2).and.(j==l_nj/2+1)) print*, Rtt(i,j,k),gama_bdf_8*dqz2w(i,j,k),a*rhst_mid(i,j,k) - b*rhst_dep(i,j,k),invT_8* ( log(tt0(i,j,k)/Cstv_Tstr_8) - (one-Cstv_Tstr_8/tt0(i,j,k) ))
!!$               tt0(i,j,k) = Cstv_Tstr_8 / (one - Buoy)
            end do
         end do
      end do
      
      call mirror ()

      call HLT_split (1, l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( tt0(l_minx,l_miny,HLT_start),l_minx,l_maxx, l_miny,l_maxy, HLT_np,-1)
      ext_t(:,:,1:l_nk) = tt0(:,:,1:l_nk)
      do k=-4,0
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               w5= (VM3%ztht(i,j,k)-VM3%ztht(i,j,2))/(VM3%ztht(i,j,2)-VM3%ztht(i,j,1))
               ext_t(i,j,k)= tt0(i,j,2)*(1+w5) -w5*tt0(i,j,1)
            end do
         end do
      end do
      do k= l_nk+1,l_nk+4
         do j= 1-G_haloy, l_nj+G_haloy
            do i= 1-G_halox, l_ni+G_halox
               w5= (VM3%ztht(i,j,k)-VM3%ztht(i,j,l_nk))/(VM3%ztht(i,j,l_nk)-VM3%ztht(i,j,l_nk-1))
               ext_t(i,j,k)= tt0(i,j,l_nk)*(1+w5) -w5*tt0(i,j,l_nk-1)
            !   w5= VM3%ztht(i,j,l_nk)-VM3%ztht(i,j,k)
            !   ext_t(i,j,k) = tt0(i,j,l_nk) + stlo_8 * w5 !??? Schuman-Newel lapse rate
            end do
         end do
      end do
      print*, 'BAC: u,v,w,t'
      i=l_ni/2
      j=l_nj/2+1
      do k=1,l_nk
         write(6,'(i3,5(1pe14.6))') k,ut0(i,j,k),vt0(i,j,k),wt0(i,j,k),tt0(i,j,k)
      end do
      call blocstat (.false.)
      
!!$      do j= ds_j0, ds_jn
!!$         do i= ds_i0, ds_in
!!$            qt0(i,j,l_nk+1)= (qt0(i,j,l_nk) - grav_8*GVM%zmom_8(i,j,l_nk)) &
!!$                            - (grav_8*Cstv_tstr_8)*(GVM%zmom_8(i,j,l_nk+1) &
!!$                            - GVM%zmom_8(i,j,l_nk))/tt0(i,j,l_nk)
!!$         end do
!!$      end do
!     
!     ---------------------------------------------------------------
!
      return
      include 'H3rd_ope.inc'
      end subroutine bacVH
