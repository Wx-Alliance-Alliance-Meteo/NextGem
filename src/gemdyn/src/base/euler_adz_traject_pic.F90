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

      subroutine euler_adz_traject_pic (F_dt_8, itpc)
      use ISO_C_BINDING
      use glb_ld
      use cstv
      use geomh
      use gmm_vt0
      use gmm_vt1
      use ver
      use dyn_fisl_options
      use HORgrid_options
      use adz_options
      use adz_mem
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      integer, intent(IN) :: itpc

      include "tricublin_f90.inc"
      integer :: iter,i,j,k,kk1,nb,k00,dim
      integer :: HLT_np, HLT_start, HLT_end
      integer,dimension(l_ni) :: kk
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8,half_dt_8,double_dt_8,pos
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64), dimension(:,:), pointer :: xchg
      type(C_PTR) :: Cpntr
!
!     ---------------------------------------------------------------
!
      k00=Adz_k0m
      if (Adz_k0>1) k00=1

      half_dt_8= 0.5d0 *Cstv_dt_8

      if(itpc.eq.1) then
!!$omp do collapse(2)
         do k= 1, l_nk
            do j= 1, l_nj
               do i= 1, l_ni
                   ut0(i,j,k)  =  ut1(i,j,k)
                   vt0(i,j,k)  =  vt1(i,j,k)
                  zdt0(i,j,k)  = zdt1(i,j,k)
               end do
            end do
         end do
!!$omp end do

      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut0, vt0,&
                                   l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call yyg_xchng_hlt (zdt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         G_nk, .false., 'CUBIC', .false.)
      end if

      endif

      call adz_prepareWinds (itpc)

!!$omp do
      !0b.Estimate departure points as x^{n-1}=x^{n+1}-dt*V^{n} using
      !current estimate of the velocity
         do k= k00, l_nk
            do j= 1, l_nj
               do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1)-half_dt_8*Adz_uu_arr(i,j,k)*geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1)-half_dt_8*Adz_vv_arr(i,j,k)*geomh_inv_hy_8
                  pos   = Ver_z_8%m(k)  -half_dt_8*Adz_ww_arr(i,j,k)
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                  kk(i) = kk1
               enddo
               do i= 1, l_ni
                  Adz_wpxyz(i,j,k,1) = xm(i)
                  Adz_pxyzm(1,i,j,k) = min(max(xm(i),Adz_iminposx),&
                                                     Adz_imaxposx)
                  Adz_wpxyz(i,j,k,2) = ym(i)
                  Adz_pxyzm(2,i,j,k) = min(max(ym(i),Adz_iminposy),&
                                                     Adz_imaxposy)
                  kk1 = min(l_nk+1,max(0,kk(i)))
                  nb  = max(min(kk1,G_nk-1),1)
                  Adz_wpz(i,j,k) = zm(i)
                  Adz_dpz(i,j,k) = zm(i) !dep = midpoint for euler
                  Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                     *Adz_odelz_m(nb) + dble(nb)
                  Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
               end do
            enddo
          enddo
!!$omp end do

!!$omp do
      !set departure points == midpoints
         do k= k00, l_nk
            do j= 1, l_nj
               do i= 1, l_ni
                  Adz_dpxyz(i,j,k,1) = Adz_wpxyz(i,j,k,1)
                  Adz_pxyzd(1,i,j,k) = Adz_pxyzm(1,i,j,k)
                  Adz_dpxyz(i,j,k,2) = Adz_wpxyz(i,j,k,2)
                  Adz_pxyzd(2,i,j,k) = Adz_pxyzm(2,i,j,k)
                  Adz_dpxyz(i,j,k,3) = Adz_wpxyz(i,j,k,3)
                  Adz_pxyzd(3,i,j,k) = Adz_pxyzm(3,i,j,k) 
               end do
            enddo
          enddo
!!$omp end do
      
!!$omp do
      do k= Adz_k0, l_nk
         Adz_pm   (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzm(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!!$omp enddo nowait

!---for departure point---      
!!$omp do
      do k= Adz_k0, l_nk
         Adz_dep  (:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)=&
         Adz_pxyzd(:,Adz_i0:Adz_in, Adz_j0:Adz_jn, k)
      end do
!!$omp enddo nowait

      dim=3*ubound(Adz_wpxyz,3)
      call HLT_split (1, dim, HLT_np, HLT_start, HLT_end)
      Cpntr= c_loc( Adz_wpxyz(-1,-1,1,1) )
      call C_F_POINTER ( Cpntr, xchg, [(l_ni+4)*(l_nj+4),dim] )
      call gem_xch_halo_8 ( xchg(1,HLT_start),&
                    -1,l_ni+2,-1,l_nj+2, HLT_np,-1)      
      dim=3*ubound(Adz_dpxyz,3)
      call HLT_split (1, dim, HLT_np, HLT_start, HLT_end)
      Cpntr= c_loc( Adz_dpxyz(-1,-1,1,1) )
      call C_F_POINTER ( Cpntr, xchg, [(l_ni+4)*(l_nj+4),dim] )
      call gem_xch_halo_8 ( xchg(1,HLT_start),&
                    -1,l_ni+2,-1,l_nj+2, HLT_np,-1)      

!like adz_interp_traj, but added code to account for the departure point
!at previous time level                 
      call BDF_interp_traj (F_dt_8)

!     ---------------------------------------------------------------
!
      return
      end subroutine euler_adz_traject_pic
