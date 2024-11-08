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

      subroutine SW_adz_traject (F_dt_8, itpc)
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
      use ptopo
      use glb_pil
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      integer, intent(IN) :: itpc

      include "tricublin_f90.inc"
      include 'mpif.h'

      integer :: iter,i,j,k,kk1,nb,k00,dim,ierr
      integer :: HLT_np, HLT_start, HLT_end
      integer,dimension(l_ni) :: kk
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8,half_dt_8,double_dt_8,pos
      real(kind=REAL64), dimension(l_ni) :: xm,ym,zm
      real(kind=REAL64), dimension(:,:), pointer :: xchg
      !real(kind=REAL64), dimension(1:l_ni,1:l_nj,1:l_nk) :: Adz_wpz, Adz_dpz
      real(kind=REAL64), parameter :: bdf1=0.75d0, bdf2=0.25d0, bdf3=-3.d0, bdf4=4.d0
      type(C_PTR) :: Cpntr

      real :: rate
      real(kind=REAL64) :: N,err2nrm
      real(kind=REAL64), dimension(l_nk) :: total_sum,err,glb
      real(kind=REAL64), dimension(Adz_i0:Adz_in) :: lenght
      real(kind=REAL64), dimension(500), save :: sucsv_err
      logical :: print_conv=.true.

!
!     ---------------------------------------------------------------
!

      if (Schm_advec == 0) then ! no advection

       half_dt_8 = 0.d0
       double_dt_8 = 0.d0

      else

       half_dt_8 = 0.5*F_dt_8
       double_dt_8 = 2.d0*F_dt_8

      end if

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

      N = ((G_ni-Glb_pil_e)-(1+Glb_pil_w)+1)*((G_nj-Glb_pil_n)-(1+Glb_pil_s)+1)

      k00=Adz_k0m
      if (Adz_k0>1) k00=1

      if(itpc.eq.1) then

!!$omp do
      !0a.Interpolate velocity V(t) at current time to estimated midpoints x^{n}: V^{n}
      do k= k00, l_nk
       call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
      end do
!!$omp end do

      endif !itpc=1

      do  iter = 1, Schm_itSL

         total_sum= 0.d0 ; err = -huge(1.) 
!!$omp do
      !1.Solve x^{n-1} = x^{n+1} - 0.5*dt*(V^{n} + V^{+})
         do k= k00, l_nk
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - half_dt_8*(Adz_uvw_dep(1,i,j,k) +  Adz_uu_arr(i,j,k))*geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - half_dt_8*(Adz_uvw_dep(2,i,j,k) +  Adz_vv_arr(i,j,k))*geomh_inv_hy_8
                  pos   = Ver_z_8%m(k)   - half_dt_8*(Adz_uvw_dep(3,i,j,k) +  Adz_ww_arr(i,j,k))
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                  kk(i) = kk1
               end do
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
                  Adz_wpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                     *Adz_odelz_m(nb) + dble(nb)
                  Adz_pxyzm(3,i,j,k) = Adz_wpxyz(i,j,k,3)
               end do
            end do
            do j= Adz_j0, Adz_jn
            do i= Adz_i0, Adz_in
               lenght(i)= sqrt((Adz_pxyzm(1,i,j,k)-Adz_pm(1,i,j,k))**2. &
                             + (Adz_pxyzm(2,i,j,k)-Adz_pm(2,i,j,k))**2. &
                             + (Adz_pxyzm(3,i,j,k)-Adz_pm(3,i,j,k))**2.)
               Adz_pm(:,i,j,k) = Adz_pxyzm(:,i,j,k)
            end do !end i
            total_sum(k)= total_sum(k) + sum(lenght)
            err(k) = max(err(k),maxval(lenght))
            end do !end j
         enddo
!!$omp enddo

!!$omp do
      !0a.Interpolate velocity V(t) at current time to estimated midpoints x^{n}: V^{n}
      do k= k00, l_nk
       call tricublin_zyx3_n ( Adz_uvw_dep(1,1,1,k),Adz_uvw_d(1,1,1,1), &
                                Adz_pxyzm(1:3,:,:,k),Adz_cpntr_q,Adz_2dnh)
      end do
!!$omp end do

!!$omp do
      !2.Solve x^{n-1} = x^{n+1} - 2*dt*(V^{n})
         do k= k00, l_nk
             do j= 1, l_nj
                do i= 1, l_ni
                  xm(i) = dble(i+l_i0-1) - double_dt_8*(Adz_uvw_dep(1,i,j,k))*geomh_inv_hx_8
                  ym(i) = dble(j+l_j0-1) - double_dt_8*(Adz_uvw_dep(2,i,j,k))*geomh_inv_hy_8
                  pos   = Ver_z_8%m(k)   - double_dt_8*(Adz_uvw_dep(3,i,j,k))
                  zm(i) = min(max(pos,Ver_zmin_8),Ver_zmax_8)

                  kk1 = (zm(i) - ver_z_8%m(0)  ) * adz_ovdzm_8 + 1.d0
                  kk1 = min(max(1,kk1),ubound(adz_search_m,1))
                  kk1 = adz_search_m(kk1)
                  if ( sig * zm(i) > sig* ver_z_8%m(min(kk1+1,l_nk+1)) ) kk1= kk1 + 1
                  if ( sig * zm(i) < sig* ver_z_8%m(kk1              ) ) kk1= kk1 - 1
                  kk(i) = kk1
               end do
               do i= 1, l_ni
                  Adz_dpxyz(i,j,k,1) = xm(i)
                  Adz_pxyzd(1,i,j,k) = min(max(xm(i),Adz_iminposx),&
                                                     Adz_imaxposx)
                  Adz_dpxyz(i,j,k,2) = ym(i)
                  Adz_pxyzd(2,i,j,k) = min(max(ym(i),Adz_iminposy),&
                                                     Adz_imaxposy)
                  kk1 = min(l_nk+1,max(0,kk(i)))
                  nb  = max(min(kk1,G_nk-1),1)
                  Adz_dpz(i,j,k) =zm(i)
                  Adz_dpxyz(i,j,k,3) = (zm(i)-ver_z_8%m(nb))&
                                     *Adz_odelz_m(nb) + dble(nb)
                  Adz_pxyzd(3,i,j,k) = Adz_dpxyz(i,j,k,3)
               end do
            end do
         enddo
!!$omp enddo

         !---check for convergence on midpoint---
         if (Schm_tolSL>0.) then
            call MPI_allreduce (err      ,glb,l_nk,MPI_DOUBLE_PRECISION,&
                                MPI_MAX,COMM_MULTIGRID,ierr)
            err= glb
            call MPI_allreduce (total_sum,glb,l_nk,MPI_DOUBLE_PRECISION,&
                                MPI_SUM,COMM_MULTIGRID,ierr)
            err2nrm = maxval(glb)/N
            sucsv_err(iter)= err2nrm
            rate=1.d0
            if (iter > 1) rate= (sucsv_err(iter-1)-sucsv_err(iter))/&
                                 sucsv_err(iter-1)
            if (print_conv) write(Lun_out,1001) iter,maxval(err),err2nrm,rate
            if ((err2nrm < Schm_tolSL).or.(rate<Schm_rateSL)) exit
         endif

      end do !Adz_niter
      
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

      call adz_int_traj (Adz_wpxyz,Adz_wpz,Adz_pmu,Adz_pmv,Adz_pt ,F_dt_8)
      call adz_int_traj (Adz_dpxyz,Adz_dpz,Adz_pdu,Adz_pdv,Adz_pdt,F_dt_8)

 1001 format (3x,"TRAJ convergence",i3,": DIFF_MAX,L_2,rate: ",2(1pe11.4),1pe10.2)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_adz_traject
