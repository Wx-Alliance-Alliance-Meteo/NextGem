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
!/@*
      subroutine itf_phy_update3 (F_apply_L)
      use phy_itf, only: phy_get, phymeta, phy_getmeta
      use itf_phy_filter, only: ipf_smooth_tend
      use itf_phy_mem
      use gmm_vt1
      use gmm_pw
      use gmm_phy
      use HORgrid_options
      use dyn_fisl_options
      use adz_options, only: Adz_slt_winds
      use glb_ld
      use cstv
      use lun
      use metric
      use tdpack, only : rgasd_8, grav_8
      use tr3d
      use mem_tstp
      use mem_tracers
      use rmn_gmm
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(in) :: F_apply_L

      logical, parameter :: SMOOTH_EXPLICIT=.false.
      integer, parameter :: SMOOTH_GWD=2
      logical :: source_ps_L

      character(len=GMM_MAXNAMELENGTH) :: trname_S
      integer istat, i,j,k,n, cnt, iteration,iend(3)
      integer :: HLT_start, HLT_end, local_np, trindx(Tr3d_ntr)
      real, dimension(:,:  ), pointer :: ptr2d
      real, dimension(:,:,:), pointer :: ptr3d
      real(kind=REAL64), dimension(:,:,:), pointer :: pm_phy_8
!
!-----------------------------------------------------------------
!
   ! The correction due to sources and sinks of specific humidity
   ! is applied only for the case of dry air conservation
      source_ps_L = (Schm_psadj == 2)
      iend = (/-1,-1,l_nk/)

   ! Make diagnosed winds at the lowest thermodynamic level available for advection
!!$omp single
      if (Adz_slt_winds) then
         ptr2d => pw_uslt(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,gmmk_pw_uslt_S,F_npath='V',F_bpath='V')
         ptr2d => pw_vslt(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
         istat = phy_get(ptr2d,gmmk_pw_vslt_S,F_npath='V',F_bpath='V')
      end if
!!$omp end single
   ! Retrieve a copy of the PW state before the physics
!!$omp do collapse(2)
      do k=1, l_nk
         do j=1,l_nj
            do i=1, l_ni
               pw_uu_plus0(i,j,k) = pw_uu_plus(i,j,k)
               pw_vv_plus0(i,j,k) = pw_vv_plus(i,j,k)
               pw_tt_plus0(i,j,k) = pw_tt_plus(i,j,k)
            end do
         end do
      end do
!!$omp end do      

      if (F_apply_L) then

      if (source_ps_L) then
!!$omp do collapse(2)
         do k=1, l_nk
            do j=1,l_nj
               do i=1, l_ni
                  qw_phy(i,j,k) = 0.
                  qw_dyn(i,j,k) = 0.
               end do
            end do
         end do
!!$omp end do      

         do n= 1, Tr3d_ntr
            trname_S = 'TR/'//trim(Tr3d_name_S(n))//':P'
            if ( (Tr3d_name_S(n)(1:2) == 'HU') .or. &
                 (Schm_wload_L.and.Tr3d_wload(n)) )then
            
!!$omp single
            ptr3d => tdu(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
            i= phy_get (ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                              F_end=iend, F_quiet=.true.)
!!$omp end single copyprivate(i)
              if ( i < 0 ) cycle
!!$omp do collapse(2)
               do k=1, l_nk
                  do j=1+pil_s,l_nj-pil_n
                     do i=1+pil_w,l_ni-pil_e
                        qw_phy(i,j,k)= qw_phy(i,j,k) + tdu(i,j,k)
                        qw_dyn(i,j,k)= qw_dyn(i,j,k) + tracers_P(n)%pntr(i,j,k)
                        tracers_P(n)%pntr(i,j,k)= tdu   (i,j,k)
                     end do
                  end do
               end do
!!$omp end do      
            else
!!$omp single
               ptr3d => tracers_P(n)%pntr(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
               istat = phy_get (ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                                           F_end=iend, F_quiet=.true. )
!!$omp end single
            end if
         end do
      else
!!$omp single
         do k= 1, Tr3d_ntr
            trname_S = 'TR/'//trim(Tr3d_name_S(k))//':P'
            ptr3d => tracers_P(k)%pntr(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
            istat = phy_get ( ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                              F_end=iend, F_quiet=.true. )
            if (Tr3d_name_S(k)(1:2) == 'HU' .and. SMOOTH_EXPLICIT) istat = ipf_smooth_tend(ptr3d,'SQE')
         end do
!!$omp end single
      end if

      ! Apply horizontal filtering on tendencies if requeted

!!$omp single
      ptr3d => pw_uu_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_uu_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      istat = ipf_smooth_tend(ptr3d,'ugwd_td1',SMOOTH_GWD)

      ptr3d => pw_vv_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_vv_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      istat = ipf_smooth_tend(ptr3d,'vgwd_td1',SMOOTH_GWD)

      ptr3d => pw_tt_plus(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
      istat = phy_get(ptr3d,gmmk_pw_tt_plus_s,F_npath='V',F_bpath='D',F_end=iend)
      if (SMOOTH_EXPLICIT) istat = ipf_smooth_tend(ptr3d,'STE')
!!$omp end single

      !Compute moisture sources for dry air conservation (GEM-P)
      !---------------------------------------------------------
      if (source_ps_L) then

         pm_phy_8 => pw_pm_plus_8

         !Obtain Pressure Momentum at TIME P (AFTER DYN)
         !----------------------------------------------
!!$omp do collapse(2)
         do k=1, l_nk+1
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  pm_dyn_8(i,j,k)= pw_pm_plus_8(i,j,k)
                  qt1i(i,j,k)    = qt1(i,j,k)
               end do
            end do
         end do
!!$omp enddo

         do iteration= 1, 4

            !Estimate Pressure Momentum at TIME P (AFTER PHY)
            !------------------------------------------------
            pm_phy_8 => pw_pm_plus_8

            !Estimate source of surface pressure due to fluxes of water:
            !-----------------------------------------------------------------------------------------------------
            !Vertical_Integral [d(p_phy)_k+1] = Vertical_Integral [ d(p_phy)_k q_phy + d(p_dyn) (1-q_dyn) based on
            !-----------------------------------------------------------------------------------------------------
            !d(ps) = Vertical_Integral [ d(qw)/(1-qw_phy)] d(pi) (Claude Girard)
            !-----------------------------------------------------------------------------------------------------
            
!!$omp do
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  p0_8(i,j)= pm_phy_8(i,j,1)
               end do
               do k=1,l_nk
                  do i=1+pil_w,l_ni-pil_e
                     p0_8(i,j) = p0_8(i,j) + (1.0-qw_dyn(i,j,k)) * (pm_dyn_8(i,j,k+1)-pm_dyn_8(i,j,k)) + &
                                                  qw_phy(i,j,k)  * (pm_phy_8(i,j,k+1)-pm_phy_8(i,j,k))
                  end do
               end do
               do i=1+pil_w,l_ni-pil_e
                  qt1(i,j,l_nk+1) = rgasd_8*Cstv_Tstr_8 * &
                            (log(p0_8(i,j))-GVM%lg_pstar_8(i,j,l_nk+1))
                  delq(i,j) = qt1(i,j,l_nk+1) - qt1i(i,j,l_nk+1)
               end do
            end do
!!$omp enddo

!!$omp do collapse(2)
            do k=l_nk,1,-1
               do j=1+pil_s,l_nj-pil_n
                  do i=1+pil_w,l_ni-pil_e
                     qt1(i,j,k) = qt1i(i,j,k) + delq(i,j)
                  end do
               end do
            end do
!!$omp enddo

            !Update Pressure PW_PLUS in accordance with TIME P
            !-------------------------------------------------
            call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,pw_log_pt, &
                            pw_pm_plus_8,pw_p0_plus_8, &
                            l_minx,l_maxx,l_miny,l_maxy,l_nk,1 )


         end do
         call pw_update_GW ()
         
!!$omp single
         if (Lun_out>0) write(Lun_out,*) ''
         if (Lun_out>0) write(Lun_out,*) '--------------------------------------'
         if (Lun_out>0) write(Lun_out,*) 'SOURCE_PS is done for DRY AIR (REAL64)'
         if (Lun_out>0) write(Lun_out,*) '--------------------------------------'
         if (Lun_out>0) write(Lun_out,*) ''
!!$omp end single

      end if

      ! Compute tendencies and reset physical world if requested

      PHY_COUPLING: if (all(phy_cplm == 1.) .and. all(phy_cplt == 1.)) then

         ! Pure split coupling
!!$omp do collapse(2)
         do k=1, l_nk
            do j=1,l_nj
               do i=1, l_ni
                  phy_uu_tend(i,j,k)= 0.
                  phy_vv_tend(i,j,k)= 0.
                  phy_tv_tend(i,j,k)= 0.
               end do
            end do
         end do
!!$omp enddo

         call tt2tvirt (tv, pw_tt_plus, 1,l_ni, 1,l_nj)

!!$omp do collapse(2)
         do k=1, l_nk
            do j=1,l_nj
               do i=1, l_ni
                  phy_tv_tend(i,j,k)= (tv(i,j,k) - tt1(i,j,k))/Cstv_dt_8
               end do
            end do
         end do
!!$omp enddo

      else

         ! Incorporate a component of phy forcing in the rhs of dyn equations
!!$omp do collapse(2)
         do k=1, l_nk
            do j=1,l_nj
               do i=1, l_ni
                  tdu(i,j,k)= pw_uu_plus(i,j,k) - pw_uu_plus0(i,j,k)
                  tdv(i,j,k)= pw_vv_plus(i,j,k) - pw_vv_plus0(i,j,k)
               end do
            end do
         end do
!!$omp enddo

         call HLT_split (1, G_nk, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( tdu(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
         call gem_xch_halo ( tdv(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
         call hwnd_stag_hlt ( phy_uu_tend,phy_vv_tend,tdu,tdv,&
                              l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)
         call tt2tvirt (tv, pw_tt_plus, 1,l_ni, 1,l_nj)

!!$omp do collapse(2)
         do k=1, l_nk
            do j=1,l_nj
               do i=1, l_ni
                  phy_tv_tend(i,j,k) = tv(i,j,k) - tt1(i,j,k)
                  phy_uu_tend(i,j,k) = phy_uu_tend(i,j,k)/Cstv_dt_8
                  phy_vv_tend(i,j,k) = phy_vv_tend(i,j,k)/Cstv_dt_8
                  phy_tv_tend(i,j,k) = phy_tv_tend(i,j,k)/Cstv_dt_8
                  pw_uu_plus (i,j,k) = pw_uu_plus0(i,j,k)+ phy_cplm(i,j)*(pw_uu_plus(i,j,k)-pw_uu_plus0(i,j,k))
                  pw_vv_plus (i,j,k) = pw_vv_plus0(i,j,k)+ phy_cplm(i,j)*(pw_vv_plus(i,j,k)-pw_vv_plus0(i,j,k))
                  pw_tt_plus (i,j,k) = pw_tt_plus0(i,j,k)+ phy_cplt(i,j)*(pw_tt_plus(i,j,k)-pw_tt_plus0(i,j,k))
               end do
            end do
         end do
!!$omp enddo

      endif PHY_COUPLING

   else

!!$omp single
      cnt = 0
      do k= 1, Tr3d_ntr
         if (trim(Tr3d_name_S(k)) == 'HU' .or.                          &
             any(NTR_Tr3d_name_S(1:NTR_Tr3d_ntr)==trim(Tr3d_name_S(k))))&
             cycle
         trname_S = 'TR/'//trim(Tr3d_name_S(k))//':P'
         ptr3d => tracers_P(k)%pntr(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:l_nk)
         if ( phy_get ( ptr3d, trim(trname_S), F_npath='V', F_bpath='D',&
                        F_end=iend, F_quiet=.true. ) < 0 ) cycle
         cnt = cnt + 1
         trindx(cnt) = k
      end do
!!$omp end single copyprivate(cnt)

      if (Grd_yinyang_L) then
         do k= 1, cnt
            call yyg_xchng_hlt (tracers_P(trindx(k))%pntr, l_minx,l_maxx,l_miny,l_maxy,&
                       l_ni,l_nj, G_nk, .true., 'CUBIC', .false.)
         end do
      endif

      if (cnt > 0) then
         call tt2tvirt (tt1(l_minx,l_miny,1), pw_tt_plus, 1,l_ni, 1,l_nj)
         if (Grd_yinyang_L) then
            call yyg_xchng_hlt (tt1(l_minx,l_miny,1), l_minx,l_maxx,l_miny,l_maxy,&
                       l_ni,l_nj, G_nk, .false., 'CUBIC', .false.)
        end if
        call pw_update_T
      end if

   end if
!
!-----------------------------------------------------------------
!
   return
   end subroutine itf_phy_update3
