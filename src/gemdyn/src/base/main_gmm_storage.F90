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

!**s/r main_gmm_storage - Allocate model gmm storage

      subroutine main_gmm_storage()
      use adz_mem
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use glb_ld
      use gmm_geof
      use rmn_gmm
      use HORgrid_options
      use lun
      use metric
      use var_gmm
      use gmm_table
      use gmm_phy
      use mem_tstp
      use psadjust
      use tr3d
      use omp_lib
      implicit none

#include "gmm_gem_flags.hf"
#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      type(gmm_metadata) :: meta, mymeta
      integer(kind=INT64) :: flag_m_f
      integer :: istat,dim,dimH
!
!-------------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,2000)

!     Initialize the time-dependent variables modules
!     -------------------------------------------------
      call heap_paint()

      call set_vt()

      if (Grd_yinyang_L) then
         call yyg_init ()
      else
         call nest_set_mem
      end if

      dimh= (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim = dimh * 6

      call gmm_build_meta1D ( meta, 1,dim,0,0,dim, 0,GMM_NULL_FLAGS )
      istat = gmm_create(gmmk_orography_s ,orography,meta,GMM_FLAG_RSTR+GMM_FLAG_IZER)
                             
      fis0  (l_minx:l_maxx,l_miny:l_maxy) => orography(       1:)
      fis0u (l_minx:l_maxx,l_miny:l_maxy) => orography(  dimh+1:)
      fis0v (l_minx:l_maxx,l_miny:l_maxy) => orography(2*dimh+1:)
      orols (l_minx:l_maxx,l_miny:l_maxy) => orography(3*dimh+1:)
      orolsu(l_minx:l_maxx,l_miny:l_maxy) => orography(4*dimh+1:)
      orolsv(l_minx:l_maxx,l_miny:l_maxy) => orography(5*dimh+1:)
      allocate (zthtu_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                zmomu_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                zthtv_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                zmomv_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )

      call gmm_build_meta3D(meta, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            1,2,0,0,2, 0,GMM_NULL_FLAGS)
      flag_m_f = FLAG_LVL_M
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)

      nullify (topo_low, topo_high, me_full, me_large)
      istat = gmm_create(gmmk_topo_low_s ,topo_low ,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_topo_high_s,topo_high,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_me_full_s, me_full  , meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_me_large_s, me_large, meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_topo_low_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='MEHi'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_topo_high_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='MElo'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_me_full_s   ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='FUME'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_me_large_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='FUML'

      dimH= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      dim = dimH*G_nk
      
!!$      allocate (LT  (5*dim+2*dimH)) ; LT =0.
!!$      Ruu  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>LT(      1:)
!!$      Rvv  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>LT(  dim+1:)
!!$      Rww  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>LT(2*dim+1:)
!!$      Rtt  (l_minx:l_maxx,l_miny:l_maxy,0:l_nk)=>LT(3*dim+1:)
!!$      Rzz  (l_minx:l_maxx,l_miny:l_maxy,0:l_nk)=>LT(4*dim+diMH+1:)
      allocate (Ruu(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (Rvv(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (Rww(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (Rtt(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3))
      allocate (Rzz(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3))
      allocate (Nuu(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (Nvv(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      Nuu= 0. ; Nvv= 0.

      !---for interpolated rhs values at midpoint
      allocate (rhs_mid(6*dim)) ; rhs_mid= 0.
      rhsu_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(      1:)
      rhsv_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(  dim+1:)
      rhst_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(2*dim+1:)
      rhsc_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(3*dim+1:)
      rhsw_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(4*dim+1:)
      rhsf_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(5*dim+1:)
          
      !---for interpolated rhs values at departure
      allocate (rhs_dep(6*dim)) ; rhs_dep= 0.
      rhsu_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(      1:)
      rhsv_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(  dim+1:)
      rhst_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(2*dim+1:)
      rhsc_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(3*dim+1:)
      rhsw_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(4*dim+1:)
      rhsf_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(5*dim+1:)
      
      ! for SW code only
      allocate (sw_rhs(6*dim)) ; sw_rhs=0.
      rhsu (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(      1:)
      rhsv (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(  dim+1:)
      rhst (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(2*dim+1:)
      rhsc (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(3*dim+1:)
      rhsw (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(4*dim+1:)
      rhsf (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_rhs(5*dim+1:)
      
      allocate (orhsu(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsv(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhst(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsc(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsw(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsf(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      orhsu=0.;orhsv=0.;orhst=0.;orhsc=0.;orhsw=0.;orhsf=0.

      !---added true_nlt, true_rhst, needed for back sub---
      allocate (nl_u(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_v(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_t(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_c(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_w(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_f(l_minx:l_maxx,l_miny:l_maxy,l_nk))

      !---for shallow-water budget
      allocate (sw_frc(12*dim))
      sw_f1 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(      1:)
      sw_f2 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(  dim+1:)
      sw_f3 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(2*dim+1:)
      sw_f4 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(3*dim+1:)
      sw_f5 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(4*dim+1:)
      sw_f6 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(5*dim+1:)
      sw_f7 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(6*dim+1:)
      sw_f8 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(7*dim+1:)
      sw_f9 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(8*dim+1:)
      sw_f10 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(9*dim+1:)
      sw_f11 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(10*dim+1:)
      sw_f12 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>sw_frc(11*dim+1:)

      dim = dimH*(G_nk+2)
      allocate (M_Jz(3*dim))
      allocate ( M_Jxozu (l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                 M_Jyozv (l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                 M_iJzq  (l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                 M_logJzu(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                 M_logJzv(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1),&
                 M_logJzq(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )
      M_Jzu(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) => M_jz(      1:)
      M_Jzv(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) => M_jz(  dim+1:)
      M_Jzq(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) => M_jz(2*dim+1:)

      allocate ( GVM%zmom_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%zmom_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%zmom_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%lg_pstar_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )
      allocate ( GVM%zthtlid_8(l_minx:l_maxx,l_miny:l_maxy))

      allocate ( GVM%mc_Jx_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Jy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_iJz_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk), &
               GVM%mc_logJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Ix_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Iy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Iz_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk) )

      allocate ( GVM%mc_css_H_8   (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_alfas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_betas_H_8 (l_minx:l_maxx,l_miny:l_maxy) )

      allocate ( GVM%mc_cst_8   (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_alfat_8 (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_cstp_8  (l_minx:l_maxx,l_miny:l_maxy) )
      GVM%mc_cst_8= 0. ; GVM%mc_alfat_8= 0. ; GVM%mc_cstp_8= 0.

      allocate (psadj_thread_sum(1:2,0:OMP_get_max_threads()-1))
      psadj_thread_sum= 0.

 2000 format( /,'INITIALIZATION OF MAIN GMM VARIABLES S/R MAIN_GMM_STORAGE', &
              /,'====================================================')
!
!     ---------------------------------------------------------------
!
      return
      end

