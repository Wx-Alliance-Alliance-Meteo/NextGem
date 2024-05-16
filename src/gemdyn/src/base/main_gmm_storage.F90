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

!     Initialize digital filter variables modules
!     --------------------------------------------

      nullify (fis0, fis0u, fis0v, orols, orolsu, orolsv)
      nullify (topo_low, topo_high, me_full, me_large)
      
      istat = gmm_create(gmmk_fis0_s ,fis0  ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create('FIS0U'     ,fis0u ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create('FIS0V'     ,fis0v ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_orols_s,orols ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create('OROLSU'    ,orolsu,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create('OROLSV'    ,orolsv,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)

      call gmm_build_meta3D(meta, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            1,2,0,0,2, 0,GMM_NULL_FLAGS)
      flag_m_f = FLAG_LVL_M
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)

      istat = gmm_create(gmmk_topo_low_s ,topo_low ,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_topo_high_s,topo_high,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)

      istat = gmm_create(gmmk_me_full_s, me_full  , meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_me_large_s, me_large, meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_fis0_s      ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='FIS0'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_orols_s     ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='LSOR'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_topo_low_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='MEHi'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_topo_high_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='MElo'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_me_full_s   ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='FUME'
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_me_large_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='SF' ; GMM_tbl%fst(gmm_cnt)='FUML'

      allocate (rhs_zero(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate (rhsb(l_minx:l_maxx,l_miny:l_maxy))
      rhsb = 0.

      dimH= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
      dim = dimH*G_nk
      allocate (rhs(6*dim))
      rhsu (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(      1:)
      rhsv (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(  dim+1:)
      rhst (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(2*dim+1:)
      rhsc (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(3*dim+1:)
      rhsw (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(4*dim+1:)
      rhsf (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(5*dim+1:)

      rhsu=0.;rhsv=0.;rhst=0.;rhsc=0.;rhsw=0.;rhsf=0.;rhs_zero=0.

      !---for interpolated rhs values at midpoint
      allocate (rhs_mid(6*dim))
      rhsu_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(      1:)
      rhsv_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(  dim+1:)
      rhst_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(2*dim+1:)
      rhsc_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(3*dim+1:)
      rhsw_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(4*dim+1:)
      rhsf_mid   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_mid(5*dim+1:)

      rhsu_mid=0.;rhsv_mid=0.;rhst_mid=0.
      rhsc_mid=0.;rhsw_mid=0.;rhsf_mid=0.
      
      allocate (var_init(6*dim+2*dimH))
      wti  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>var_init(      1:)
      zdti (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>var_init(  dim+1:)
      uti  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>var_init(2*dim+1:)
      vti  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>var_init(3*dim+1:)
      tti  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>var_init(4*dim+1:)
      qti  (l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=>var_init(5*dim+1:)
      sti  (l_minx:l_maxx,l_miny:l_maxy)=>var_init(6*dim+dimH+1:)

      !---for interpolated rhs values at departure
      allocate (rhs_dep(6*dim))
      rhsu_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(      1:)
      rhsv_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(  dim+1:)
      rhst_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(2*dim+1:)
      rhsc_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(3*dim+1:)
      rhsw_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(4*dim+1:)
      rhsf_dep   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs_dep(5*dim+1:)

      rhsu_dep=0.;rhsv_dep=0.;rhst_dep=0.
      rhsc_dep=0.;rhsw_dep=0.;rhsf_dep=0.

      !---for interpolated nonlinear values at midpoint
      allocate (nl_mid(5*dim))
      nlu_mid (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_mid(      1:)
      nlv_mid (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_mid(  dim+1:)
      nlt_mid (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_mid(2*dim+1:)
      nlw_mid (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_mid(3*dim+1:)
      nlq_mid (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_mid(4*dim+1:)

      nlu_mid=0.;nlv_mid=0.;nlt_mid=0.;nlw_mid=0.;nlq_mid=0.

      !---for interpolated nonlinear values at departure
      allocate (nl_dep(5*dim))
      nlu_dep (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_dep(      1:)
      nlv_dep (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_dep(  dim+1:)
      nlt_dep (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_dep(2*dim+1:)
      nlw_dep (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_dep(3*dim+1:)
      nlq_dep (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>nl_dep(4*dim+1:)

      nlu_dep=0.;nlv_dep=0.;nlt_dep=0.;nlw_dep=0.;nlq_dep=0.
      
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
                nl_f(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_b(l_minx:l_maxx,l_miny:l_maxy)     ,&
                true_nlt(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                true_rhst(l_minx:l_maxx,l_miny:l_maxy,l_nk) )
      nl_b = 0.

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

      allocate ( GVM%zmom_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%zmom_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%zmom_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%ztht_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                 GVM%lg_pstar_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )

      allocate ( GVM%mc_Jx_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Jy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_iJz_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk), &
               GVM%mc_logJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Ix_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Iy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                 GVM%mc_Iz_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk) )

      allocate ( GVM%mc_css_H_8   (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_alfas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_betas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                 GVM%mc_cssp_H_8  (l_minx:l_maxx,l_miny:l_maxy) )

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

