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
      subroutine set_vt ()
      use adz_options
      use gmm_table
      use gmm_contiguous
      use gmm_vt2
      use gmm_vt1
      use gmm_vt0
      use mem_tracers
      use gmm_tracers
      use gmm_geof
      use gmm_pw
      use gmm_smag
      use gmm_phy
      use glb_ld
      use lun
      use tr3d
      use var_gmm
      use, intrinsic :: iso_fortran_env
      implicit none

#include <rmn/msg.h>
#include "gmm_gem_flags.hf"

      type(gmm_metadata) :: mymeta

      integer, parameter :: maxlenght= 32
      character(len=maxlenght) :: nvar
      integer :: i,istat,dim,dim1,dim2,dimH
      integer :: flag_n, flag_r_n
!
!     ---------------------------------------------------------------
!
      flag_n   = GMM_FLAG_IZER
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

      dimh= (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim = dimH * 6*(l_nk + 6)

      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )

      nullify(dynt2,dynt1,dynt0) ; istat=0
      istat = min(gmm_create('DYNT2',dynt2,mymeta,flag_r_n),istat)
      istat = min(gmm_create('DYNT1',dynt1,mymeta,flag_r_n),istat)
      istat = min(gmm_create('DYNT0',dynt0,mymeta,flag_r_n),istat)

      allocate (timlvl2(7),timlvl1(7),timlvl0(7))

      dim = dimH * (l_nk+6)

      timlvl2(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(      1:)
      timlvl2(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(  dim+1:)
      timlvl2(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(2*dim+1:)
      timlvl2(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(3*dim+1:)
      timlvl2(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(4*dim+1:)
      timlvl2(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt2(5*dim+1:)

      istat= gmm_create(gmmk_wt2_s  ,timlvl2(1)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_zdt2_s ,timlvl2(2)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_ut2_s  ,timlvl2(3)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_vt2_s  ,timlvl2(4)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_tt2_s  ,timlvl2(5)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_qt2_s  ,timlvl2(6)%pntr_3d,meta3d_nk3, 0)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_ut2_s  ; GMM_tbl%ara(gmm_cnt)='UU' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_ut2_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_vt2_s  ; GMM_tbl%ara(gmm_cnt)='VV' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_vt2_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_tt2_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_tt2_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_qt2_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MQ' ; GMM_tbl%fst(gmm_cnt)=gmmk_qt2_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_wt2_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_wt2_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_zdt2_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_zdt2_s

      timlvl1(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(      1:)
      timlvl1(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(  dim+1:)
      timlvl1(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(2*dim+1:)
      timlvl1(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(3*dim+1:)
      timlvl1(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(4*dim+1:)
      timlvl1(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt1(5*dim+1:)

      istat= gmm_create(gmmk_wt1_s  ,timlvl1(1)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_zdt1_s ,timlvl1(2)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_ut1_s  ,timlvl1(3)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_vt1_s  ,timlvl1(4)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_tt1_s  ,timlvl1(5)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_qt1_s  ,timlvl1(6)%pntr_3d,meta3d_nk3, 0)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_ut1_s  ; GMM_tbl%ara(gmm_cnt)='UU' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_ut1_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_vt1_s  ; GMM_tbl%ara(gmm_cnt)='VV' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_vt1_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_tt1_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_tt1_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_qt1_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MQ' ; GMM_tbl%fst(gmm_cnt)=gmmk_qt1_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_wt1_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_wt1_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_zdt1_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_zdt1_s

      timlvl0(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(      1:)
      timlvl0(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(  dim+1:)
      timlvl0(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(2*dim+1:)
      timlvl0(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(3*dim+1:)
      timlvl0(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(4*dim+1:)
      timlvl0(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,-2:l_nk+3) => dynt0(5*dim+1:)

      istat= gmm_create(gmmk_wt0_s  ,timlvl0(1)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_zdt0_s ,timlvl0(2)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_ut0_s  ,timlvl0(3)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_vt0_s  ,timlvl0(4)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_tt0_s  ,timlvl0(5)%pntr_3d,meta3d_nk3 ,0)
      istat= gmm_create(gmmk_qt0_s  ,timlvl0(6)%pntr_3d,meta3d_nk3 ,0)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_ut0_s  ; GMM_tbl%ara(gmm_cnt)='UU' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_ut0_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_vt0_s  ; GMM_tbl%ara(gmm_cnt)='VV' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_vt0_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_tt0_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_tt0_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_qt0_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MQ' ; GMM_tbl%fst(gmm_cnt)=gmmk_qt0_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_wt0_s  ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_wt0_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_zdt0_s ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_zdt0_s

      istat= gmm_get (gmmk_ut0_s , ut0)
      istat= gmm_get (gmmk_vt0_s , vt0)
      istat= gmm_get (gmmk_tt0_s , tt0)
      istat= gmm_get (gmmk_wt0_s , wt0)
      istat= gmm_get (gmmk_qt0_s , qt0)
      istat= gmm_get (gmmk_zdt0_s, zdt0)
      istat= gmm_get (gmmk_ut1_s , ut1)
      istat= gmm_get (gmmk_vt1_s , vt1)
      istat= gmm_get (gmmk_tt1_s , tt1)
      istat= gmm_get (gmmk_wt1_s , wt1)
      istat= gmm_get (gmmk_qt1_s , qt1)
      istat= gmm_get (gmmk_zdt1_s, zdt1)
      istat= gmm_get (gmmk_ut2_s , ut2)
      istat= gmm_get (gmmk_vt2_s , vt2)
      istat= gmm_get (gmmk_tt2_s , tt2)
      istat= gmm_get (gmmk_wt2_s , wt2)
      istat= gmm_get (gmmk_qt2_s , qt2)
      istat= gmm_get (gmmk_zdt2_s, zdt2)

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * (l_nk+3) * Tr3d_ntr
      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )

      istat = GMM_OK
      nullify(trt0,trt1,trt2,trdf)
      istat= gmm_create('TRACERS:t0',trt0,mymeta,flag_r_n)
      istat= gmm_create('TRACERS:t1',trt1,mymeta,flag_r_n)
      istat= gmm_create('TRACERS:t2',trt2,mymeta,flag_r_n)
      istat= gmm_create('TRACERS:df',trdf,mymeta,flag_r_n)

      allocate (tracers_P(Tr3d_ntr), tracers_M(Tr3d_ntr), &
                tracers_t2(Tr3d_ntr))
      allocate (sumq_8(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (air_mass(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate (w_tr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      
      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * (l_nk)
      do i=1,Tr3d_ntr
         nullify(tracers_P(i)%pntr)
         tracers_P(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt1((i-1)*dim+1:)
         tracers_M(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((i-1)*dim+1:)
         tracers_t2(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt2((i-1)*dim+1:)
!     necessary to create those GMM tracers variables for phy_input
         nvar(1:maxlenght)='' ; nvar= 'TR/'//trim(Tr3d_name_S(i))//':M'
         istat= gmm_create(nvar,tracers_M(i)%pntr,meta3d_nk,0)
         nvar(1:maxlenght)='' ; nvar= 'TR/'//trim(Tr3d_name_S(i))//':P'
         istat= gmm_create(nvar,tracers_P(i)%pntr,meta3d_nk,0)
         nvar(1:maxlenght)='' ; nvar= 'TR/'//trim(Tr3d_name_S(i))//':t2'
         istat= min(gmm_create(nvar,tracers_t2(i)%pntr,meta3d_nk,0),istat)
      end do

      !Allocation if Bermejo-Conde LAM ZFL
      !-----------------------------------
      if (adz_BC_LAM_flux==2) then

         istat = GMM_OK
         nullify(trtb)
         istat= gmm_create('TRACERS:B',trtb,mymeta,flag_r_n)
         allocate (tracers_B(Tr3d_ntr))

         dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * (l_nk)
         do i=1,Tr3d_ntr
            nullify(tracers_B(i)%pntr)
            tracers_B(i)%pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trtb((i-1)*dim+1:)
            istat= gmm_create('TR/'//trim(Tr3d_name_S(i))//':B',tracers_B(i)%pntr,meta3d_nk,0)
         end do

      end if

      istat = GMM_OK

      istat= gmm_create(gmmk_airm1_s,airm1,meta3d_nk,flag_r_n)
      istat= gmm_create(gmmk_airm0_s,airm0,meta3d_nk,flag_r_n)
      istat= gmm_create(gmmk_pkps_s ,pkps ,meta3d_nk,flag_r_n)

      istat= gmm_create(gmmk_dgzm_s   ,dgzm  ,meta3d_nk  ,0)
      istat= gmm_create(gmmk_dgzt_s   ,dgzt  ,meta3d_nk  ,0)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_dgzm_s; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)=gmmk_dgzm_s
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)=gmmk_dgzt_s; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='TH' ; GMM_tbl%fst(gmm_cnt)=gmmk_dgzt_s

      dimh= (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim = dimH * (7*(l_nk+1) + 2)

      call gmm_build_meta1D ( mymeta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )

      nullify(pwPLUS, pwMOINS)
      istat= gmm_create('pwPLUS' ,pwPLUS ,mymeta,flag_r_n)
      istat= gmm_create('pwMOINS',pwMOINS,mymeta,flag_r_n)

      allocate (pw_pnt_P(7),pw_pnt_M(7))

      dim1 = dimH * (l_nk+1)
      dim2 = 7*dim1

      pw_pnt_P(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwPLUS(      1:)
      pw_pnt_P(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwPLUS(  dim1+1:)
      pw_pnt_P(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwPLUS(2*dim1+1:)
      pw_pnt_P(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwPLUS(3*dim1+1:)
      pw_pnt_P(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> pwPLUS(4*dim1+1:)
      pw_pnt_P(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> pwPLUS(5*dim1+1:)
      pw_pnt_P(7)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwPLUS(6*dim1+1:)
      pw_pnt_P(1)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => pwPLUS(     dim2+1:)
      pw_pnt_P(2)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => pwPLUS(dimH+dim2+1:)

      pw_pnt_M(1)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwMOINS(      1:)
      pw_pnt_M(2)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwMOINS(  dim1+1:)
      pw_pnt_M(3)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwMOINS(2*dim1+1:)
      pw_pnt_M(4)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwMOINS(3*dim1+1:)
      pw_pnt_M(5)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> pwMOINS(4*dim1+1:)
      pw_pnt_M(6)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk+1)=> pwMOINS(5*dim1+1:)
      pw_pnt_M(7)%pntr_3d(l_minx:l_maxx,l_miny:l_maxy,1:l_nk)  => pwMOINS(6*dim1+1:)
      pw_pnt_M(1)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => pwMOINS(     dim2+1:)
      pw_pnt_M(2)%pntr_2d(l_minx:l_maxx,l_miny:l_maxy) => pwMOINS(dimH+dim2+1:)
      
      istat= gmm_create(gmmk_pw_uu_plus_s  ,pw_pnt_P(1)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_vv_plus_s  ,pw_pnt_P(2)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_tt_plus_s  ,pw_pnt_P(3)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_wz_plus_s  ,pw_pnt_P(4)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_pt_plus_s  ,pw_pnt_P(5)%pntr_3d, meta3d_nk1,0)
      istat= gmm_create(gmmk_pw_pm_plus_s  ,pw_pnt_P(6)%pntr_3d, meta3d_nk1,0)
      istat= gmm_create(gmmk_pw_gz_plus_s  ,pw_pnt_P(7)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_p0_plus_s  ,pw_pnt_P(1)%pntr_2d, meta2d    ,0)
      istat= gmm_create(gmmk_pw_me_plus_s  ,pw_pnt_P(2)%pntr_2d, meta2d    ,0)
      istat= gmm_get (gmmk_pw_uu_plus_s , pw_uu_plus)
      istat= gmm_get (gmmk_pw_vv_plus_s , pw_vv_plus)
      istat= gmm_get (gmmk_pw_tt_plus_s , pw_tt_plus)
      istat= gmm_get (gmmk_pw_wz_plus_s , pw_wz_plus)
      istat= gmm_get (gmmk_pw_pt_plus_s , pw_pt_plus)
      istat= gmm_get (gmmk_pw_pm_plus_s , pw_pm_plus)
      istat= gmm_get (gmmk_pw_gz_plus_s , pw_gz_plus)
      istat= gmm_get (gmmk_pw_p0_plus_s , pw_p0_plus)
      istat= gmm_get (gmmk_pw_me_plus_s , pw_me_plus)

      istat= gmm_create(gmmk_pw_uu_moins_s  ,pw_pnt_M(1)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_vv_moins_s  ,pw_pnt_M(2)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_tt_moins_s  ,pw_pnt_M(3)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_pt_moins_s  ,pw_pnt_M(5)%pntr_3d, meta3d_nk1,0)
      istat= gmm_create(gmmk_pw_pm_moins_s  ,pw_pnt_M(6)%pntr_3d, meta3d_nk1,0)
      istat= gmm_create(gmmk_pw_gz_moins_s  ,pw_pnt_M(7)%pntr_3d, meta3d_nk ,0)
      istat= gmm_create(gmmk_pw_p0_moins_s  ,pw_pnt_M(1)%pntr_2d, meta2d    ,0)
      istat= gmm_create(gmmk_pw_me_moins_s  ,pw_pnt_M(2)%pntr_2d, meta2d    ,0)
      istat= gmm_get (gmmk_pw_uu_moins_s , pw_uu_moins)
      istat= gmm_get (gmmk_pw_vv_moins_s , pw_vv_moins)
      istat= gmm_get (gmmk_pw_tt_moins_s , pw_tt_moins)
      istat= gmm_get (gmmk_pw_pt_moins_s , pw_pt_moins)
      istat= gmm_get (gmmk_pw_pm_moins_s , pw_pm_moins)
      istat= gmm_get (gmmk_pw_gz_moins_s , pw_gz_moins)
      istat= gmm_get (gmmk_pw_p0_moins_s , pw_p0_moins)
      istat= gmm_get (gmmk_pw_me_moins_s , pw_me_moins)
      
      nullify(pw_uu_copy ,pw_vv_copy, pw_log_pm, pw_log_pt, smag)
      nullify(pw_pm_plus_8,pw_p0_plus_8,pw_pm_moins_8,pw_p0_moins_8)

      istat= gmm_create(gmmk_pw_log_pm_s    ,pw_log_pm    ,meta3d_nk1 ,flag_r_n)
      istat= gmm_create(gmmk_pw_log_pt_s    ,pw_log_pt    ,meta3d_nk1 ,flag_r_n)
      istat= gmm_create(gmmk_pw_pm_plus_8_s ,pw_pm_plus_8 ,meta3d_nk1 ,flag_r_n)
      istat= gmm_create(gmmk_pw_p0_plus_8_s ,pw_p0_plus_8 ,meta2d     ,flag_r_n)
      istat= gmm_create(gmmk_pw_pm_moins_8_s,pw_pm_moins_8,meta3d_nk1 ,flag_r_n)
      istat= gmm_create(gmmk_pw_p0_moins_8_s,pw_p0_moins_8,meta2d     ,flag_r_n)
      istat= gmm_create(gmmk_pw_uslt_s      ,pw_uslt      ,meta2d     ,flag_r_n)
      istat= gmm_create(gmmk_pw_vslt_s      ,pw_vslt      ,meta2d     ,flag_r_n)
      istat= gmm_create(gmmk_pw_uu_copy_s   ,pw_uu_copy   ,meta3d_nk  ,flag_r_n)
      istat= gmm_create(gmmk_pw_vv_copy_s   ,pw_vv_copy   ,meta3d_nk  ,flag_r_n)
      istat= gmm_create(gmmk_smag_s         ,smag         ,meta3d_nk  ,flag_n  )

      istat= gmm_create(gmmk_phy_cplm_s   , phy_cplm   , meta2d   , flag_r_n)
      istat= gmm_create(gmmk_phy_cplt_s   , phy_cplt   , meta2d   , flag_r_n)
      istat= gmm_create(gmmk_phy_uu_tend_s, phy_uu_tend, meta3d_nk, flag_r_n)
      istat= gmm_create(gmmk_phy_vv_tend_s, phy_vv_tend, meta3d_nk, flag_r_n)
      istat= gmm_create(gmmk_phy_tv_tend_s, phy_tv_tend, meta3d_nk, flag_r_n)

      call canonical_cases ("SET_VT")
!
!     ---------------------------------------------------------------
!
      return
      end subroutine set_vt
