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

      subroutine adz_set()
      use adz_mem
      use adz_options
      use ctrl
      use dyn_fisl_options
      use gem_options
      use glb_pil
      use HORgrid_options
      use geomh
      use lam_options
      use init_options
      use step_options
      use rmn_gmm
      use gmm_pw
      use mem_tstp
      use ptopo
      use rstr
      use tr3d
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

#include "gmm_gem_flags.hf"
#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      include 'mpif.h'
      include "tricublin_f90.inc"
      
      integer i,j,k, n, k0, BCS_BASE, pnz, ext, halom, dim, ub, ierr
      real(kind=REAL64) :: ra,rb,rc,rd,rx,re
      real(kind=REAL64), parameter :: EPS_8= 2.D-5
      real(kind=REAL64), dimension(:), allocatable :: lat
      real :: verysmall
      real(kind=REAL128) :: smallest_dz,posz
      real   , dimension(:,:), contiguous, pointer :: temp_rptr
      integer, dimension(:,:), contiguous, pointer :: temp_iptr

      type(gmm_metadata) :: meta, mymeta
      integer(kind=INT64) :: flag_m_f
      integer :: kp1, flag_r_n, istat,istatu,istatv
      integer :: DISP_UNIT, INFO
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      type(C_PTR), save :: BASEPTR_q,BASEPTR_u,BASEPTR_v,BASEPTR_t,&
                           BASEtrajq,BASEtraju,BASEtrajv,BASEtrajt,&
                           BASEcor, BASECOFFS, BASEPTR_list, BASEPTR_pos
!
!---------------------------------------------------------------------
!
      ADZ_OD_L= .false.

      if (ADZ_OD_L) then
         Adz_halox = G_halox
         Adz_haloy = G_haloy
      else
         Adz_maxcfl= max(1,Grd_maxcfl)
         Adz_halox = Adz_maxcfl + 1
         Adz_haloy = Adz_halox

         Adz_lminx = 1    - Adz_halox
         Adz_lmaxx = l_ni + Adz_halox
         Adz_lminy = 1    - Adz_haloy
         Adz_lmaxy = l_nj + Adz_haloy

         Adz_nit = Adz_lmaxx - Adz_lminx + 1
         Adz_njt = Adz_lmaxy - Adz_lminy + 1
         Adz_nij = Adz_nit * Adz_njt
      endif

      Adz_2dnh = l_ni * l_nj
      Adz_3dnh = Adz_2dnh * l_nk

      Adz_ioff = l_i0 - 2
      Adz_joff = l_j0 - 2

      ext=3
      if (Grd_yinyang_L) ext=4

      Adz_i0=   1  + pil_w - (ext+1)*west
      Adz_in= l_ni - pil_e +  ext   *east
      Adz_j0=   1  + pil_s - (ext+1)*south
      Adz_jn= l_nj - pil_n +  ext   *north

      Adz_i0u = Adz_i0 ;  Adz_inu = Adz_in
      Adz_j0v = Adz_j0 ;  Adz_jnv = Adz_jn
      if (.not. Grd_yinyang_L) then
         Adz_i0u= 1 + pil_w     ; Adz_j0v= 1 + pil_s
         Adz_inu= l_niu - pil_e ; Adz_jnv= l_njv - pil_n
      end if

      Adz_k0 = Lam_gbpil_t+1
      Adz_k0t= max(1,Lam_gbpil_t)
      Adz_k0m= max(Adz_k0t-2,1)
      if (adz_BC_LAM_flux/=0) Adz_k0m=1

      Adz_num_u= (Adz_inu-Adz_i0u+1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0 +1)
      Adz_num_v= (Adz_in -Adz_i0 +1)*(Adz_jnv-Adz_j0v+1)*(l_nk-Adz_k0 +1)
      Adz_num_q= (Adz_in -Adz_i0 +1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0 +1)
      Adz_num_t= (Adz_in -Adz_i0 +1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0t+1)

      Adz_kkmax= l_nk-2

      if (ADZ_OD_L) then
         Adz_cpntr_q= vsearch_setup_plus(Ver_z_8%m(1:G_nk), G_nk,&
                    l_maxx-l_minx+1,l_maxy-l_miny+1,1-l_i0,1-l_j0)
         Adz_cpntr_t= vsearch_setup_plus(Ver_z_8%t(1:G_nk), G_nk,&
                    l_maxx-l_minx+1,l_maxy-l_miny+1,1-l_i0,1-l_j0)
         Adz_iminposx = l_i0+1    - Adz_halox   + EPS_8
         Adz_imaxposx = l_i0+l_ni + Adz_halox-2 - EPS_8
         Adz_iminposy = l_j0+1    - Adz_haloy   + EPS_8
         Adz_imaxposy = l_j0+l_nj + Adz_haloy-2 - EPS_8
      else
         Adz_cpntr_q= vsearch_setup_plus(Ver_z_8%m(1:G_nk), G_nk,&
                                 Adz_nit,Adz_njt,1-l_i0,1-l_j0)
         Adz_cpntr_t= vsearch_setup_plus(Ver_z_8%t(1:G_nk), G_nk,&
                                 Adz_nit,Adz_njt,1-l_i0,1-l_j0)
         Adz_iminposx = l_i0+adz_lminx   + EPS_8
         Adz_imaxposx = l_i0+adz_lmaxx-2 - EPS_8
         Adz_iminposy = l_j0+adz_lminy   + EPS_8
         Adz_imaxposy = l_j0+adz_lmaxy-2 - EPS_8
      endif
      BCS_BASE= 4
      if (Grd_yinyang_L) BCS_BASE = 3

      allocate (Adz_gindx_alongX(0:Ptopo_npex), Adz_gindx_alongY(0:Ptopo_npey))
      do i=0, Ptopo_npex-2
         Adz_gindx_alongX(i) = Ptopo_gindx(2,Ptopo_colrow(Ptopo_couleur,i,0)+1)
      end do
      Adz_gindx_alongX(Ptopo_npex-1) = G_ni - BCS_BASE - 1
      Adz_gindx_alongX(Ptopo_npex  ) = 2*G_ni
      do i=0, Ptopo_npey-2
         Adz_gindx_alongY(i) = Ptopo_gindx(4,Ptopo_colrow(Ptopo_couleur,0,i)+1)
      end do
      Adz_gindx_alongY(Ptopo_npey-1) = G_nj - BCS_BASE - 1
      Adz_gindx_alongY(Ptopo_npey  ) = 2*G_nj

! ceci est a revoir dans le contexte ou on peut placer des donnees dans le halo externe
! et calculer les RHS sur potentiellement plus large et ainsi repousser la frontiere
! de troncature

      if (l_west ) Adz_iminposx = l_i0+     BCS_BASE   + EPS_8
      if (l_east ) Adz_imaxposx = l_i0+l_ni-BCS_BASE-1 - EPS_8
      if (l_south) Adz_iminposy = l_j0+     BCS_BASE   + EPS_8
      if (l_north) Adz_imaxposy = l_j0+l_nj-BCS_BASE-1 - EPS_8
      Adz_glbminpos= dble(BCS_BASE+1) + EPS_8

      Adz_yyminposx = 1+ BCS_BASE + EPS_8
      Adz_yymaxposx = 1+ G_ni - BCS_BASE -1 - EPS_8
      Adz_yyminposy = 1+ BCS_BASE + EPS_8
      Adz_yymaxposy = 1+ G_nj - BCS_BASE -1 - EPS_8

      Adz_i0b =    1 + BCS_BASE*west
      Adz_inb = l_ni - BCS_BASE*east
      Adz_j0b =    1 + BCS_BASE*south
      Adz_jnb = l_nj - BCS_BASE*north

      Adz_num_b = (Adz_inb -Adz_i0b +1)*(Adz_jnb -Adz_j0b +1)*l_nk

      allocate (  Adz_delz_m(0:l_nk),  Adz_delz_t(0:l_nk), &
                 Adz_odelz_m(0:l_nk), Adz_odelz_t(0:l_nk) )

! Vert coord: sig= 1.d0 for P coord, sig= -1.d0 for H coord
      sig=       (Ver_z_8%m(l_nk)-Ver_z_8%m(1)) / &
            ( abs(Ver_z_8%m(l_nk)-Ver_z_8%m(1)) )

      Adz_delz_m= 9.e33 ; Adz_delz_t= 9.e33
      verysmall= tiny(verysmall)*1.e5
      do k = 0,l_nk
         kp1 = k+1
         if (abs(Ver_z_8%m(kp1)-Ver_z_8%m(k))>verysmall) &
         Adz_delz_m(k) = Ver_z_8%m(kp1)-Ver_z_8%m(k)
         if (abs(Ver_z_8%t(kp1)-Ver_z_8%t(k))>verysmall) &
         Adz_delz_t(k) = Ver_z_8%t(kp1)-Ver_z_8%t(k)
         Adz_odelz_m(k) = 1.0d0 / Adz_delz_m(k)
         Adz_odelz_t(k) = 1.0d0 / Adz_delz_t(k)
      end do

      smallest_dz= minval(sig*Adz_delz_m(1:l_nk))
      adz_ovdzm_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+ sig*(Ver_z_8%m(l_nk+1)-Ver_z_8%m(0))*adz_ovdzm_8)
      adz_ovdzm_8 = sig/smallest_dz
      allocate ( Adz_search_m(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = sig*Ver_z_8%m(0) + dble(k-1) * smallest_dz
         if (posz > sig*Ver_z_8%m(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_m(k) = k0
      end do

      smallest_dz= minval(sig*Adz_delz_t(1:l_nk))
      adz_ovdzt_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+ sig*(Ver_z_8%t(l_nk+1)-Ver_z_8%t(0))*adz_ovdzt_8)
      adz_ovdzt_8 = sig/smallest_dz
      allocate ( Adz_search_t(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = sig*Ver_z_8%t(0) + dble(k-1) * smallest_dz
         if (posz > sig*Ver_z_8%t(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_t(k) = k0
      end do

      if (.not.ADZ_OD_L) then

      allocate ( Adz_cy_8(Adz_lminy:Adz_lmaxy) )
      halom= max(Adz_haloy,G_haloy)
      allocate ( lat(1-halom:G_nj+halom+1))
      lat(1-G_haloy:G_nj+G_haloy+1)= G_yg_8(1-G_haloy:G_nj+G_haloy+1)
      do j= lbound(lat,dim=1), -G_haloy
         lat(j)= G_yg_8(1-G_haloy) + (j-1+G_haloy) * geomh_hy_8
      end do
      do j = G_nj+G_haloy+2, ubound(lat,dim=1)
         lat(j)= G_yg_8(G_nj+G_haloy) + (j-G_nj-G_haloy) * geomh_hy_8
      end do
      do j = Adz_lminy, Adz_lmaxy
         Adz_cy_8(j) = 1.d0 / cos(lat(l_j0+j-1))
      end do
      deallocate(lat)

      else

      allocate ( Adz_cy_8(l_miny:l_maxy) )
      halom= G_haloy
      allocate ( lat(min(l_j0+l_miny-1,1-halom):max(G_nj+halom+1,l_j0+l_maxy-1)))
      lat(1-G_haloy:G_nj+G_haloy+1)= G_yg_8(1-G_haloy:G_nj+G_haloy+1)
      do j= lbound(lat,dim=1), -G_haloy
         lat(j)= G_yg_8(1-G_haloy) + (j-1+G_haloy) * geomh_hy_8
      end do
      do j = G_nj+G_haloy+2, ubound(lat,dim=1)
         lat(j)= G_yg_8(G_nj+G_haloy) + (j-G_nj-G_haloy) * geomh_hy_8
      end do
      do j = l_miny, l_maxy
         Adz_cy_8(j) = 1.d0 / cos(lat(l_j0+j-1))
      end do
      deallocate(lat)

      end if

      allocate (adz_vwt2m(4,l_nk))
      allocate (adz_vw1t(l_nk),adz_vw2t(l_nk),adz_vw3t(l_nk),adz_vw4t(l_nk))
      allocate( Adz_zabcd_8%t(l_nk),Adz_zbacd_8%t(l_nk),Adz_zcabd_8%t(l_nk),Adz_zdabc_8%t(l_nk))
      allocate( Adz_zxabcde_8%t(l_nk),Adz_zaxbcde_8%t(l_nk),Adz_zbxacde_8%t(l_nk),&
                Adz_zcxabde_8%t(l_nk),Adz_zdxabce_8%t(l_nk),Adz_zexabcd_8%t(l_nk))
      allocate( Adz_zabcd_8%m(l_nk),Adz_zbacd_8%m(l_nk),Adz_zcabd_8%m(l_nk),Adz_zdabc_8%m(l_nk))
      allocate( Adz_zxabcde_8%m(l_nk),Adz_zaxbcde_8%m(l_nk),Adz_zbxacde_8%m(l_nk),&
                Adz_zcxabde_8%m(l_nk),Adz_zdxabce_8%m(l_nk),Adz_zexabcd_8%m(l_nk))

      do k=2, l_nk-1
         rx = Ver_z_8%t(k)
         ra = Ver_z_8%m(k-1)
         rb = Ver_z_8%m(k  )
         rc = Ver_z_8%m(k+1)
         rd = Ver_z_8%m(k+2)
         adz_vw1t(k) = lag3(rx, ra, rb, rc, rd)
         adz_vw2t(k) = lag3(rx, rb, ra, rc, rd)
         adz_vw3t(k) = lag3(rx, rc, ra, rb, rd)
         adz_vw4t(k) = lag3(rx, rd, ra, rb, rc)
      end do
      adz_vw5 = (Ver_z_8%t(0)-Ver_z_8%m(1)) / (Ver_z_8%t(0)-Ver_z_8%t(1))

      do k= 2,l_nk-2
         ra = Ver_z_8%t(k-1)
         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         rd = Ver_z_8%t(k+2)
         Adz_zabcd_8%t(k) = 1.0/triprod(ra,rb,rc,rd)
         Adz_zbacd_8%t(k) = 1.0/triprod(rb,ra,rc,rd)
         Adz_zcabd_8%t(k) = 1.0/triprod(rc,ra,rb,rd)
         Adz_zdabc_8%t(k) = 1.0/triprod(rd,ra,rb,rc)
         ra = Ver_z_8%m(k-1)
         rb = Ver_z_8%m(k)
         rc = Ver_z_8%m(k+1)
         rd = Ver_z_8%m(k+2)
         Adz_zabcd_8%m(k) = 1.0/triprod(ra,rb,rc,rd)
         Adz_zbacd_8%m(k) = 1.0/triprod(rb,ra,rc,rd)
         Adz_zcabd_8%m(k) = 1.0/triprod(rc,ra,rb,rd)
         Adz_zdabc_8%m(k) = 1.0/triprod(rd,ra,rb,rc)
      end do

      do k = 3,l_nk-3
         rx = Ver_z_8%t(k-2)
         ra = Ver_z_8%t(k-1)
         rb = Ver_z_8%t(k)
         rc = Ver_z_8%t(k+1)
         rd = Ver_z_8%t(k+2)
         re = Ver_z_8%t(k+3)
         Adz_zxabcde_8%t(k) = 1.0/quiprod(rx,ra,rb,rc,rd,re)
         Adz_zaxbcde_8%t(k) = 1.0/quiprod(ra,rx,rb,rc,rd,re)
         Adz_zbxacde_8%t(k) = 1.0/quiprod(rb,rx,ra,rc,rd,re)
         Adz_zcxabde_8%t(k) = 1.0/quiprod(rc,rx,ra,rb,rd,re)
         Adz_zdxabce_8%t(k) = 1.0/quiprod(rd,rx,ra,rb,rc,re)
         Adz_zexabcd_8%t(k) = 1.0/quiprod(re,rx,ra,rb,rc,rd)
         rx = Ver_z_8%m(k-2)
         ra = Ver_z_8%m(k-1)
         rb = Ver_z_8%m(k)
         rc = Ver_z_8%m(k+1)
         rd = Ver_z_8%m(k+2)
         re = Ver_z_8%m(k+3)
         Adz_zxabcde_8%m(k) = 1.0/quiprod(rx,ra,rb,rc,rd,re)
         Adz_zaxbcde_8%m(k) = 1.0/quiprod(ra,rx,rb,rc,rd,re)
         Adz_zbxacde_8%m(k) = 1.0/quiprod(rb,rx,ra,rc,rd,re)
         Adz_zcxabde_8%m(k) = 1.0/quiprod(rc,rx,ra,rb,rd,re)
         Adz_zdxabce_8%m(k) = 1.0/quiprod(rd,rx,ra,rb,rc,re)
         Adz_zexabcd_8%m(k) = 1.0/quiprod(re,rx,ra,rb,rc,rd)
      end do
      do k= 2, l_nk
         rx = Ver_z_8%m(k)
         ra = Ver_z_8%t(k-2)
         rb = Ver_z_8%t(k-1)
         rc = Ver_z_8%t(k)
         rd = Ver_z_8%t(k+1)
         adz_vwt2m(1,k) = lag3(rx, ra, rb, rc, rd)
         adz_vwt2m(2,k) = lag3(rx, rb, ra, rc, rd)
         adz_vwt2m(3,k) = lag3(rx, rc, ra, rb, rd)
         adz_vwt2m(4,k) = lag3(rx, rd, ra, rb, rc)
      end do
      
      allocate (Adz_uu_arr(l_ni,l_nj,l_nk), Adz_vv_arr   (l_ni,l_nj,l_nk),&
                Adz_ww_arr(l_ni,l_nj,l_nk), Adz_uvw_dep(3,l_ni,l_nj,Adz_k0m:l_nk), &
                Adz_uvw_lastl(3,l_ni,l_nj,Adz_k0m:l_nk))
      
      if (ADZ_OD_L) then
         allocate ( Adz_uvw_d(3,l_minx:l_maxx,l_miny:l_maxy,l_nk))
         allocate ( Adz_uvw_d_dep(3,l_minx:l_maxx,l_miny:l_maxy,l_nk))
      else
         allocate (orhs_extended(Adz_nij*6*l_nk)) ; orhs_extended=0.
         orhsu_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(1:)
         orhsv_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(Adz_nij*l_nk  +1:)
         orhst_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(Adz_nij*l_nk*2+1:)
         orhsc_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(Adz_nij*l_nk*3+1:)
         orhsw_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(Adz_nij*l_nk*4+1:)
         orhsf_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>orhs_extended(Adz_nij*l_nk*5+1:)
         allocate (nl_terms(Adz_nij*10*l_nk)) ; nl_terms=0.
         nlu_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(1:)
         nlv_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk  +1:)
         nlt_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*2+1:)
         nlw_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*3+1:)
         nlq_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*4+1:)
         nlu_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*5+1:)
         nlv_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*6+1:)
         nlt_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*7+1:)
         nlw_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*8+1:)
         nlq_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>nl_terms(Adz_nij*l_nk*9+1:)

         ub=0 ; dim= Adz_nij*l_nk
         allocate (rhs_bdf(dim*12+2*Adz_nij)) ; rhs_bdf=0.
         rhst_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsf_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsu_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsv_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsw_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsc_bdf_t1(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk+1)=>rhs_bdf(ub+1:); ub=ub+dim+Adz_nij
         rhst_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsf_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsu_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsv_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsw_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk)=>rhs_bdf(ub+1:); ub=ub+dim
         rhsc_bdf_t2(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk+1)=>rhs_bdf(ub+1:); ub=ub+dim+Adz_nij
         
         allocate (Adz_uvw_d(3,Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk))
         allocate (Adz_uvw_d_dep(3,Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk))
      endif

      allocate(Adz_pm (3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pmu(3,Adz_i0u:Adz_inu, Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pmv(3,Adz_i0 :Adz_in , Adz_j0v:Adz_jnv,Adz_k0 :l_nk),&
               Adz_pt (3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0t:l_nk),&
               Adz_pb (3,Adz_i0b:Adz_inb, Adz_j0b:Adz_jnb,      1:l_nk),&
               Adz_dep(3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pdu(3,Adz_i0u:Adz_inu, Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pdv(3,Adz_i0 :Adz_in , Adz_j0v:Adz_jnv,Adz_k0 :l_nk),&
               Adz_pdt(3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0t:l_nk),&
               Adz_pb2(3,Adz_i0b:Adz_inb, Adz_j0b:Adz_jnb,      1:l_nk) ) !<--for departure point

!for Adz_BC_LAM_Aranami_hlt
      allocate(Adz_o(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk,1:MAXTR3D+1),&
               Adz_i(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk,1:MAXTR3D+1))

      flag_m_f = FLAG_LVL_M
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

     ! dim= Adz_nij*6*l_nk
      dim= Adz_nij*3*l_nk
      call gmm_build_meta1D ( meta, 1,dim,0,0,dim, &
                              0,GMM_NULL_FLAGS )
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)
      istat = gmm_create('ADZ_WINDS',Adz_winds_ext,meta,flag_r_n)
      Adz_uu_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk) => Adz_winds_ext(1:)
      Adz_vv_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk) => Adz_winds_ext(Adz_nij*l_nk  +1:)
      Adz_ww_ext(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk) => Adz_winds_ext(Adz_nij*l_nk*2+1:)

      call gmm_build_meta4D (meta,  1,3   ,0,0,   3, &
                                    1,l_ni,0,0,l_ni, &
                                    1,l_nj,0,0,l_nj, &
                                    Adz_k0m,l_nk,0,0,l_nk, &
                                    0,GMM_NULL_FLAGS)
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)
      istat= gmm_create('ADZ_PXYZM',Adz_pxyzm, mymeta, flag_r_n)
      istat= gmm_create('ADZ_PXYZD',Adz_pxyzd, mymeta, flag_r_n)
      
      allocate (Adz_wpxyz(-1:l_ni+2,-1:l_nj+2,l_nk+1,3),Adz_dpxyz(-1:l_ni+2,-1:l_nj+2,l_nk+1,3))
      allocate (Adz_wpz(l_ni,l_nj,l_nk+1),Adz_dpz(l_ni,l_nj,l_nk+1))

      call gmm_build_meta4D (meta,  -1,l_ni+2,0,0,l_ni+4, &
                                    -1,l_nj+2,0,0,l_nj+4, &
                                    1,l_nk,0,0,l_nk, &
                                    1,3   ,0,0,   3, &
                                    0,GMM_NULL_FLAGS)
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)
      !for displacement array in trajectories
      istat= gmm_create('ADZ_DISP',Adz_disp, mymeta, flag_r_n)
      
!     Wind information at the lowest thermodynamic level
!     from the physics surface layer scheme
      nullify(Adz_uslt,Adz_vslt)
      if ( .not. Ctrl_phyms_L ) Adz_slt_winds= .false.
      if (Adz_slt_winds) then
         istatu= gmm_get(gmmk_pw_uslt_s, Adz_uslt)
         istatv= gmm_get(gmmk_pw_vslt_s, Adz_vslt)
         if ((istatu/=0).or.(istatv/=0)) Adz_slt_winds= .false.
      end if

      allocate (pnt_stack(max(10,Tr3d_ntr)))
      allocate (Adz_stack(max(10,Tr3d_ntr)))

      !Allocation Conservation Postprocessing
      !--------------------------------------
      allocate (adz_flux_3CWP(max(1,Tr3d_ntrTRICUB_WP)))
      allocate (adz_flux_BQWP(max(1,Tr3d_ntrBICHQV_WP)))
      allocate (adz_flux_3CWP_PS(1))

      do n = 1,max(1,Tr3d_ntrTRICUB_WP)

         allocate(adz_flux_3CWP(n)%fo(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_flux_3CWP(n)%fi(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      do n = 1,max(1,Tr3d_ntrBICHQV_WP)

         allocate(adz_flux_BQWP(n)%fo(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_flux_BQWP(n)%fi(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      allocate(adz_flux_3CWP_PS(1)%fo(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      allocate(adz_flux_3CWP_PS(1)%fi(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      allocate (adz_post_3CWP(max(1,Tr3d_ntrTRICUB_WP)))
      allocate (adz_post_BQWP(max(1,Tr3d_ntrBICHQV_WP)))

      do n = 1,max(1,Tr3d_ntrTRICUB_WP)

         allocate(adz_post_3CWP(n)%lin(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_3CWP(n)%min(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_3CWP(n)%max(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      do n = 1,max(1,Tr3d_ntrBICHQV_WP)

         allocate(adz_post_BQWP(n)%lin(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_BQWP(n)%min(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_BQWP(n)%max(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      allocate (adz_bc(1+Tr3d_ntrTRICUB_WP+Tr3d_ntrBICHQV_WP))

      do n = 1,1+Tr3d_ntrTRICUB_WP+Tr3d_ntrBICHQV_WP

         allocate(adz_bc(n)%w(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      Adz_niter = Adz_itraj
      if ( .not. Rstri_rstn_L ) call adz_inittraj
!      
!---------------------------------------------------------------------
!
      return

 1000 format(/,'=================================================',/,A31,I2,/,&
               '=================================================' )

contains

      real(kind=REAL64) function triprod(za,zb,zc,zd)
         real(kind=REAL64), intent(in) :: za,zb,zc,zd
         triprod = ((za-zb)*(za-zc)*(za-zd))
      end function triprod

      real(kind=REAL64) function quiprod(zx,za,zb,zc,zd,ze)
         real(kind=REAL64), intent(in) :: zx,za,zb,zc,zd,ze
         quiprod = ((zx-za)*(zx-zb)*(zx-zc)*(zx-zd)*(zx-ze))
      end function quiprod

      real(kind=REAL64) function lag3(zz, z1, z2, z3, z4)
      implicit none
      real(kind=REAL64), intent(in)  :: zz, z1, z2, z3, z4
      
      lag3 = ( ((zz - z2) * (zz - z3) * (zz - z4) ) / &
               ((z1 - z2) * (z1 - z3) * (z1 - z4) ) )
      end function lag3

      end subroutine adz_set
