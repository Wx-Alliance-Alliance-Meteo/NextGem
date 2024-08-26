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

!**s/r set_sol

      subroutine set_sol
      use gem_options
      use HORgrid_options
      use glb_pil
      use lam_options
      use dyn_fisl_options
      use glb_ld
      use lun
      use ldnh
      use mem_tstp
      use sol_mem
      use sol_options
      use ldnh
      use ptopo
      use gmm_vt1
      use gmm_table
      use opr
      use prec
      use yyg_param
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none

      type(gmm_metadata) :: meta
      integer i,j,k,ni,nj,istat,dim,dimH
      integer :: f1,f2,f3,f4, nx1, nx2, ny1, ny2
      integer :: HLT_start, HLT_end, local_np
      real(kind=REAL64) yg_8(G_nj)
!     __________________________________________________________________
!
      ni=ldnh_maxx-ldnh_minx+1
      nj=ldnh_maxy-ldnh_miny+1
      call gmm_build_meta4D (meta,&
                             l_minx,l_maxx,G_halox,G_halox,l_ni, &
                             l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                             0,l_nk+1,0,0,l_nk+2,&
                             0,0,0,0,0,0,GMM_NULL_FLAGS)
      istat= gmm_create('SOL_LHS',Sol_lhs,meta, GMM_FLAG_RSTR+GMM_FLAG_IZER)
      gmm_cnt=gmm_cnt+1 ; GMM_tbl%vname(gmm_cnt)='SOL_LHS' ; GMM_tbl%ara(gmm_cnt)='QQ' ; GMM_tbl%cn(gmm_cnt)='MM' ; GMM_tbl%fst(gmm_cnt)='SOLS'
      do k=1, l_nk+1
         Sol_lhs(:,:,k)= qt1(:,:,k)
      end do
      if ( Grd_yinyang_L) then
         call yyg_xchng_8 (Sol_lhs, YYG_HALO_q2q,l_minx,l_maxx,l_miny,l_maxy, &
                           l_ni,l_nj, l_nk+2, .false., 'CUBIC', .true.)
      else
         call HLT_split (0, G_nk+1, local_np, HLT_start, HLT_end)
         call gem_xch_halo_8 ( Sol_lhs(l_minx,l_miny,HLT_start),l_minx,l_maxx,&
                               l_miny,l_maxy,local_np,-1 )
      endif
      
      allocate (Sol_rhs(ni,nj,l_nk)) ; Sol_rhs= 0.
         
      if (Lun_out > 0) write (Lun_out,1002) trim(Sol_krylov3D_S), trim(Sol_precond3D_S)

      sol_pil_w= pil_w ; sol_pil_e= pil_e
      sol_pil_s= pil_s ; sol_pil_n= pil_n
      Sol_i0 = 1    + pil_w
      Sol_in = l_ni - pil_e
      Sol_j0 = 1    + pil_s
      Sol_jn = l_nj - pil_n      
      Sol_k0 = 1+Lam_gbpil_T
      Sol_nk = l_nk-Sol_k0+1
!! Bloc extension == restrictive Jacobi preconditionner

      if (allocated(Prec_xevec_8)) deallocate (Prec_xevec_8)
      if (allocated(Prec_xeval_8)) deallocate (Prec_xeval_8)
      if (allocated(Prec_ai_8)) deallocate (Prec_ai_8)
      if (allocated(Prec_bi_8)) deallocate (Prec_bi_8)
      if (allocated(Prec_invbi_8)) deallocate (Prec_invbi_8)
      if (allocated(Prec_ci_8)) deallocate (Prec_ci_8)
      
      Sol_ii0  = 1    - Sol_ovlpx
      Sol_iin  = l_ni + Sol_ovlpx
      Sol_jj0  = 1    - Sol_ovlpy
      Sol_jjn  = l_nj + Sol_ovlpy
      
      Sol_imin = Sol_ii0
      Sol_imax = Sol_iin
      Sol_jmin = Sol_jj0
      Sol_jmax = Sol_jjn
      Sol_dimx = Sol_iin - Sol_ii0 + 1
      sol_dimy = Sol_jjn - Sol_jj0 + 1
         
      f1 = G_ni/Ptopo_npex + min(1,mod(G_ni,Ptopo_npex))
      f2 = G_ni-f1*(Ptopo_npex-1)
      f3 = G_nj/Ptopo_npey + min(1,mod(G_nj,Ptopo_npey))
      f4 = G_nj-f3*(Ptopo_npey-1)
         
      nx1 = l_ni - glb_pil_w 
      nx2 = f2 - glb_pil_e 
      ny1 = l_nj - glb_pil_s 
      ny2 = f4 - glb_pil_n 
         
      if (Ptopo_mycol==1)  then
         Sol_ii0  = 1 -  min(Sol_ovlpx,nx1)
      endif
      if (Ptopo_myrow==1)  then
         Sol_jj0  = 1 -  min(Sol_ovlpy,ny1)
      endif
      
      if (Ptopo_mycol.eq.Ptopo_npex-2) then
         Sol_iin = l_ni + min(Sol_ovlpx,nx2)
      endif
      if (Ptopo_myrow.eq.Ptopo_npey-2) then
         Sol_jjn = l_nj + min(Sol_ovlpy,ny2)
      endif

         
      if (l_west)  Sol_ii0 = 1    + sol_pil_w
      if (l_east)  Sol_iin = l_ni - sol_pil_e
      if (l_south) Sol_jj0 = 1    + sol_pil_s
      if (l_north) Sol_jjn = l_nj - sol_pil_n
      
      sol_niloc=Sol_iin-Sol_ii0+1
      sol_njloc=Sol_jjn-Sol_jj0+1
      sol_nloc = sol_niloc*sol_njloc*Schm_nith
      
      allocate (Prec_xevec_8(sol_niloc,sol_niloc)  ,&
                Prec_xeval_8(sol_niloc)            ,&
                Prec_ai_8(sol_niloc,sol_njloc,G_nk),&
                Prec_bi_8(sol_niloc,sol_njloc,G_nk),&
                Prec_invbi_8(sol_niloc,sol_njloc,G_nk),&
                Prec_ci_8(sol_niloc,sol_njloc,G_nk))
       call eigenabc_local (Prec_xeval_8,Prec_xevec_8,Prec_ai_8,&
                            Prec_bi_8,Prec_invbi_8,Prec_ci_8   ,&
                            sol_niloc,sol_njloc,Schm_nith)

       allocate (gg(1:sol_im+1),rot_cos(1:sol_im+1), rot_sin(1:sol_im+1), IPIV_arr(1:sol_im+1))
       allocate (v_lcl_sum(1:sol_im+1,1:2),rr(1:sol_im+1,1:sol_im+1),&
                 tt(1:sol_im+1,1:sol_im+1), hessenberg(1:sol_im+1, 1:sol_im))!,&
                ! N_mat(1:sol_im+1,1:sol_im+1), M_mat(1:sol_im+1, 1:sol_im+1),&
                ! T_mat1(1:sol_im+1,1:sol_im+1))
       allocate (work_space(sol_imin:sol_imax,sol_jmin:sol_jmax,1:l_nk)   ,&
                 vv(sol_imin:sol_imax,sol_jmin:sol_jmax,1:l_nk,1:sol_im+1),&
                 wint_8(sol_ii0:sol_iin,sol_jj0:sol_jjn,1:l_nk,1:sol_im+1))
       allocate (thread_s   (1:4,0:OMP_get_max_threads()-1),&
                 thread_s128(1:4,0:OMP_get_max_threads()-1),&
                 thread_s2(1:2,1:sol_im+1,0:OMP_get_max_threads()-1))

      !needed for matvec lateral boundary conditions
       allocate (m_west(l_ni),m_east(l_ni),m_south(l_nj),m_north(l_nj))
       m_west=1.d0 ; m_east=1.d0 ; m_south=1.d0 ; m_north=1.d0
       if (.not.Grd_yinyang_L) then
          if (l_west ) m_west(1+pil_w    ) = 0.
          if (l_east ) m_east(l_ni-pil_e ) = 0.
          if (l_south) m_south(1+pil_s   ) = 0.
          if (l_north) m_north(l_nj-pil_n) = 0.
       endif

       dimH= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)
       dim = dimH*G_nk
       allocate (sol_ws(5*dim))
       dqdzu(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => sol_ws(      1:)
       dqdzv(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => sol_ws(  dim+1:)
       Qu   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => sol_ws(2*dim+1:)
       Qv   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => sol_ws(3*dim+1:)
       Qw   (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => sol_ws(4*dim+1:)

       ni= Sol_iin-Sol_ii0+1
       nj= Sol_jjn-Sol_jj0+1
       allocate (fdg(ni,nj,l_nk),w2_8(ni,nj,l_nk),w3_8(ni,nj,l_nk))
       allocate (ext_q(l_minx:l_maxx,l_miny:l_maxy,0:l_nk+1))
                 
       allocate (fdg2(l_minx:l_maxx,l_miny:l_maxy,l_nk+1))

       work_space= 0.
       wint_8=0.
       vv= 0.
       ext_q=0. ; fdg2=0.
       thread_s= 0.
       thread_s2= 0.

 1002 format(/,'WILL USE ',a,' 3D ITERATIVE SOLVER WITH ',a, &
               ' PRECONDITIONNER' &
             /,'=============================================')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine set_sol
