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
!
!	same as adz_main_h, expect the source for the interpolation is different.
!	Trajectories code and staggering is the same.
!
!
!------------------------------------------------------------------------------
      subroutine SW_adz_will (F_dt_8)
      use ISO_C_BINDING
      use dyn_fisl_options
      use cstv
      use gem_options
      use adz_mem
      use adz_interp_mod
      use dynkernel_options
      use mem_tstp
      use gmm_vt0
      use gmm_vt1
      use gmm_vt2
      use step_options
      implicit none

      real(kind=REAL64), intent(IN) :: F_dt_8

      integer :: n,i,j,k
      integer :: i0, j0, k0, in, jn, k0t
      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      real(kind=REAL64) :: dtA_8,dtzA_8,dtD_8,dtzD_8,invT_n_8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
      type(Adz_pntr_stack), dimension(3), target :: stack
!
!     ---------------------------------------------------------------

      i0= 1   +pil_w
      in= l_ni-pil_e
      j0= 1   +pil_s
      jn= l_nj-pil_n
      k0= 1
      k0t=k0
      invT_n_8 = 3.0/(2.0*F_dt_8)
    
      call SW_disp_traject (Cstv_dt_8)

!!$omp do collapse(2)
!      do k=k0, l_nk
!         do j= j0-1, jn    !following the indices of nli terms
!            do i= i0-1, in !following the indices of nli terms
!            orhsw_ext(i,j,k)  = qt1(i,j,k)
!            orhsf_ext(i,j,k)  = qt2(i,j,k)
!            end do
!         end do
!      end do
!!$omp enddo


      HLT_j0 = 1
      HLT_jn = l_nj
      HLT_nk = l_nk
      HLT_nj = HLT_jn - HLT_j0 + 1
      call HLT_split (1, HLT_nj*HLT_nk, HLT_np, HLT_start, HLT_end)
      do n= HLT_start, HLT_end
         k= (n-1)/HLT_nj
         j= n - k*HLT_nj + HLT_j0 - 1
         k= k+1

         do i= 1, l_ni
            orhsw_ext(i,j,k)  = qt1(i,j,k)
            orhsf_ext(i,j,k)  = qt2(i,j,k)
        end do
      end do


!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo (  orhsu_ext(Adz_lminx,Adz_lminy,HLT_start),&
                Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, HLT_np,-1)

!==========RHS terms advection======

!---midpoint---
      stack(1)%src => orhsw_ext
      stack(1)%dst => rhsc_mid
      call adz_tricub ( stack,1,Adz_pm ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

!---departure---
      stack(1)%src => orhsf_ext
      stack(1)%dst => rhsc_dep
      call adz_tricub ( stack,1,Adz_dep ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

!$OMP BARRIER


      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_mid(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_dep(l_minx,l_miny,HLT_start),&
                  l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)



!***********************************************************
! Compute the rhs terms that come from the bdf method      *
!***********************************************************

!!$omp do collapse(2)
      do k=k0, l_nk
         do j= j0-1, jn    
            do i= i0-1, in 
              !---rhsc---
               qt0(i,j,k) = (4.0/3.0)*rhsc_mid(i,j,k) - (one/3.0)*rhsc_dep(i,j,k) 
              !qt0(i,j,k) = rhsc_mid(i,j,k)
            end do
         end do
      end do
!!$omp enddo


!
!     ---------------------------------------------------------------
!
      return
      end subroutine SW_adz_will
