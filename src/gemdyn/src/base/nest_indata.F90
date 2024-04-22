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

!**s/r nest_indata - Read and process nesting data during LAM
!                    integration for LBC.

      subroutine nest_indata ( F_u , F_v, F_w, F_t , F_q, F_zd, F_s, &
                               F_tr, F_topo, F_stag_L, F_datev_S    ,&
                               Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr )
      use dyn_fisl_options
      use gem_options
      use inp_mod
      use mem_tstp
      use mem_nest
      use glb_ld
      use tr3d
      implicit none
      
      character(len=*), intent(in):: F_datev_S
      logical, intent(in) :: F_stag_L
      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk),   intent(out) :: F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), intent(out) :: F_q
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy),      intent(out) :: F_s
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,2),    intent(out) :: F_topo
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk*Ntr),intent(out) :: F_tr

      integer :: i,j,k,dim
      integer :: HLT_start, HLT_end, local_np
      real, dimension(:,:  ), pointer :: nest_orols
      real, dimension(:,:,:), pointer :: uu, vv, tt, sumpqj
!     
!     ---------------------------------------------------------------
!
      dim=  (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
      uu    (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(      1:)
      vv    (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(  dim+1:)
      tt    (l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(2*dim+1:)
      sumpqj(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => WS1(3*dim+1:)
      nest_orols(l_minx:l_maxx,l_miny:l_maxy)    => WS1(4*dim+1:)
      
!!$omp single
      Inp_src_GZ_L= .false. ; Inp_gtmg= (/85,90/)
      call inp_data (F_u , F_v, F_w, F_t , F_q, F_zd, F_s, F_tr     ,&
                     F_topo(l_minx,l_miny,1),F_topo(l_minx,l_miny,2),&
              F_stag_L,F_datev_S,l_minx,l_maxx,l_miny,l_maxy,G_nk,Ntr)
!!$omp end single

      if (Schm_sleve_L) then
!!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               nest_orols(i,j) = F_topo(i,j,2)
            enddo
         end do
!!$omp enddo

      else
!!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               nest_orols(i,j) = 0.
            enddo
         end do
!!$omp enddo
      endif

      call vertical_metric_omp (nest_metric, F_topo, nest_orols, &
                                l_minx,l_maxx,l_miny,l_maxy)
                   
      call canonical_indata()
      
      if (.not.F_stag_L) then

!!$omp do collapse(2)
         do k=1, l_nk
            do j=1-G_haloy, l_nj+G_haloy
               do i=1-G_halox, l_ni+G_halox
                  uu(i,j,k) = F_u(i,j,k)
                  vv(i,j,k) = F_v(i,j,k)
                  tt(i,j,k) = F_t(i,j,k)
                  sumpqj(i,j,k) = 0.
            enddo
         end do
         end do
!!$omp enddo

!!$omp single
         call mfottvh2 (tt, F_t, F_tr(l_minx,l_minx,(Tr3d_hu-1)*l_nk+1),&
                        sumpqj,l_minx, l_maxx, l_miny, l_maxy, l_nk    ,&
                   1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,.true.)
!!$omp end single
         call HLT_split (1, 2*G_nk, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( uu(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
                    
         call hwnd_stag2 ( F_u, F_v,uu,vv                  ,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk  ,&
                         1-G_halox*west ,l_niu+G_halox*east,&
                         1-G_haloy*south,l_njv+G_haloy*north, .true. )
      endif
                       
      call derivate_data ( F_zd, F_w, F_u, F_v, F_t , F_s, F_q           ,&
                           F_topo(l_minx,l_miny,1),nest_orols,nest_metric,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk             ,&
                           .not.Inp_zd_L, .not.Inp_w_L, .not. Inp_qt_L )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine nest_indata

      subroutine nest_indata_svr ( F_u , F_v, F_w, F_t , F_q, F_zd, F_s, &
                               F_tr, F_topo, F_datev_S    ,&
                               Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr )
      use, intrinsic :: iso_fortran_env
      use dyn_fisl_options
      use gem_options
      use step_options
      use inp_mod
      use lun
      use mem_tstp
      use mem_nest
      use glb_ld
      use svri_mod
      use ptopo
      implicit none
      
      character(len=*), intent(in):: F_datev_S
      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk),   intent(out) :: F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), intent(out) :: F_q
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy),      intent(out) :: F_s
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,2),    intent(out) :: F_topo
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk*Ntr),intent(out) :: F_tr

      include 'mpif.h'
      character(len=16) :: Next_pilot_S
      integer :: i,j,k,err
      integer :: HLT_start, HLT_end, local_np
      real, dimension(:,:  ), pointer :: nest_orols
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: one=1.0d0, sid=86400.0d0, rsid=one/sid
!     
!     ---------------------------------------------------------------
!             
!!$omp single
      Inp_src_GZ_L= .false. ; Inp_gtmg= (/85,90/)
      if (INs_server_L) then
         call gtmg_start (81, 'Recv_data', 30 )
         call gemtime ( Lun_out, 'waitall', .false. )
         call MPI_waitall (size(INs_irecv),INs_irecv,&
                           MPI_STATUSES_IGNORE,err)
         call gemtime ( Lun_out, 'BARRIER DONE', .false. )
         call itf_Iserv_GZ3d ()
         if (Lun_out > 0) print*, 'INs-status1: ',Inp_src_GZ_L
         call gtmg_stop (81)
         dayfrac = Step_nesdt*rsid
         call incdatsd (Next_pilot_S, F_datev_S, dayfrac)
         if (Next_pilot_S<=Step_runend_S) then
            call itf_Iserv_request (Next_pilot_S,INs_list_S,1001,INs_nrequests)
         endif
      endif
      
      call inp_data (F_u , F_v, F_w, F_t , F_q, F_zd, F_s, F_tr     ,&
                     F_topo(l_minx,l_miny,1),F_topo(l_minx,l_miny,2),&
                .true.,F_datev_S,l_minx,l_maxx,l_miny,l_maxy,G_nk,Ntr)
!!$omp end single

      nest_orols(l_minx:l_maxx,l_miny:l_maxy) => WS1(1:)
      if (Schm_sleve_L) then
!!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               nest_orols(i,j) = F_topo(i,j,2)
            enddo
         end do
!!$omp enddo

      else
!!$omp do
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               nest_orols(i,j) = 0.
            enddo
         end do
!!$omp enddo
      endif

      call vertical_metric_omp (nest_metric, F_topo, nest_orols, &
                                l_minx,l_maxx,l_miny,l_maxy)
                   
      call canonical_indata()
                             
      call derivate_data ( F_zd, F_w, F_u, F_v, F_t , F_s, F_q           ,&
                           F_topo(l_minx,l_miny,1),nest_orols,nest_metric,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk             ,&
                           .not.Inp_zd_L, .not.Inp_w_L, .not. Inp_qt_L )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine nest_indata_svr

