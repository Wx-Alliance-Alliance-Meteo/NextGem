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

!**s/r nest_bcs_t0 - Apply horizontal boundary conditions

      subroutine nest_bcs_t0  ( F_invT, F_rhsu, F_rhsv, &
                                Minx, Maxx, Miny, Maxy, Nk )
      use theo_options
      use lam_options
      use mem_nest
      use mem_tracers
      use gmm_vt0
      use sol_mem
      use glb_ld
      use tr3d
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: &
                                               F_rhsu, F_rhsv
      real(kind=REAL64), intent(IN) :: F_invT

      integer j,k,n,deb
!
!----------------------------------------------------------------------
!
      if (Theo_periodicX_L) then
         call periodicX ()
         do k=1,G_nk
            if (l_north) F_rhsv (1+pil_w:l_ni-pil_e,l_nj-pil_n:l_nj,k) = F_invT * nest_v(1+pil_w:l_ni-pil_e,l_nj-pil_n:l_nj,k)
            if (l_south) F_rhsv (1+pil_w:l_ni-pil_e,1:pil_s,k) = F_invT * nest_v(1+pil_w:l_ni-pil_e,1:pil_s,k)
            if (l_east ) F_rhsu (l_ni-pil_e:l_ni,1+pil_s:l_nj-pil_n,k) = F_invT * ut0(l_ni-pil_e:l_ni,1+pil_s:l_nj-pil_n,k)
            if (l_west ) F_rhsu (1:pil_w:l_ni,1+pil_s:l_nj-pil_n,k) = F_invT * ut0(1:pil_w:l_ni,1+pil_s:l_nj-pil_n,k)
         end do
      else
!!$omp do
      do n=1,Tr3d_ntr
         deb= (n-1)*G_nk + 1

         if (l_north) &
         tracers_M(n)%pntr (1:l_ni ,l_nj-pil_n+1:l_nj ,1:G_nk) = &
         nest_tr(1:l_ni,l_nj-pil_n+1:l_nj,deb:deb+G_nk-1)

         if (l_east ) &
         tracers_M(n)%pntr (l_ni-pil_e+1:l_ni ,1:l_nj ,1:G_nk) = &
         nest_tr(l_ni-pil_e+1:l_ni,1:l_nj,deb:deb+G_nk-1)

         if (l_south) &
         tracers_M(n)%pntr (1:l_ni ,1:pil_s ,1:G_nk) = &
         nest_tr (1:l_ni,1:pil_s,deb:deb+G_nk-1)

         if (l_west ) &
         tracers_M(n)%pntr (1:pil_w ,1:l_nj ,1:G_nk) = &
         nest_tr (1:pil_w,1:l_nj,deb:deb+G_nk-1)

         if (Schm_opentop_L) then
            tracers_M(n)%pntr (1:l_ni, 1:l_nj, 1:Lam_gbpil_t-1) = &
            nest_tr (1:l_ni,1:l_nj,deb:deb+Lam_gbpil_t-2)
         end if

      enddo
!!$omp enddo nowait

!!$omp do
      do k=1, l_nk
         if (l_north) then
            ut0  (1:l_niu,l_nj-pil_n+1:l_nj ,k) = nest_u (1:l_niu,l_nj-pil_n+1:l_nj ,k)
            vt0  (1:l_ni ,l_nj-pil_n  :l_njv,k) = nest_v (1:l_ni ,l_nj-pil_n  :l_njv,k)
            wt0  (1:l_ni ,l_nj-pil_n+1:l_nj ,k) = nest_w (1:l_ni ,l_nj-pil_n+1:l_nj ,k)
            tt0  (1:l_ni ,l_nj-pil_n+1:l_nj ,k) = nest_t (1:l_ni ,l_nj-pil_n+1:l_nj ,k)
            qt0  (1:l_ni ,l_nj-pil_n+1:l_nj ,k) = nest_q (1:l_ni ,l_nj-pil_n+1:l_nj ,k)
            zdt0 (1:l_ni ,l_nj-pil_n+1:l_nj ,k) = nest_zd(1:l_ni ,l_nj-pil_n+1:l_nj ,k)
            Sol_lhs(1:l_ni,l_nj-pil_n+1:l_nj,k) =     qt0(1:l_ni ,l_nj-pil_n+1:l_nj ,k)
            F_rhsv (1+pil_w:l_ni-pil_e,l_nj-pil_n:l_nj,k) = F_invT * nest_v(1+pil_w:l_ni-pil_e,l_nj-pil_n:l_nj,k)
         endif
         if (l_south) then
            ut0 (1:l_niu,1:pil_s ,k) = nest_u (1:l_niu,1:pil_s ,k)
            vt0 (1:l_ni ,1:pil_s ,k) = nest_v (1:l_ni ,1:pil_s ,k)
            wt0 (1:l_ni ,1:pil_s ,k) = nest_w (1:l_ni ,1:pil_s ,k)
            tt0 (1:l_ni ,1:pil_s ,k) = nest_t (1:l_ni ,1:pil_s ,k)
            qt0 (1:l_ni ,1:pil_s ,k) = nest_q (1:l_ni ,1:pil_s ,k)
            zdt0(1:l_ni ,1:pil_s ,k) = nest_zd(1:l_ni ,1:pil_s ,k)
            Sol_lhs(1:l_ni,1:pil_s,k)=     qt0(1:l_ni ,1:pil_s ,k)
            F_rhsv (1+pil_w:l_ni-pil_e,1:pil_s,k) = F_invT * nest_v(1+pil_w:l_ni-pil_e,1:pil_s,k)
         endif
         if (l_east) then
            ut0 (l_ni-pil_e  :l_niu,1:l_nj ,k) = nest_u (l_ni-pil_e  :l_niu,1:l_nj ,k)
            vt0 (l_ni-pil_e+1:l_ni ,1:l_njv,k) = nest_v (l_ni-pil_e+1:l_ni ,1:l_njv,k)
            tt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,k) = nest_t (l_ni-pil_e+1:l_ni ,1:l_nj ,k)
            wt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,k) = nest_w (l_ni-pil_e+1:l_ni ,1:l_nj ,k)
            qt0 (l_ni-pil_e+1:l_ni ,1:l_nj ,k) = nest_q (l_ni-pil_e+1:l_ni ,1:l_nj ,k)
            zdt0(l_ni-pil_e+1:l_ni ,1:l_nj ,k) = nest_zd(l_ni-pil_e+1:l_ni ,1:l_nj ,k)
            Sol_lhs(l_ni-pil_e+1:l_ni,1:l_nj,k)=     qt0(l_ni-pil_e+1:l_ni ,1:l_nj ,k)
            F_rhsu (l_ni-pil_e:l_ni,1+pil_s:l_nj-pil_n,k) = F_invT * nest_u(l_ni-pil_e:l_ni,1+pil_s:l_nj-pil_n,k)
         endif
         if (l_west) then
            ut0 (1:pil_w, 1:l_nj , k) = nest_u (1:pil_w, 1:l_nj , k)
            vt0 (1:pil_w, 1:l_njv, k) = nest_v (1:pil_w, 1:l_njv, k)
            tt0 (1:pil_w, 1:l_nj , k) = nest_t (1:pil_w, 1:l_nj , k)
            wt0 (1:pil_w, 1:l_nj , k) = nest_w (1:pil_w, 1:l_nj , k)
            qt0 (1:pil_w, 1:l_nj , k) = nest_q (1:pil_w, 1:l_nj , k)
            zdt0(1:pil_w, 1:l_nj , k) = nest_zd(1:pil_w, 1:l_nj , k)
            Sol_lhs(1:pil_w,1:l_nj,k) =     qt0(1:pil_w, 1:l_nj , k)
            F_rhsu (1:pil_w,1+pil_s:l_nj-pil_n,k) = F_invT * nest_u(1:pil_w,1+pil_s:l_nj-pil_n,k)
         endif
      end do
!!$omp enddo nowait
      
      if (Lam_toptt_L) then
!        Pilot the temperature for the whole top level
!$OMP BARRIER
!!$omp do
         do j=1, l_nj
            tt0(1:l_ni,j,1) = nest_t(1:l_ni,j,1)
         end do
!!$omp enddo  nowait
      end if
!$OMP BARRIER

      endif
!
!----------------------------------------------------------------------
!
      return
      end
