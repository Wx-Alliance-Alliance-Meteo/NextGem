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
      subroutine adz_main (F_dt_8, itpc, F_euler_L, F_type_S)
      use adz_mem
      use mem_tstp
      implicit none

      character(len=*), intent(IN) :: F_type_S
      logical, intent(IN) :: F_euler_L
      integer, intent(IN) :: itpc
      real(kind=REAL64), intent(IN) :: F_dt_8
      
      integer :: n,i,j,k, HLT_np, HLT_start, HLT_end, order
!
!     ---------------------------------------------------------------
!
      call adz_traject_BDF2 (F_dt_8,itpc,F_euler_L)

      if ( F_type_S == 'TURBO' ) then
         call adz_SL_turbo ()
      else
         read (F_type_S, '(i1)') order
         call adz_SL_intrp (order)
      endif

!$OMP BARRIER

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_mid(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)

      call HLT_split (1, 6*l_nk, HLT_np, HLT_start, HLT_end)
      call gem_xch_halo ( rhsu_dep(l_minx,l_miny,HLT_start),&
                    l_minx,l_maxx,l_miny,l_maxy, HLT_np,-1)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_main
