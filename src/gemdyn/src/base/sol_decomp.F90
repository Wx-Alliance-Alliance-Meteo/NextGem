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

!**s/r sol_transpose - Establish memory layout for the elliptic solver

      integer function sol_decomp ( F_npex, F_npey, F_checkparti_L )
      use glb_ld
      use ldnh
      use lun
      implicit none

      logical, intent(in) :: F_checkparti_L
      integer, intent(in) :: F_npex, F_npey

      logical, external :: decomp
      integer, parameter :: lowest = 2
      integer, dimension(F_npex) :: lnis, nis
      integer, dimension(F_npey) :: lnjs
      integer :: npartiel
!
!     ---------------------------------------------------------------
!
      sol_decomp = 0

! Establishing local dimensions and computing arena (data topology) for:
!          G_ni distributed on Ptopo_npex PEs and
!          G_nj distributed on Ptopo_npey PEs both without halo

      if (Lun_out > 0) then
         write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                             ' G_ni distributed on F_npex PEs', G_ni,F_npex
      end if

      if (.not. decomp (G_ni, ldnh_minx, ldnh_maxx, lnis, npartiel, 0, ldnh_i0, &
               .true., .true., F_npex, lowest, F_checkparti_L, 0 )) sol_decomp = -1
      ldnh_ni= lnis(1)

      if (Lun_out > 0) then
         write(Lun_out,1002) ' Memory layout for SOLVER (no halo):', &
                             ' G_nj distributed on F_npey PEs', G_nj,F_npey
      end if

      if (.not. decomp (G_nj, ldnh_miny, ldnh_maxy, lnjs, npartiel, 0, ldnh_j0,&
               .false.,.true., F_npey, lowest, F_checkparti_L, 0 )) sol_decomp = -1
      ldnh_nj= lnjs(1)

 1002 format (a/a45,i6,' /',i5)
!
!     ---------------------------------------------------------------
!
      return
      end function sol_decomp

