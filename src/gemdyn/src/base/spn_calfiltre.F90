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
!---------------------------------- LICENCE END --------------------------------

!*s/r spn_calfiltre - compute a filter for spectral nudging
! Updated 2023 May, Christopher Subich -- refactor for separate transpose module
!                                         This refactor pre-removes the piloting region before computing
!                                         the parallel transpose, and it uses a maximally-flexible process
!                                         decomposition in the transposes.  The net effect is that the size
!                                         of the doubly-transposed array has changed, requiring adjustments
!                                         to the bounds of (but not the meaning of) the spectral filter 
!                                         coefficients

      subroutine spn_calfiltre
      use dcst
      use glb_ld
      use glb_pil
      use HORgrid_options
      use spn_options
      use tdpack
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none

      integer i,j,il,jl
      integer ni_trunc, ni_truncx, nj_trunc
      real(kind=REAL64) nix, njx, nkx, nk_cut
      real(kind=REAL64) WXL, WXS, DX, DY
!
!----------------------------------------------------------------------
!
      ! The spectral filter here works on the Cartesian lat/lon grid, rather
      ! than in a projected space with physical distances
      DX = Grd_dx*pi_8*Dcst_rayt_8/(180.*1000.) ! delta-x at equator, in kilometers
      DY = Grd_dy*pi_8*Dcst_rayt_8/(180.*1000.) ! delta-y at equator, in kilometers
      ! Config-specified cutoff scales, in kilometers
      WXL= Spn_cutoff_scale_large
      WXS= Spn_cutoff_scale_small

      if (Lun_out > 0) write(Lun_out,1000) WXL,WXS

      ! Integral wavenumber corresponding to the large-scale cutoff, in x (i) and y (j),
      ! used for normalization
      ni_trunc = int(DX*(G_ni-2*Grd_extension)/WXL)
      nj_trunc = int(DY*(G_nj-2*Grd_extension)/WXL)

      ! Integral wavenumber corresponding to the small-scale cuttof, in x (i) only
      ni_truncx= int(DX*(G_ni-2*Grd_extension)/WXS)

      ! Normalized cutoff frequency
      nk_cut = float(ni_truncx)/float(ni_trunc)

      ! Spn_flt is allocated in spn_init, so it should be prepared already.
      ! allocate (Spn_flt(Spn_22n,G_nj))
      ! Spn_flt= 0.

      ! The filter coefficients are defined in wavenumber space, and they ''live''
      ! on the y-contiguous grid (after the double transpose), with (x,y) order.
      do j = 1, Spn_ygrid_lny
         jl = j-1 ! Integral j-wavenumber bucket, with 0 being the constant frequency
         do i = 1, Spn_ygrid_lnx
            ! The integral i-wavenumber bucket is offset by Spn_ygrid_llbx, since
            ! this process probably doesn't whon the whole x-grid.
            il = i - 1 + Spn_ygrid_llbx 
            ! Compute scaled i and j wavenumbers
            nix = dble(il)/dble(ni_trunc)
            njx = dble(jl)/dble(nj_trunc)
            ! Compute the distance from the origin, in wavenumber space
            nkx = sqrt(nix*nix + njx*njx)
            ! nkx = 1.0 when the scaled wavenumber is equal to the low-frequency cutoff,
            ! and it is equal to the (higher) nk_cut value at the high-frequency cutoff.

            if (nkx > nk_cut) then
               ! If this value is above the high-frequency cutoff, the filter coefficient is 0
               Spn_flt(i,j) = 0.0
            elseif (nkx > 1.0) then
               ! If this value is between the low and high cutoff scales, smoothly interpolate
               ! between the two with a cos^2
               Spn_flt(i,j) = cos((pi_8/2.0) * (nkx-1.0)/(nk_cut-1.0))**2
            else
               ! Otherwise, the filter coefficient is 1.0
               Spn_flt(i,j) = 1.0
            end if
         end do
      end do

      ! do j=1+Grd_extension,G_nj-Grd_extension
      !    jl= j-Grd_extension-1
      !    do i=1+Spn_22pil_w,Spn_22n-Spn_22pil_e
      !       il= i+Spn_22n0-Grd_extension-2
      !       nix = dble(il)/dble(ni_trunc)
      !       njx = dble(jl)/dble(nj_trunc)
      !       nkx = sqrt(nix*nix + njx*njx)
      !       if ( nkx > nk_cut ) then
      !          Spn_flt(i,j)= 0.0
      !       else if ( nkx > 1.0 ) then
      !          Spn_flt(i,j)=(cos( (pi_8/2.0)* ((nkx-1.)/(nk_cut-1.))))**2
      !       else
      !          Spn_flt(i,j)= 1.0
      !       end if
      !    end do
      ! end do
      
 1000 format(/' Spn_calfiltre, Large,Small cutoff_scales= ',2f7.2)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_calfiltre


