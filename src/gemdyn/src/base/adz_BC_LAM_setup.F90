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

!**   s/r adz_BC_LAM_setup - Set F_mask_o/F_mask_i for Flux calculations
!                            based on Aranami et al. (2015)

      subroutine adz_BC_LAM_setup (F_amask_o,F_amask_i,&
                                   F_aminx,F_amaxx,F_aminy,F_amaxy,F_nk)
      use adz_mem
      implicit none

      !arguments
      !---------
      integer, intent(in) :: F_aminx,F_amaxx,F_aminy,F_amaxy,F_nk
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,F_nk), &
            intent(out) :: F_amask_o,F_amask_i

      integer :: i,j,k
!
!-------------------------------------------------------------------
!
      !For Flux_out: Set mixing ratio = 0 on NEST-HV
      !---------------------------------------------
      F_amask_o = 0.
      do k=Adz_k0t,F_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               F_amask_o(i,j,k) = 1.
            end do
         end do
      end do

      !For Flux_in: Set mixing ratio = 0 on CORE-HV
      !--------------------------------------------
      F_amask_i = 1.
      do k=Adz_k0t,F_nk
         do j=1+pil_s,l_nj-pil_n
            do i=1+pil_w,l_ni-pil_e
               F_amask_i(i,j,k) = 0.
            end do
         end do
      end do

      call gem_xch_halo ( F_amask_o,F_aminx,F_amaxx,F_aminy,F_amaxy,&
                          F_nk,-1 )
      call gem_xch_halo ( F_amask_i,F_aminx,F_amaxx,F_aminy,F_amaxy,&
                          F_nk,-1 )
!
!-------------------------------------------------------------------
!
      return
      end subroutine adz_BC_LAM_setup
