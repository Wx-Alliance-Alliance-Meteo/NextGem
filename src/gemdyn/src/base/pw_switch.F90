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

!**s/r pw_switch -  Rename time level M -> P
!
      subroutine pw_switch ()
      use rmn_gmm
      use gmm_contiguous
      use gmm_pw
      implicit none

      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_total  = ['pwPLUS ','pwMOINS']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_uulist = ['PW_UU:P','PW_UU:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_vvlist = ['PW_VV:P','PW_VV:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_ttlist = ['PW_TT:P','PW_TT:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_pmlist = ['PW_PM:P','PW_PM:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_ptlist = ['PW_PT:P','PW_PT:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_gzlist = ['PW_GZ:P','PW_GZ:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_melist = ['PW_ME:P','PW_ME:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_p0list = ['PW_P0:P','PW_P0:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_pm8list = ['PW_PM8:P','PW_PM8:M']
      character(len=GMM_MAXNAMELENGTH), dimension(2), parameter :: pw_p08list = ['PW_P08:P','PW_P08:M']
      integer :: istat
!
!     ---------------------------------------------------------------
!
!!$omp single
      istat = gmm_shuffle(pw_total )
      istat = gmm_shuffle(pw_uulist)
      istat = gmm_shuffle(pw_vvlist)
      istat = gmm_shuffle(pw_ttlist)
      istat = gmm_shuffle(pw_pmlist)
      istat = gmm_shuffle(pw_ptlist)
      istat = gmm_shuffle(pw_gzlist)
      istat = gmm_shuffle(pw_melist)
      istat = gmm_shuffle(pw_p0list)
      istat = gmm_shuffle(pw_pm8list)
      istat = gmm_shuffle(pw_p08list)

      istat = gmm_get('pwPLUS'  ,pwPLUS )
      istat = gmm_get('pwMOINS' ,pwMOINS)
      istat = gmm_get(gmmk_pw_uu_plus_s  ,pw_uu_plus)
      istat = gmm_get(gmmk_pw_vv_plus_s  ,pw_vv_plus)
      istat = gmm_get(gmmk_pw_tt_plus_s  ,pw_tt_plus)
      istat = gmm_get(gmmk_pw_pm_plus_s  ,pw_pm_plus)
      istat = gmm_get(gmmk_pw_pt_plus_s  ,pw_pt_plus)
      istat = gmm_get(gmmk_pw_gz_plus_s  ,pw_gz_plus)
      istat = gmm_get(gmmk_pw_me_plus_s  ,pw_me_plus)
      istat = gmm_get(gmmk_pw_p0_plus_s  ,pw_p0_plus)
      istat = gmm_get(gmmk_pw_pm_plus_8_s  ,pw_pm_plus_8)
      istat = gmm_get(gmmk_pw_p0_plus_8_s  ,pw_p0_plus_8)

      istat = gmm_get(gmmk_pw_uu_moins_s  ,pw_uu_moins)
      istat = gmm_get(gmmk_pw_vv_moins_s  ,pw_vv_moins)
      istat = gmm_get(gmmk_pw_tt_moins_s  ,pw_tt_moins)
      istat = gmm_get(gmmk_pw_pt_moins_s  ,pw_pt_moins)
      istat = gmm_get(gmmk_pw_gz_moins_s  ,pw_gz_moins)
      istat = gmm_get(gmmk_pw_pm_moins_s  ,pw_pm_moins)
      istat = gmm_get(gmmk_pw_me_moins_s  ,pw_me_moins)
      istat = gmm_get(gmmk_pw_p0_moins_s  ,pw_p0_moins)
      istat = gmm_get(gmmk_pw_pm_moins_8_s ,pw_pm_moins_8)
      istat = gmm_get(gmmk_pw_p0_moins_8_s ,pw_p0_moins_8)
!!$omp end single
      call adz_wnds_ext ()
!
!     ---------------------------------------------------------------
!
      return
      end subroutine pw_switch
