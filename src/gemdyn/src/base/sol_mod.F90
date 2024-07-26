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

module sol_mem
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   integer :: Sol_pil_w,Sol_pil_e,Sol_pil_n,Sol_pil_s
   integer :: Sol_niloc,Sol_njloc,Sol_nloc,Sol_nk,Sol_sock_nk
   integer :: Sol_i0,Sol_in,Sol_j0,Sol_jn,Sol_dimx,Sol_dimy
   integer :: Sol_miny,Sol_maxy,Sol_mink,Sol_maxk,Sol_k0
   integer :: Sol_ii0,Sol_iin,Sol_jj0,Sol_jjn,Sol_imin,Sol_imax,Sol_jmin,Sol_jmax

   real(kind=REAL64), pointer, dimension (:,:,:) :: Sol_rhs, Sol_lhs

   real(kind=REAL64) :: isol_i, isol_d
   real(kind=REAL64) :: norm_residual, relative_tolerance, nu , r0, rr2, ro2, lcl_sum(2)
   real(kind=REAL64),dimension(:      ), allocatable :: gg, rot_cos, rot_sin, IPIV_arr
   real(kind=REAL64),dimension(:,:    ), allocatable :: v_lcl_sum,rr,tt, hessenberg,thread_s
   real(kind=REAL64),dimension(:,:,:  ), allocatable :: work_space,thread_s2,fdg,w2_8,w3_8
   real(kind=REAL128),dimension(:,:   ), allocatable :: thread_s128
   real(kind=REAL64),dimension(:,:,:,:), contiguous, pointer :: vv, wint_8
   real(kind=REAL64),dimension(:), allocatable :: m_west,m_east,m_south,m_north
   real             ,dimension (:,:,: ), allocatable :: fdg2,ext_q

end module sol_mem
