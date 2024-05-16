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

!**s/r eqspng - apply vertical diffusion at and near model

      subroutine eqspng ()
      use hvdif_options
      use glb_ld
      use gmm_vt1
      use lun
      implicit none

   !    There are nlev levels that will be diffused.
   !    There are nlev-1 coefficent passed in namelist by user.
   !
   !    The boundary conditions are flux=0. This is achieved by
   !    putting coef=0 at boundaries (see drawing below).
   !    Also the index km and kp make the wind derivatives zero
   !    at boundary.
   !
   !    ---- this u equals u(1) (see km in loop)   \
   !                                               |
   !    ==== coef(1)=0, Top boundary.              | this derivative = 0
   !                                               |
   !    ---- u(1) first levels diffused            /
   !
   !    ==== coef(2)=first coef passed by user
   !
   !    ---- u(2)
   !
   !        ...
   !
   !    ---- u(nlev-1)
   !
   !    ==== coef(nlev)=last coef passed by user
   !
   !    ---- u(nlev) last level diffused           \
   !                                               |
   !    ==== coef(nlev+1)=0, Bottom boundary.      | this derivative = 0
   !                                               |
   !    ---- this u equal u(nlev) (see kp in loop) /
   !
   !_____________________________________________________________________

      real, dimension(l_ni,eq_nlev) :: u,v,w,zd
      integer :: i,j,k,km,kp
!
!-------------------------------------------------------------------
!
!!$omp do
   do j=1,l_nj
      do k=1,eq_nlev
         kp=min(eq_nlev,k+1)
         km=max(1,k-1)
         do i=1,l_ni
            u(i,k)=ut1(i,j,k)+eponmod(i,j)*(cp(k)*(ut1(i,j,kp)-ut1(i,j,k )) &
                                           -cm(k)*(ut1(i,j,k )-ut1(i,j,km)))
            v(i,k)=vt1(i,j,k)+eponmod(i,j)*(cp(k)*(vt1(i,j,kp)-vt1(i,j,k )) &
                                           -cm(k)*(vt1(i,j,k )-vt1(i,j,km)))
            w(i,k)=wt1(i,j,k)+eponmod(i,j)*(cp(k)*(wt1(i,j,kp)-wt1(i,j,k )) &
                                           -cm(k)*(wt1(i,j,k )-wt1(i,j,km)))
            zd(i,k)=zdt1(i,j,k)+eponmod(i,j)*(cp(k)*(zdt1(i,j,kp)-zdt1(i,j,k )) &
                                           -cm(k)*(zdt1(i,j,k )-zdt1(i,j,km)))
         end do
      end do
      ut1(1:l_ni,j,1:eq_nlev)=u(1:l_ni,1:eq_nlev)
      vt1(1:l_ni,j,1:eq_nlev)=v(1:l_ni,1:eq_nlev)
      wt1(1:l_ni,j,1:eq_nlev)=w(1:l_ni,1:eq_nlev)
      zdt1(1:l_ni,j,1:eq_nlev)=zd(1:l_ni,1:eq_nlev)
   end do
!!$omp end do
!
!-------------------------------------------------------------------
!
      return
      end subroutine eqspng
