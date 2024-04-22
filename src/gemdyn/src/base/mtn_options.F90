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
module mtn_options
   use dyn_fisl_options
   use HORgrid_options
   use VERgrid_options
   use tdpack
   use lun
   use dcst
   use glb_ld
   use cstv
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !#
   integer :: mtn_ni = 401
   namelist /mtn_cfgs/ mtn_ni
   !#
   integer :: mtn_nj = 1
   namelist /mtn_cfgs/ mtn_nj
   !#
   integer :: mtn_nk = 65
   namelist /mtn_cfgs/ mtn_nk
   !#
   real :: mtn_dx = 500.
   namelist /mtn_cfgs/ mtn_dx
   !#
   real :: mtn_dz = 300.
   namelist /mtn_cfgs/ mtn_dz
   !#
   real(kind=REAL64) :: mtn_tzero = 303.15
   namelist /mtn_cfgs/ mtn_tzero
   !#
   real :: mtn_flo = 10.
   namelist /mtn_cfgs/ mtn_flo
   !#
   real :: mtn_hwx = 10.
   namelist /mtn_cfgs/ mtn_hwx
   !#
   real :: mtn_hwx1 = 8.
   namelist /mtn_cfgs/ mtn_hwx1
   !#
   real :: mtn_hght = 250.
   namelist /mtn_cfgs/ mtn_hght
   !#
   real :: mtn_nstar = 0.01
   namelist /mtn_cfgs/ mtn_nstar
   !#
   real :: mtn_zblen_thk = 0.
   namelist /mtn_cfgs/ mtn_zblen_thk
   !#
   real :: mtn_wind_seed = 0.
   namelist /mtn_cfgs/ mtn_wind_seed
   !#
   real, dimension(4):: mtn_rcoef = [ 1., 1., -1., -1. ]
   namelist /mtn_cfgs/ mtn_rcoef
   !#
   integer, dimension(2) :: mtn_pos_seed = [ 1, 1 ]
   namelist /mtn_cfgs/ mtn_pos_seed

contains

      integer function mtn_nml (F_unf)
      use, intrinsic :: iso_fortran_env
      implicit none

      integer F_unf

      logical nml_must
      character(len=64) :: nml_S
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         mtn_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=mtn_cfgs)
         return
      end if

      mtn_nml= -1 ; nml_must= .true. ; nml_S= 'mtn_cfgs'

      rewind(F_unf)
      read (F_unf, nml=mtn_cfgs, end= 1001, err=1003)
      mtn_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         mtn_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (mtn_nml < 0 ) return
      if ((Lun_out>=0).and.(mtn_nml==0)) write (Lun_out, 6004) trim(nml_S)
      ! establish horizontal grid configuration
      ! (must absolutely be done here)
      Grd_typ_S='LU'
      Grd_ni = mtn_ni ; Grd_nj = mtn_nj
      Grd_dx = (mtn_dx/Dcst_rayt_8)*(180./pi_8)
      Grd_dy = Grd_dx
      Grd_latr = 0.
      Grd_lonr = (mtn_ni/2 + 20) * Grd_dx
      Grd_maxcfl = 3

      mtn_nml=0

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      return
      end function mtn_nml
!
!-------------------------------------------------------------------
!
      integer function mtn_cfg()
      use, intrinsic :: iso_fortran_env
      implicit none

      integer k
      real(kind=REAL64) c1_8,Exner_8,height_8,pres_8,ptop_8,pref_8,htop_8
!
!     ---------------------------------------------------------------
!
      mtn_cfg = -1

      ! establish vertical grid configuration
      pref_8 = 1.d5

      G_nk = mtn_nk
      htop_8 = (mtn_nk+1)*mtn_dz
      Hyb_rcoef = mtn_rcoef

      if (hyb_H(1) < 0 ) then
         do k=1,G_nk
            height_8=htop_8*(1.d0-(dble(k)-.5d0)/G_nk)
            hyb_H(k)=height_8
         end do
      else
         do k=1024,1,-1
            if(hyb_H(k) < 0 ) G_nk=k-1
         end do
      end if
      
      mtn_cfg = 1

      return
      end function mtn_cfg
!
!     ---------------------------------------------------------------
!
      subroutine mtn_data ( F_u, F_v, F_t, F_s, F_q, F_topo, F_topo_ls, F_orols,&
                            Mminx,Mmaxx,Mminy,Mmaxy,Nk,F_theocase_S )
      use gmm_vt1
      use gmm_geof
      use geomh
      use gem_options
      use glb_pil
      use glb_ld
      use ptopo
      use type_mod
      use ver
      use mem_nest
      use gmm_pw
      use rmn_gmm
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=*) F_theocase_S
      integer Mminx,Mmaxx,Mminy,Mmaxy,Nk
      real F_u    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_v    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_t    (Mminx:Mmaxx,Mminy:Mmaxy,Nk), &
           F_s    (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_topo (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_topo_ls(Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_orols(Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q    (Mminx:Mmaxx,Mminy:Mmaxy,Nk+1)

      integer :: i,j,k,i00
      real    :: a00, a01, a02, xcntr, zdi, zfac, zfac1, capc1, psurf
      real    :: hauteur, press, theta, tempo, dx, slp, slpmax, exner
      real(kind=REAL64)  :: temp1, temp2, oneoRT
      real(kind=REAL64), parameter :: one=1.d0
      real, allocatable, dimension(:,:) :: log_pstar, hm
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,9000)

      allocate (log_pstar(l_minx:l_maxx,l_miny:l_maxy), &
                        hm(l_minx:l_maxx,l_miny:l_maxy) )

!---------------------------------------------------------------------
!     Initialize orography
!---------------------------------------------------------------------

      xcntr = int(float(Grd_ni-1)*0.5)+1
      do j=1-G_haloy,l_nj+G_haloy
         do i=1-G_halox,l_ni+G_halox
            i00 = i + l_i0 - 1
            zdi  = float(i00)-xcntr
            zfac = (zdi/mtn_hwx)**2
            if (      F_theocase_S == 'MTN_SCHAR' &
                .or.  F_theocase_S == 'MTN_SCHAR2' ) then
            zfac1= pi_8 * zdi / mtn_hwx1
            topo_high(i,j,1)= mtn_hght* exp(-zfac) * cos(zfac1)**2
         else if ( F_theocase_S == 'NOFLOW' ) then
            topo_high(i,j,1)= mtn_hght* exp(-zfac)
         else
            topo_high(i,j,1)= mtn_hght/(zfac + 1.)
         end if
      end do
      end do
      call mc2_topols (topo_high(l_minx,l_miny,2),&
                       topo_high(l_minx,l_miny,1),&
                       l_minx,l_maxx,l_miny,l_maxy,Schm_orols_np)

      topo_high= topo_high * grav_8
      if (Vtopo_L) then
         topo_low= 0.
      else
         topo_low = topo_high
      endif

      Nest_current_dayfrac_8= 0.d0
      call var_topo ( F_topo, F_topo_LS, &
                      l_minx,l_maxx,l_miny,l_maxy )      
      F_topo_ls= 0.
      F_topo   =  F_topo    / grav_8
      
         if (     F_theocase_S == 'MTN_SCHAR2' &
             .or. F_theocase_S == 'MTN_PINTY'  &
             .or. F_theocase_S == 'MTN_PINTY2' &
             .or. F_theocase_S == 'NOFLOW' ) then
!
!---------------------------------------
!        Initialize temperature = const.
!---------------------------------------
!
             F_t=mtn_tzero
!
!-----------------------------------------------------
!        Initialize pressure variable  q=RTstar*log(p/pstar)
!-----------------------------------------------------
!
             F_q(:,:,G_nk+1)=-grav_8*(Cstv_tstr_8-mtn_tzero)/mtn_tzero*F_topo(:,:)
             do k=G_nk,1,-1
                F_q(:,:,k)=F_q(:,:,k+1)+grav_8*(Cstv_tstr_8-mtn_tzero)/mtn_tzero &
                *(Ver_dz_8%t(k)+(Ver_b_8%m(k+1)-Ver_b_8%m(k))*F_topo(:,:) +&
                (Ver_c_8%m(k+1)-Ver_c_8%m(k))*F_topo_ls(:,:) )
                !F_q(:,:,k)=F_q(:,:,G_nk+1)
         !       print*,k,'F_q==',F_q(1,1,k)
             end do
         !     stop
          elseif (     F_theocase_S == 'MTN_SCHAR'   &
               .or. F_theocase_S == 'MTN_PINTYNL' ) then
!
!-----------------------------------------------------------------------------
!        Initialize temperature and pressure variable assuming mtn_nstar=const.
!-----------------------------------------------------------------------------
!
               a00 = mtn_nstar**2/grav_8
               capc1 = grav_8**2/(mtn_nstar**2*cpd_8*mtn_tzero)

        !SURFACE
               k=G_nk+1
               do j=1-G_haloy,l_nj+G_haloy
                  do i=1-G_halox,l_ni+G_halox
                     hauteur=F_topo(i,j)
                     log_pstar(i,j)=log(1.d5)+grav_8*hauteur/(rgasd_8*Cstv_Tstr_8)
                     hm(i,j)=hauteur
                     exner=1.d0-capc1*(1.d0-exp(-a00*hauteur))
                     press=1.d5*exner**(1.d0/cappa_8)
                     F_q(i,j,k)=rgasd_8*Cstv_Tstr_8*(log(press)-log_pstar(i,j))
                  end do
               end do
               
               do k=G_nk,1,-1
                  do j=1-G_haloy,l_nj+G_haloy
                     do i=1-G_halox,l_ni+G_halox
                        hauteur=Ver_z_8%m(k)+Ver_b_8%m(k)*F_topo(i,j)+Ver_c_8%m(k)*F_topo_ls(i,j)
                        log_pstar(i,j)=log_pstar(i,j)+grav_8*(hm(i,j)-hauteur)/(rgasd_8*Cstv_Tstr_8)
                        hm(i,j)=hauteur
                        
                        exner=1.d0-capc1*(1.d0-exp(-a00*hauteur))
                        press=1.d5*exner**(1.d0/cappa_8)
                        F_q(i,j,k)=rgasd_8*Cstv_Tstr_8*(log(press)-log_pstar(i,j))
                        
                        hauteur=Ver_z_8%t(k)+Ver_b_8%t(k)*F_topo(i,j)+Ver_c_8%t(k)*F_topo_ls(i,j)
                        exner=1.d0-capc1*(1.d0-exp(-a00*hauteur))
                        theta=mtn_tzero*exp(a00*hauteur)
                        F_t(i,j,k)=theta*exner
                     end do
                  end do
               end do
               
            else
               print*,'error in mtn_case'
               stop
            end if

      if ( F_theocase_S == 'NOFLOW' ) then
!        calculate maximum mountain slope
         slpmax=0
         dx=Dcst_rayt_8*Grd_dx*pi_8/180.
         do j=1-G_haloy,l_nj+G_haloy
            do i=1-G_halox,l_ni+G_halox
               slp=abs(F_topo(i,j)-F_topo(i-1,j))/dx
               slpmax=max(slp,slpmax)
            end do
         end do
         slpmax=(180.d0/pi_8)*atan(slpmax)
         print*,"SLPMAX=",slpmax," DEGREES"
      end if

!-----------------------------------------------------------------------
!     Transform orography from geometric to geopotential height
!-----------------------------------------------------------------------
      F_topo    = grav_8 * F_topo
      F_topo_ls = 0.
!
!---------------------------------------------------------------------
!     Set winds (u,v)
!---------------------------------------------------------------------
!
      F_u(:,:,1:G_nk)  = mtn_flo
      F_v(:,:,1:G_nk)  = 0.0
!
!-------------------------------------------------
!        Introduce an horizontal wind perturbation
!-------------------------------------------------
!
      if ( F_theocase_S == 'NOFLOW' ) then
         do i=1-G_halox,l_ni+G_halox
            if(i+l_i0-1==mtn_pos_seed(1)) F_u(i,:,mtn_pos_seed(2))=mtn_wind_seed
         end do
      end if

      deallocate(log_pstar,hm)
!
 9000 format(/,'CREATING INPUT DATA FOR MOUNTAIN WAVE THEORETICAL CASE' &
            /,'======================================================')
!
!     -----------------------------------------------------------------
!
      return
      end subroutine mtn_data

end module mtn_options
