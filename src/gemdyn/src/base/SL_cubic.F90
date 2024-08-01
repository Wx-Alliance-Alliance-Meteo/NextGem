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

!**s/r SL_cubic

!Arguments for function:
!	1.  F_qut: out array, will be set to src in SL_interp
!	2.  F_lin: out array; not sure of purpose
!	3.  F_min: out array; not sure of purpose
!	4.  F_max: out array ; not sure of purpose
!	5.  F_in:  values used for the interpolation, true q vals
!	6.  F_xyz: values to interpolate to
!	7.  F_num: total points but for the slice???
!	8.  F_minx
!	9.  F_maxx
!	10. F_miny
!	11. F_maxy
!	12. F_ns: dimext in SL_interp = total d.o.f = ni*ni*nk
!	13. F_nptr: 1, not sure what this means
!	14. F_k0: start of k index, i.e. 1
!	15. F_kn: end k index, i.e. l_nk
!	16. F_lev: flag for thermo or mom; i.e. 'm' or 't'


      subroutine SL_cubic ( F_qut, F_lin , F_min , F_max , F_in, F_xyz, &
                            F_num, F_minx, F_maxx, F_miny, F_maxy     , &
                            F_ns , F_nptr, F_k0  , F_kn  , F_lev)
      use adz_mem
      use dynkernel_options
      use ver
      use, intrinsic :: iso_fortran_env

      implicit none

      character(len=1) , intent(in) :: F_lev
      integer, intent(in) :: F_num,F_minx,F_maxx,F_miny,F_maxy,F_ns,F_nptr,F_k0,F_kn
      real,dimension(F_ns,F_nptr) , intent(in ) :: F_in
      real,dimension(F_num,F_nptr), intent(out) :: F_qut,F_lin,F_min,F_max
      real,dimension(3,F_num), intent(in) :: F_xyz

      logical :: zcubic_L,zqutic_L
      integer :: n,o1,indxij,ii1,jj1,kk1,nit,nij,i
      real    :: xyz_1,xyz_2
      real(kind=REAL64), parameter :: cp167 = 1.0/6.0, cm167 = -1.0/6.0, cp5 = 0.5, cm5 = -0.5
      real(kind=REAL64), dimension(:), pointer ::&
                         p_zabcd_8, p_zbacd_8,&
                         p_zcabd_8, p_zdabc_8, posz, dz, odz,&
                         p_zxabcde_8, p_zaxbcde_8, p_zbxacde_8,&
                         p_zcxabde_8, p_zdxabce_8, p_zexabcd_8
      real(kind=REAL64) :: px(0:3),py(0:3),qz(4),del,del123,del210,&
                           ra,rb,rc,rd,rx,re,lx(2),ly(2),lz(2),zt
!
!---------------------------------------------------------------------
!
      if (F_lev == 'm') then
      posz => Ver_z_8%m
      dz   => Adz_delz_m
      odz  => Adz_odelz_m

      p_zabcd_8 => Adz_zabcd_8%m
      p_zbacd_8 => Adz_zbacd_8%m
      p_zcabd_8 => Adz_zcabd_8%m
      p_zdabc_8 => Adz_zdabc_8%m

      p_zxabcde_8 => Adz_zxabcde_8%m
      p_zaxbcde_8 => Adz_zaxbcde_8%m
      p_zbxacde_8 => Adz_zbxacde_8%m
      p_zcxabde_8 => Adz_zcxabde_8%m
      p_zdxabce_8 => Adz_zdxabce_8%m
      p_zexabcd_8 => Adz_zexabcd_8%m
      endif
      if (F_lev == 't') then
      posz => Ver_z_8%t
      dz   => Adz_delz_t
      odz  => Adz_odelz_t

      p_zabcd_8 => Adz_zabcd_8%t
      p_zbacd_8 => Adz_zbacd_8%t
      p_zcabd_8 => Adz_zcabd_8%t
      p_zdabc_8 => Adz_zdabc_8%t

      p_zxabcde_8 => Adz_zxabcde_8%t
      p_zaxbcde_8 => Adz_zaxbcde_8%t
      p_zbxacde_8 => Adz_zbxacde_8%t
      p_zcxabde_8 => Adz_zcxabde_8%t
      p_zdxabce_8 => Adz_zdxabce_8%t
      p_zexabcd_8 => Adz_zexabcd_8%t
      endif
      
      nit = (F_maxx - F_minx + 1)
      nij = (F_maxx - F_minx + 1)*(F_maxy - F_miny + 1)

      !for each point in the slice...
      do n=1,F_num

         !---x direction---
         xyz_1 = max(real(Adz_iminposx)+1,F_xyz(1,n))
         xyz_1 = min(real(Adz_imaxposx)-1,xyz_1)
         ii1   = xyz_1
         del   = xyz_1 - dble(ii1)

         px(0) = del*(del-1.d0)*(del-2.d0)*cm167
         px(1) = (del+1.d0)*(del-1.d0)*(del-2.d0)*cp5
         px(2) = del*(del+1.d0)*(del-2.d0)*cm5
         px(3) = del*(del+1.d0)*(del-1.d0)*cp167

         lx(1) = 1.d0 - del ; lx(2) = del

         !---y direction---
         xyz_2 = max(real(Adz_iminposy)+1,F_xyz(2,n))
         xyz_2 = min(real(Adz_imaxposy)-1,xyz_2)
         jj1   = xyz_2
         del   = xyz_2 - dble(jj1)

         py(0) = del*(del-1.d0)*(del-2.d0)*cm167
         py(1) = (del+1.d0)*(del-1.d0)*(del-2.d0)*cp5
         py(2) = del*(del+1.d0)*(del-2.d0)*cm5
         py(3) = del*(del+1.d0)*(del-1.d0)*cp167

         ly(1) = 1.d0 - del ; ly(2) = del

         ii1 = ii1 - F_minx - 1
         jj1 = jj1 - F_miny - 1
         indxij = (jj1-Adz_joff-1)*nit + (ii1-Adz_ioff) - nit - 1

         !---z direction---
         kk1 = F_xyz(3,n)

         !zqutic_L = (kk1 >= F_k0+3-1) .and.  (kk1 <= F_kn-3)  
         !zcubic_L =((kk1 == F_k0+2-1) .or.   (kk1 == F_kn-2))
  
         !cubic everywhere expect first 2 and last 2 levels
         !first 2 and last 2 use linear, see tricublin
         zcubic_L = (kk1 < F_kn-2) .and. (kk1 > F_k0+2)  

         if (.not.zcubic_L) kk1 = min(F_kn-1,max(F_k0,kk1))

         zt = posz(kk1) + (F_xyz(3,n)-dble(kk1))*dz(kk1)
         lz(2) = (zt-posz(kk1))*odz(kk1); lz(1) = 1.d0 - lz(2)

         if (zcubic_L) then

            o1 = (kk1-2)*nij + indxij

            ra = posz(kk1-1)
            rb = posz(kk1  )
            rc = posz(kk1+1)
            rd = posz(kk1+2)

            qz(1) = triprd(zt,rb,rc,rd)*p_zabcd_8(kk1)
            qz(2) = triprd(zt,ra,rc,rd)*p_zbacd_8(kk1)
            qz(3) = triprd(zt,ra,rb,rd)*p_zcabd_8(kk1)
            qz(4) = triprd(zt,ra,rb,rc)*p_zdabc_8(kk1)

         else

            o1 = (kk1-1)*nij + indxij

            qz(1) = 0.d0
            qz(4) = 0.d0
            qz(3) = (zt-posz(kk1))*odz(kk1)
            qz(2) = 1.d0 - qz(3)

         end if

         do i=1,F_nptr
            call bicubHQV_lin_min_max_zyx (&
               F_qut(n,i),F_lin(n,i),F_min(n,i),F_max(n,i),F_in(o1,i),&
               px,py,qz,lx,ly,lz,nit,nij,zcubic_L)
         end do

      end do
!
!---------------------------------------------------------------------
!
      return

contains

      real(kind=REAL64) function triprd(za,zb,zc,zd)
         real(kind=REAL64), intent(in) :: za,zb,zc,zd
         triprd = ((za-zb)*(za-zc)*(za-zd))
      end function triprd

      real(kind=REAL64) function quiprd(zx,za,zb,zc,zd,ze)
         real(kind=REAL64), intent(in) :: zx,za,zb,zc,zd,ze 
         quiprd = (zx-za)*(zx-zb)*(zx-zc)*(zx-zd)*(zx-ze)
      end function quiprd

      end subroutine SL_cubic
