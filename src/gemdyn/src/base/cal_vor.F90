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

!** s/r cal_vor - Computes horizontal vorticity and relative vorticity

      subroutine cal_vor ( F_QR,F_QQ, F_uu,F_vv          , &
                           F_filtqq, F_coefqq, F_absvor_L, &
                           Minx,Maxx,Miny,Maxy,Nk )
      use dcst
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use gmm_vt0
      use gmm_geof
      use glb_ld
      use tdpack
      use ptopo
      use lun
      implicit none


      logical  F_absvor_L
      integer  F_filtqq, Minx,Maxx,Miny,Maxy,Nk
      real     F_QR (Minx:Maxx,Miny:Maxy,Nk), &
               F_QQ (Minx:Maxx,Miny:Maxy,Nk), &
               F_uu (Minx:Maxx,Miny:Maxy,Nk), &
               F_vv (Minx:Maxx,Miny:Maxy,Nk), F_coefqq

      integer i, j, k, i0, in, j0, jn
      include 'mpif.h'
      integer :: err,comm
      real(kind=REAL64) sum_8, summa_8, sumet_8, sumen_8
      real(kind=REAL64) :: gathS(Ptopo_numproc*Ptopo_ncolors)
      real deg2rad
!
!----------------------------------------------------------------------
!
      i0 = 1
      in = l_niu
      j0 = 1
      jn = l_njv

      do k = 1 , Nk
         do j = j0, jn
            do i = i0, in
               F_QR(i,j,k) = ((F_vv(i+1,j,k) - F_vv(i,j,k)) * geomh_invDXv_8(j)) &
                           - ( (F_uu(i,j+1,k)*geomh_cy_8(j+1) - F_uu(i,j  ,k)*geomh_cy_8(j  )) &
                           * geomh_invDY_8 * geomh_invcyv_8(j))
            end do
         end do
         F_QR(1:i0-1,:,k) = 0. ; F_QR(in+1:l_ni,:,k)= 0.
         F_QR(:,1:j0-1,k) = 0. ; F_QR(:,jn+1:l_nj,k)= 0.
      end do

      if (F_filtqq > 0) call filter ( F_QR, F_filtqq,F_coefqq, &
                                  l_minx,l_maxx,l_miny,l_maxy,Nk )

      if (F_absvor_L)then
         deg2rad= pi_8/180.d0
         do k =  1, Nk
            do j = j0, jn
               do i = i0, in
                  F_QQ(i,j,k)= F_QR(i,j,k) + 2.0*Dcst_omega_8 &
                             * sin(geomh_latrx(i,j)*deg2rad)
               end do
            end do
            F_QQ(1:i0-1,:,k) = 0. ; F_QQ(in+1:l_ni,:,k)= 0.
            F_QQ(:,1:j0-1,k) = 0. ; F_QQ(:,jn+1:l_nj,k)= 0.
         end do
      end if

      comm = COMM_multigrid

!potential enstrophy
      sum_8=0.
      do k = 1 , Nk
         do j = j0, jn
            do i = i0, in
            sum_8 = sum_8 + 0.5*F_QQ(i,j,k)**2/(qt0(i,j,k)-fis0(i,j)+1./Cstv_h0inv_8) * geomh_mask_8(i,j)
            end do
         end do
      end do

      call MPI_Allgather(sum_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

      sumen_8 = sum(gathS) 

!total energy 
      sum_8=0.
      do k = 1 , Nk
         do j = j0, jn
            do i = i0, in
            sum_8 = sum_8 + 0.5*((qt0(i,j,k)-fis0(i,j)+1./Cstv_h0inv_8)*(F_uu(i,j,k)**2 + F_vv(i,j,k)**2) + &
                            grav_8*((qt0(i,j,k)+1./Cstv_h0inv_8)**2-fis0(i,j)**2) ) * geomh_mask_8(i,j)
            end do
         end do
      end do

      call MPI_Allgather(sum_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

      sumet_8 = sum(gathS) 
      
!total mass 
      sum_8=0.
      do k = 1 , Nk
         do j = j0, jn
            do i = i0, in
            sum_8 = sum_8 + (qt0(i,j,k)-fis0(i,j)+1./Cstv_h0inv_8)*geomh_mask_8(i,j)
            end do
         end do
      end do

      call MPI_Allgather(sum_8,1,MPI_DOUBLE_PRECISION,gathS,1,MPI_DOUBLE_PRECISION,comm,err)

      summa_8 = sum(gathS) 

      if(Ptopo_myproc.eq.0) then
         write(Lun_out,1001) 'summa,sumet,sumen',summa_8,sumet_8,sumen_8
      endif

 1001 format(1X,A24,1X,1E19.8,1X,1E19.8,1X,1E19.8)
!
!----------------------------------------------------------------------
!
      return
      end
