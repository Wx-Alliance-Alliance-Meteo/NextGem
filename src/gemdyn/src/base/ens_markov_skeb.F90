!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                          Environnement Canada
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

!**   s/r ens_marfield_skeb define a markov chain field for SKEB

      subroutine ens_markov_skeb()
      use cstv
      use ens_gmm_dim
      use ens_gmm_var
      use ens_options
      use ens_param
      use mem_tstp
      use glb_ld
      use rmn_gmm
      use HORgrid_options
      use gem_options
      use init_options
      use lun
      use path
      use gem_fft
      use ptopo
      use out_mod
      use tdpack
      use step_options
      use omp_lib
      use, intrinsic :: iso_fortran_env
      implicit none

!      real, dimension(l_ni,l_nj), intent(out) :: fgem
!
#include <rmnlib_basics.hf>

!
! nlat, nlon                 dimension of the Gaussian grid
! idum                       Semence du générateur de nombres aléatoires
! dt   Pas de temps du modèle (secondes)
! tau  Temps de décorrélation du champ aléatoire f(i,j) (secondes)
! epsi  EXP(-dt/tau/2.146)
! paidum   pointer vers l'etat du generateur sauvegarde idum
!
      character(len=1024) fn0, fn1
      integer :: unf0, unf1, err, errop, gmmstat
      integer :: sig, l ,m, n, i, j, dim,  indx, ier, grd
      integer, pointer :: paiv,paiy,paiset,pagset,paidum
      real    :: fstd, tau
      real,    external :: gasdev
      real(kind=REAL64)  :: pl
      real(kind=REAL64)  :: epsi, pri_8, fact1 

!
!---------------------------------------------------------------------
!

!  Begin Markov chains

! Save random numbers and coefficient ar,ai,br,bi
!!$omp single
      if (write_markov_l) then
         if (ptopo_couleur == 0  .and. ptopo_myproc == 0) then
            fn1=trim(Out_dirname_S)//'/'// 'MRKV_SKEB.bin'
            unf1=0
            open ( unf1,file=trim(fn1),status='NEW', &
                 form='unformatted',iostat=errop )
            if ( errop == 0 ) then
               write(output_unit,1000) 'WRITING', trim(fn1)
               do i=1,36
                  write (unf1)dumdum(i,1)
               end do
               do l=lmin,lmax
                  do m=1,l+1
                     write(unf1)ar_s(lmax-l+1,m)
                     write(unf1)ai_s(lmax-l+1,m)
                     write(unf1)br_s(lmax-l+1,m)
                     write(unf1)bi_s(lmax-l+1,m)
                  end do
               end do
               close(unf1)
            else
               write(output_unit,2000) 'WRITING', trim(fn1)
            endif
         endif
      endif
!!$omp end single 

!!$omp single
      lmin = Ens_skeb_trnl
      lmax = Ens_skeb_trnh
      nlat = Ens_skeb_nlat
      nlon = Ens_skeb_nlon
      fstd = Ens_skeb_std

      tau  = Ens_skeb_tau/2.146
      epsi  = exp(-dt/tau)

      ! Random generator function
      paiv  => dumdum(1,1)
      paiy  => dumdum(33,1)
      pagset=> dumdum(34,1)
      paiset=> dumdum(35,1)
      paidum=> dumdum(36,1)

      do l=lmin,lmax
         fact1=fstd*SQRT(4.*pi_8/real((2*l+1)) &
                   *pspec2(l))*SQRT((1.-epsi*epsi))
         br_s(lmax-l+1,1) = epsi*br_s(lmax-l+1,1)  &
                          + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1
         ar_s(lmax-l+1,1) = epsi*ar_s(lmax-l+1,1)  + br_s(lmax-l+1,1)*fact3
         do m=2,l+1
            br_s(lmax-l+1,m) = epsi*br_s(lmax-l+1,m) &
                             + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1/SQRT(2.)
            ar_s(lmax-l+1,m) = epsi*ar_s(lmax-l+1,m)+br_s(lmax-l+1,m)*fact3
            bi_s(lmax-l+1,m) = epsi*bi_s(lmax-l+1,m) &
                             + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1/SQRT(2.)
            ai_s(lmax-l+1,m) = epsi*ai_s(lmax-l+1,m)+bi_s(lmax-l+1,m)*fact3
         end do
      end do
!!$omp end single 

!!$omp do
      do m=1,lmax+1
         do j=1,nlat
            cc2(1,j,m)=0.d0
            cc2(2,j,m)=0.d0
            cc2(1,j,m)=cc2(1,j,m) + Dot_product(pls(j,1:lmax-lmin+1,m),ar_s(1:lmax-lmin+1,m))
            cc2(2,j,m)=cc2(2,j,m) + Dot_product(pls(j,1:lmax-lmin+1,m),ai_s(1:lmax-lmin+1,m))
         end do
      end do
!!$omp end do

! Fourier Transform (inverse)

      n=-1
      wrk2=0.
!!$omp single 
      do i=1,nlat
         do j=1,lmax+1
            n = n + 2
            wrk2(n)   = cc2(1,i,j)
            wrk2(n+1) = cc2(2,i,j)
         end do
         n=n+nlon-2*lmax
      end do
!!$omp end single

      Fft_type_S = 'PERIODIC'
      Fft_n      = nlon
!!$omp single 
      call itf_fft_drv(wrk2,1,nlon+2,nlat,1)
!!$omp end single

      n=0
!!$omp single
      do j=1,nlat
         do i=1,nlon+2
            n = n + 1
             if(i<=nlon)  f2(i,j) = wrk2(n)
         end do
      end do
!!$omp end single

      ! Interpolation to the processors grids
!!$omp single 
      grd = ezqkdef(nlon,nlat,'A', 0,0,0,0,0)
      ier = ezdefset(Grd_local_gid, grd)
      ier = ezsetopt('INTERP_DEGREE', 'LINEAR')
      ier = ezsint(markov2,f2(1:nlon,1:nlat))
!!$omp end single

      if(Ens_stat)then
!!$omp single 
         call glbstat (markov2,'MCSK','',&
           1,l_ni,1,l_nj,1,1,1,G_ni,1,G_nj,1,1)
!!$omp end single
      end if


 1000 format (/' MARKOV: ',a,' FILE ',a)

 2000 format(' S/R  ENS_MARKOV_SKEB : problem with registering MRKV_PARAM_SKEB file)', &
              /,'======================================================')


      end subroutine ens_markov_skeb

