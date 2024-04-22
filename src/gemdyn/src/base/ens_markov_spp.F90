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

!**   s/r ens_markov_spp - define a markov chain field for PTP
!
      subroutine ens_markov_spp()
      use ens_param
      use phy_itf
      use ens_gmm_dim
      use step_options
      use ens_gmm_var
      use ens_options
      use HORgrid_options
      use gem_options
      use init_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use rmn_gmm
      use path
      use gem_fft
      use ptopo
      use out_mod
      use geomh, only: geomh_latrx
      use gmm_phy, only: phy_cplm, phy_cplt
      use ens_spp, only: spp_list, spp_ncha, spp_lmax, spp_mmax
      use clib_itf_mod, only: clib_toupper
      use wb_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none
!
!author: Rabah Aider R.P.N.-A
!
#include <rmnlib_basics.hf>
!

!
! nlat, nlon                 dimension of the Gaussian grid
! idum                       Semence du generateur de nombres al?atoires
! paidum   pointer vers l'etat du generateur sauvegarde idum
! dt   Pas de temps du mod?le (secondes)
! tau  Temps de decorrelation du champ aleatoire f(i,j) (secondes)

      character(len=1024) :: fn0, fn1
      integer, pointer :: paiv,paiy,paiset,pagset,paidum
      integer :: l ,m, n, nc,np, i, j, indx, ier, istat, grd
      integer ::  err, errop, ierr ,unf0, unf1, nch2d, spp_indx, stat,gmmstat 

      real    :: fstd,  fstr, sumsp , fact
      real,    external ::  gasdev
      real,    dimension(:,:,:),pointer   ::  ptr3d

      real(kind=REAL64)  :: fact1,  pri_8
      real(kind=REAL64)   ::  fmax, fmin , fmean
!
!-------------------------------------------------------------------
!
      nch2d = Ens_ptp_ncha + spp_ncha
      if (nch2d == 0) return

! Begin Markov chains

      CONSTRUCT_CHAINS: do nc=1,nch2d
!!$omp single
         lmin = vlmin(nc)
         lmax = vlmax(nc)
         nlon = vnlon(nc)
         nlat = vnlat(nc)
         fmin = vfmin(nc)
         fmax = vfmax(nc)
         fstd = vfstd(nc)
         fstr = vfstr(nc)
         eps = veps(nc)
         np = nc+1

! Random generator function
         paiv  => dumdum(1,np)
         paiy  => dumdum(33,np)
         paiset=> dumdum(34,np)
         pagset=> dumdum(35,np)
         paidum=> dumdum(36,np)


! Save random numbers and coefficient ar,ai,br,bi
         if (write_markov_l) then
            if (ptopo_couleur == 0  .and. ptopo_myproc == 0) then
               fn1=trim(Out_dirname_S)//'/'// 'MRKV_SPP_'//trim(vname(nc))//'.bin'
               unf1=0
               open ( unf1,file=trim(fn1),status='NEW', &
                     form='unformatted',iostat=errop )
               if ( errop == 0 ) then
                  write(output_unit,1000) 'WRITING', trim(fn1)
                  do i=1,36
                     write (unf1)dumdum(i,np)
                  end do
                  do l=lmin,lmax
                     do m=1,l+1
                        write(unf1)ar_p(lmax-l+1,m,nc)
                        write(unf1)ai_p(lmax-l+1,m,nc)
                        write(unf1)br_p(lmax-l+1,m,nc)
                        write(unf1)bi_p(lmax-l+1,m,nc)
                     end do
                  end do
                  close(unf1)
               else
                  write(output_unit,2000) 'WRITING', trim(fn1)
               endif
            endif
         endif

         do l=lmin,lmax
            fact1=fstd*SQRT(4.*pi_8/real((2*l+1))*pspec1(l,nc))*SQRT((1.-eps*eps))
            br_p(lmax-l+1,1,nc) = eps*br_p(lmax-l+1,1,nc)  &
                                + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1
            ar_p(lmax-l+1,1,nc) = eps*ar_p(lmax-l+1,1,nc)  + br_p(lmax-l+1,1,nc)*fact2(nc)
            do m=2,l+1
               br_p(lmax-l+1,m,nc) = eps*br_p(lmax-l+1,m,nc) &
                                   + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1/SQRT(2.)
               ar_p(lmax-l+1,m,nc) = eps*ar_p(lmax-l+1,m,nc)+br_p(lmax-l+1,m,nc)*fact2(nc)
               bi_p(lmax-l+1,m,nc) = eps*bi_p(lmax-l+1,m,nc) &
                                   + gasdev(paiv,paiy,paiset,pagset,paidum)*fact1/SQRT(2.)
               ai_p(lmax-l+1,m,nc) = eps*ai_p(lmax-l+1,m,nc)+bi_p(lmax-l+1,m,nc)*fact2(nc)
            end do
         end do
!!$omp end single 

!!$omp do
         do m=1,lmax+1
            do j=1,nlat
	       cc1(1,j,m,nc)= Dot_product(plp(j,1:lmax-lmin+1,m,nc),ar_p(1:lmax-lmin+1,m,nc))
               cc1(2,j,m,nc)= Dot_product(plp(j,1:lmax-lmin+1,m,nc),ai_p(1:lmax-lmin+1,m,nc))
            end do
         end do
!!$omp end do

!  Fourier Transform (inverse)

!!$omp single 
         n=-1
         wrk1=0.
         do i=1,nlat
            do j=1,lmax+1
               n = n + 2
               wrk1(n,nc)   = cc1(1,i,j,nc)
               wrk1(n+1,nc) = cc1(2,i,j,nc)
            end do
            n=n+nlon-2*lmax
         end do
!!$omp end single

         Fft_type_S = 'PERIODIC'
         Fft_n      = nlon
!!$omp single 
         call itf_fft_drv(wrk1(:,nc),1,nlon+2,nlat,1)
!!$omp end single

         n=0
!!$omp single
         do j=1,nlat
            do i=1,nlon+2
               n = n + 1
               if (i <= nlon) f1(i,j,nc) = wrk1(n,nc)
            end do
         end do
!!$omp end single

!  Interpolation to the processors grids and fill in perbus

!!$omp single 
         grd = ezqkdef(nlon,nlat,'A', 0,0,0,0,0)
         ier = ezdefset(Grd_local_gid, grd)
         ier = ezsetopt('INTERP_DEGREE', 'LINEAR')
!!$omp end single

!  Check the limits, stretch, and add mean if stretching asked
!  for the physics perturbation

         if(fstr /= 0.0)then
            fmean=(fmin+fmax)/2.
!!$omp do 
            do j=1,nlat
            do i=1,nlon          
 	           f1_str(i,j,nc)=ERF(f1(i,j,nc)/(fstr*fstd)/SQRT(2.)) *(fmax-fmin)/2. + fmean
               enddo
            enddo
!!$omp enddo 
 
!!$omp single 
            ier = ezsint(markov1(:,:,nc),f1_str(:,:,nc))
!!$omp end single
         else
!!$omp do 
            do j=1,l_nj
               do i=1,l_ni          
                  markov1(i,j,nc)=1.0
               enddo
            enddo
!!$omp enddo 
         end if

         if(Ens_stat)then
!!$omp single
            call glbstat (markov1(:,:,nc),'MSPP',vname(nc),&
            1,l_ni,1,l_nj,1,1,1,G_ni,1,G_nj,1,1)
!!$omp end single
         end if

         ! Dynamics perturbations are intercepted here
         if (vname(nc) == 'ADV_RHSINT') then 
!!$omp do
            do j=1,l_nj
               do i=1,l_ni
                  mcrhsint(i,j) = tropwt(i,j) * markov1(i,j,nc)
               enddo
            enddo
!!$omp enddo 
         endif
         if (vname(nc) == 'PHYCPL') then
!!$omp do
            do j=1,l_nj
               do i=1,l_ni
            	  phy_cplm(i,j) = markov1(i,j,nc)
                  phy_cplt(i,j) = markov1(i,j,nc)
               enddo
            enddo
!!$omp enddo 
         endif

      end do CONSTRUCT_CHAINS

      ptr3d => markov1(Grd_lphy_i0:Grd_lphy_in, &
                        Grd_lphy_j0:Grd_lphy_jn, 1:nch2d)
      istat = phy_put(ptr3d,'mrk2',F_npath='V',F_bpath='P')

1000 format (/' MARKOV: ',a,' FILE ',a)

2000 format(' S/R  ENS_MARFIELD_PTP : problem in WRITING  MRKV_PARAM_PTP_SPP file)', &
              /,'======================================================')

end subroutine ens_markov_spp
