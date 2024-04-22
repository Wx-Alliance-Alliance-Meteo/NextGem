
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


      subroutine ens_markov_skeb_init
      use cstv
      use dcst
      use geomh
      use ens_gmm_dim
      use ens_gmm_var
      use ens_options
      use ens_param
      use glb_ld
      use rmn_gmm
      use gem_options
      use init_options
      use lun
      use path
      use ptopo
      use out_mod
      use step_options
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

!      real, dimension(l_ni,l_nj), intent(out) :: fgem
!
#include <rmnlib_basics.hf>

       real,    external :: gasdev
       real(kind=REAL64)  :: pl
!
! nlat, nlon                 dimension of the Gaussian grid
! idum                       Semence du générateur de nombres aléatoires
! paidum   pointer vers l'etat du generateur sauvegarde idum
! dt   Pas de temps du modèle (secondes)
! tau  Temps de décorrélation du champ aléatoire f(i,j) (secondes)
!  epsi  EXP(-dt/tau/2.146)

      character(len=1024) fn0, fn1
      logical, save :: init_done=.false.
      integer, pointer :: paiv,paiy,paiset,pagset,paidum
      integer :: sig, l ,m, n, i, j, dim, ier
      integer :: unf0, unf1, err, errop, gmmstat
      real    :: fstd, tau, fact 
      real(kind=REAL64)  :: rad2deg_8, deg2rad_8 
      real(kind=REAL64),    dimension(:), allocatable :: fact1
      real(kind=REAL64)   epsi, sigma, w1,w2,sumsp

!---------------------------------------------------------------------
!
      dt=real(Cstv_dt_8)
      rad2deg_8=180.0d0/pi_8
      itstep_s=step_dt*step_kount
      iperiod_iau = iau_period
      write_markov_l=(itstep_s==iperiod_iau)

      if (.not.init_done) then
        if (Lun_out > 0) then
        write( Lun_out,1000)
      end if
        init_done=.true.
      end if

      lmin = Ens_skeb_trnl
      lmax = Ens_skeb_trnh
      nlat = Ens_skeb_nlat
      nlon = Ens_skeb_nlon

      fstd = Ens_skeb_std
      tau  = Ens_skeb_tau/2.146
      epsi  = exp(-dt/tau)

      allocate (pspec2(lmin:lmax))

! Compute Associated Legendre polynomials to use in each timestep
         pls=0.D0
         do l=lmin,lmax
            fact=DSQRT((2.D0*DBLE(l)+1.D0)/(4.D0*pi_8))
            do m=0,l
               sig=(-1.D0)**(l+m)

               do j=1,nlat/2
                  call pleg (l, m, j, nlat, pl)
                  pls(j,lmax-l+1,m+1)=pl*fact
                  pls(nlat-j+1,lmax-l+1,m+1)=pl*fact*sig
               end do
            end do
         end do

         !Initialise spectral coeffs and stochastic params
         if(Ens_recycle_mc) then
            !Read saved stochastic numbers and spectral coeffs ar,br.ai,bi
            if (ptopo_myproc==0 .and. ptopo_couleur==0) then
                  unf0=1
                  fn0= trim(Path_input_S)//'/MODEL_INPUT'//'/MRKV_SKEB.bin'
                  open ( unf0,file=trim(fn0),status='OLD', &
                             form='unformatted',iostat=errop )
                  if (errop == 0) then
                     write(output_unit,2000) 'READING', trim(fn0)
                     do i=1,36
                        read (unf0)dumdum(i,1)
                     end do
                     do l=lmin,lmax
                        do m=1,l+1
                           read(unf0)ar_s(lmax-l+1,m)
                           read(unf0)ai_s(lmax-l+1,m)
                           read(unf0)br_s(lmax-l+1,m)
                           read(unf0)bi_s(lmax-l+1,m)
                        end do
                     end do
                     close(unf0)
                  else
                     write (output_unit, 3000) trim(fn0)
                !     call gem_error ( err,'read_markov_skeb', 'problem reading file' )
                  endif
            endif
            dim=ens_skeb_l*ens_skeb_m
            call RPN_COMM_bcast (dumdum,36,"MPI_INTEGER",0,"MULTIGRID", err)
            call RPN_COMM_bcast (ar_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (ai_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (br_s,dim,"MPI_REAL",0, "MULTIGRID", err)
            call RPN_COMM_bcast (bi_s,dim,"MPI_REAL",0, "MULTIGRID", err)
         else
            fstd=Ens_skeb_std
            tau = Ens_skeb_tau/2.146
            epsi=exp(-dt/tau)
! Bruit blanc en nombre d'ondes
            allocate (  fact1(lmin:lmax) )
            do l=lmin,lmax
               pspec2(l)=1.D0
            end do
!Normalisation du spectre pour que la variance du champ aléatoire soit std**2
            sumsp=0.D0
            do l=lmin,lmax
               sumsp=sumsp+pspec2(l)
            end do
               pspec2=pspec2/sumsp

               do l=lmin,lmax
! fact1(l)=fstd*SQRT(4.*pi/real((2*l+1))*pspec2(l))
                  w1= 4.d0*pi
                  w2= dble(2*l+1)*pspec2(l)
                  fact1(l)=fstd*SQRT(w1/w2)
            end do

            fact3 =(1.-epsi*epsi)/SQRT(1.+epsi*epsi)
! Random fonction generator
            dumdum(:,1)=0
            paiv  => dumdum(1,1)
            paiy  => dumdum(33,1)
            paiset=> dumdum(34,1)
            pagset=> dumdum(35,1)
            paidum=> dumdum(36,1)
            paidum=-Ens_mc_seed
! Valeurs initiales des coefficients spectraux
            ar_s(:,:)=0.d0
            br_s(:,:)=0.d0
            ai_s(:,:)=0.d0
            bi_s(:,:)=0.d0

            do l=lmin,lmax
               br_s(lmax-l+1,1)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)
               ar_s(lmax-l+1,1)=br_s(lmax-l+1,1)
               do m=2,l+1
                  br_s(lmax-l+1,m)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                  ar_s(lmax-l+1,m)=br_s(lmax-l+1,m)
                  bi_s(lmax-l+1,m)=fact1(l)*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                  ai_s(lmax-l+1,m)=bi_s(lmax-l+1,m)
               end do
            end do

            deallocate (fact1)

         end if

!cpdi --  specific heat
      cpdi=1./real(cpd_8)
! Deltax -- typical model gridlength
      deltax= 0.5d0 * (geomh_hx_8 + geomh_hy_8)*Dcst_rayt_8
      
      allocate(cc2(2 , nlat, lmax+1))
      allocate(wrk2( nlat * (nlon+2)))
      allocate(f2(nlon+2, nlat) )
      allocate( markov2(1:l_ni,1:l_nj))

      allocate( dsp_local(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate( dsp_dif(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate( dsp_gwd(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate( psi_local(l_minx:l_maxx,l_miny:l_maxy,l_nk))

      allocate(fg(7), fg2(7))
      allocate(fg1(l_nj,7))

      wrk2=0.d0 ; cc2=0.d0
      dsp_local=0.0;  dsp_dif =0.0; dsp_gwd=0.0; psi_local=0.0 

      sigma= dble(Ens_skeb_lam)*sqrt(-2.d0*log(dble(Ens_skeb_bfc)))/(2*pi_8);

      do j=1,l_nj
         fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hx_8*geomh_cy_8(j);fg(3)=-fg(5)
         fg(2)=2*fg(3);fg(1)=3*fg(3);fg(6)=2*fg(5);fg(7)=3*fg(5)
         fg=1/(sigma*sqrt(2*pi_8))*exp(-fg*fg/(2*sigma*sigma))
         fg=fg/sum(fg)
         fg1(j,1:7)=fg(1:7)
      end do

      fg2(4)=0.;fg2(5)=Dcst_rayt_8*geomh_hy_8;fg2(3)=-fg2(5)
      fg2(2)=2*fg2(3);fg2(1)=3*fg2(3);fg2(6)=2*fg2(5);fg2(7)=3*fg2(5)
      fg2=1/(sigma*sqrt(2*pi_8))*exp(-fg2*fg2/(2*sigma*sigma))
      fg2=fg2/sum(fg2)

 1000 format( &
           /,'INITIALIZE SCHEMES CONTROL PARAMETERS (S/R ENS_MARFIELD_SKEB)', &
           /,'======================================================')

 2000 format (/' MARKOV: ',a,' FILE ',a)

 3000 format(' S/R  ENS_MARFIELD_SKEB : problem in opening MRKV_PARAM_SKEB file)', &
              /,'======================================================')


      return

contains

      subroutine pleg(l, m, jlat, nlat, pls )
      use, intrinsic :: iso_fortran_env
      implicit none

      integer l,m ,i,j ,jlat ,nlat
      real(kind=REAL64)   pls
      real(kind=REAL64)  factor , x  ,lat, theta
      real(kind=REAL64) , dimension(0:l+1) :: pl
      real(kind=REAL64), parameter :: ZERO=0.0D0  , ONE_8=1.0d0 , TWO_8=2.0d0

!-------------------------------------------------------------------------
!
      if ( m < 0 .OR. m > l ) then
         print*, ' error :  m must non-negative and m <=l '
         stop
      end if

      lat=(-90.D0+90.D0/DBLE(nlat)+DBLE(jlat-1)*180.D0/DBLE(nlat))*pi_8/180.D0
      theta=pi_8/2.D0-lat
      x=DCOS(theta)

      pl=ZERO
      if ( m <= l ) then
         pl(m) = ONE_8
         factor = ONE_8
         do i = 1, m
            pl(m) = -pl(m)*factor*sqrt(1.d0 - x**2)/ &
                   dsqrt(dble((l+i)*(l-m+i)))
            factor = factor + 2.d0
         end do
         pls=pl(m)
      end if

      if ( m + 1 <= l ) then
         pls = x * dble ( 2 * m + 1 ) * pl(m)
         pl(m+1)=pls
      end if

      do j = m + 2, l
          pl(j) = ( x * dble (2*j-1) * pl(j-1) &
                     + dble (-j-m+1) * pl(j-2) ) &
                     / dble (j-m)
      end do

      pls=pl(l)

  end subroutine pleg


end subroutine ens_markov_skeb_init 
