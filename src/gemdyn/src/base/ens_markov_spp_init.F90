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

!**   s/r ens_markov_spp_ininitialize SPP/PTP markov chain
!
      subroutine ens_markov_spp_init()
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
#include <rmnlib_basics.hf>
!
       real,    external ::  gasdev
       real(kind=REAL64) :: polg

!
! nlat, nlon                 dimension of the Gaussian grid
! idum                       Semence du g?n?rateur de nombres al?atoires
!
      integer ::dlen, l ,m, n, nc,np, i, j, indx, ier, istat
      real    :: fstd, fstr, tau, sumsp , fact
      real    :: xfi(l_ni),yfi(l_nj)
      real(kind=REAL64)  :: rad2deg_8,  deg2rad_8, pri_8
      logical, save :: init_done=.false.
!
! paidum   pointer vers l'etat du generateur sauvegarde idum
      integer, pointer :: paiv,paiy,paiset,pagset,paidum
!
! dt   Pas de temps du mod?le (secondes)
! tau  Temps de d?corr?lation du champ al?atoire f(i,j) (secondes)
! eps  EXP(-dt/tau/2.146)
      real(kind=REAL64)   :: fact1, fmax, fmin , fmean
      real, dimension(2) :: spp_range
      integer, dimension(2) :: spp_trn
      integer ::  err, errop, ierr ,unf0, unf1, nch2d, spp_indx, stat,gmmstat 
      integer :: lmx,mmx,mch2d
      character(len=WB_MAXNAMELENGTH) :: prefix, key, spp_type
      character(len=1024) :: fn0, fn1
!
!-------------------------------------------------------------------
!
      nch2d = Ens_ptp_ncha + spp_ncha
      if (nch2d == 0) return
      dt=real(Cstv_dt_8)
      rad2deg_8=180.0d0/pi_8
      deg2rad_8=1d0/rad2deg_8
      itstep_s=step_dt*step_kount
      iperiod_iau = iau_period
      write_markov_l=(itstep_s==iperiod_iau)

      allocate(  vname(Ens_ptp_ncha+spp_ncha ))
      allocate ( vlmin(Ens_ptp_ncha+spp_ncha), vlmax(Ens_ptp_ncha+spp_ncha), &
                vnlon(Ens_ptp_ncha+spp_ncha), vnlat(Ens_ptp_ncha+spp_ncha), &
                vfmin(Ens_ptp_ncha+spp_ncha), vfmax(Ens_ptp_ncha+spp_ncha), &
                vfstd(Ens_ptp_ncha+spp_ncha), vfstr(Ens_ptp_ncha+spp_ncha), &
                vtau(Ens_ptp_ncha+spp_ncha) ,  veps(Ens_ptp_ncha+spp_ncha) )
      allocate(markov1(l_ni,l_nj,nch2d) )

      ! Generate tropically focused weights
      if (.not.associated(tropwt)) then
         allocate(tropwt(l_ni,l_nj))
         if (Ens_spp_rhsint_lat > 0.) then
            where (abs(geomh_latrx(1:l_ni,1:l_nj)) > Ens_spp_rhsint_lat)
               tropwt(:,:) = 0.
            elsewhere
               tropwt(:,:) = cos( (90./Ens_spp_rhsint_lat) * &
                    deg2rad_8*geomh_latrx(1:l_ni,1:l_nj) )**2
            end where
         else
            tropwt(:,:) = 1.
         endif
      endif

      if (.not.init_done) then
         if (Lun_out > 0) then
            write( Lun_out,1000)
         end if
         init_done=.true.
      end if
      spp_indx = 0
      do nc=1,nch2d
         if (nc <= Ens_ptp_ncha) then
            ! Add chains for PTP
            vname(nc) = 'PTP'
            vlmin(nc) = Ens_ptp_trnl(nc)
            vlmax(nc) = Ens_ptp_trnh(nc)
            vnlon(nc) = Ens_ptp_nlon(nc)
            vnlat(nc) = Ens_ptp_nlat(nc)
            vfmin(nc) = Ens_ptp_min(nc)
            vfmax(nc) = Ens_ptp_max(nc)
            vfstd(nc) = Ens_ptp_std(nc)
            vfstr(nc) = Ens_ptp_str(nc)
            vtau(nc) = Ens_ptp_tau(nc)
         else
            ! Add chains for SPP
            spp_indx = spp_indx + 1
            vname(nc) = spp_list(spp_indx)
            prefix = 'spp/'//trim(vname(nc))//'/'
            stat = clib_toupper(vname(nc))
            stat = WB_OK
            key = trim(prefix)//'spp_trn'
            stat = min(wb_get(key, spp_trn, dlen), stat)
            vlmin(nc) = minval(spp_trn)
            vlmax(nc) = maxval(spp_trn)
            key = trim(prefix)//'spp_nlon'
            stat = min(wb_get(key, vnlon(nc)), stat)
            key = trim(prefix)//'spp_nlat'
            stat = min(wb_get(key, vnlat(nc)), stat)
            key = trim(prefix)//'spp_std'
            stat = min(wb_get(key, vfstd(nc)), stat)
            key = trim(prefix)//'spp_str'
            stat = min(wb_get(key, vfstr(nc)), stat)
            key = trim(prefix)//'spp_tau'
            stat = min(wb_get(key, vtau(nc)), stat)
            key = trim(prefix)//'spp_type'
            stat = min(wb_get(key, spp_type), stat)
            if (spp_type == 'DISCRETE') then
               vfmin(nc) = 0. ; vfmax(nc) = 1.
            else
               key = trim(prefix)//'spp_range'
               stat = min(wb_get(key, spp_range, dlen), stat)
               vfmin(nc) = minval(spp_range)
               vfmax(nc) = maxval(spp_range)
            endif
            if ( WB_IS_ERROR(stat)) then
               write(Lun_out, *) 'Error retrieving chain specifications for '// &
                    trim(spp_list(spp_indx))
               return
            endif
         endif
         vtau(nc) = vtau(nc) / 2.146
         veps(nc) = exp(-dt/vtau(nc))
      enddo

      lmx = max(Ens_ptp_lmax, spp_lmax)
      mmx = max(Ens_ptp_mmax, spp_mmax)
      mch2d = 2*(MAX2DC+MAX_NSPP)

         do nc=1,nch2d
            lmin = vlmin(nc)
            lmax = vlmax(nc)
            nlat = vnlat(nc)
            nlon = vnlon(nc)

            ! Associated Legendre polynomials
            plp(:,:,:,nc)=0.D0
            do l=lmin,lmax
               fact=DSQRT((2.D0*DBLE(l)+1.D0)/(4.D0*pi_8))
               do m=0,l
                  do j=1,nlat/2
                     call pleg (l, m, j, nlat, polg)
                     plp(j,lmax-l+1,m+1,nc)=polg*fact
                     plp(nlat-j+1,lmax-l+1,m+1,nc)=polg*fact*(-1.D0)**(l+m)
                  end do
               end do
            end do

            ! Memory allocation
            allocate(cc1(2 , nlat, lmax+1 ,nch2d))
            allocate(wrk1( nlat * (nlon+2),nch2d))
            allocate(f1(    nlon, nlat    ,nch2d))
            allocate(f1_str(nlon, nlat    ,nch2d))
            cc1=0.0 ; wrk1=0.0 
         enddo

         !Initialise spectral coeffs and stochastic params
         if(Ens_recycle_mc) then
            !Read saved stochastic numbers and spectral coeffs ar,br.ai,bi
            err=0
            np=0
            if (ptopo_myproc==0 ) then
               do nc=1,nch2d
                 unf0=0
                 np=nc+1
                 fn0= trim(Path_input_S)//'/MODEL_INPUT'//'/MRKV_SPP_'//trim(vname(nc))//'.bin'
                  open ( unf0,file=trim(fn0),status='OLD', &
                             form='unformatted',iostat=errop )
                  if (errop == 0) then
                     err=0
                     write(output_unit,2000) 'READING', trim(fn0)
                     do i=1,36
                        read (unf0)dumdum(i,np)
                     end do
                     do l=lmin,lmax
                        do m=1,l+1
                           read(unf0)ar_p(lmax-l+1,m,nc)
                           read(unf0)ai_p(lmax-l+1,m,nc)
                           read(unf0)br_p(lmax-l+1,m,nc)
                           read(unf0)bi_p(lmax-l+1,m,nc)
                        end do
                     end do
                     close(unf0)
                  else
                     err = -1
                     write(output_unit, *) 'Error: problem opening Markov file '//trim(fn0)
                  endif
               enddo
            endif
        !    call gem_error(err,'Ens_marfield_ptp','Error in reading Markov files')


            call RPN_COMM_bcast (dumdum,36*mch2d,"MPI_INTEGER",0,"MULTIGRID", ierr)
            call RPN_COMM_bcast (ar_p, lmx*mmx*nch2d, "MPI_REAL",0, "MULTIGRID", ierr)
            call RPN_COMM_bcast (ai_p, lmx*mmx*nch2d, "MPI_REAL",0, "MULTIGRID", ierr)
            call RPN_COMM_bcast (br_p, lmx*mmx*nch2d, "MPI_REAL",0, "MULTIGRID", ierr)
            call RPN_COMM_bcast (bi_p, lmx*mmx*nch2d, "MPI_REAL",0, "MULTIGRID", ierr)

         else
            np=0
            do nc=1,nch2d
               lmin = vlmin(nc)
               lmax = vlmax(nc)
               nlon = vnlon(nc)
               nlat = vnlat(nc)
               fstd = vfstd(nc)
               tau  = vtau(nc)
               eps  = veps(nc)

               !  White noise in wave number

               if(.not. allocated(pspec1)) allocate ( pspec1(lmin:lmax,nch2d))
               if(.not. allocated(fact2) ) allocate ( fact2(nch2d) )

               do l= lmin,lmax
                  pspec1(l,nc)=1.D0
               end do

	       !Normalization of the spectrum so that the variance of the random field is
               sumsp=0.D0
               do l=lmin,lmax
                  sumsp=sumsp+pspec1(l,nc)
               end do
               do l=lmin,lmax
               	   pspec1(l,nc)=pspec1(l,nc)/sumsp
               end do

               fact2(nc) =(1.-eps*eps)/SQRT(1.+eps*eps)

               ! Random function generator
               np=nc+1
               dumdum(:,np)=0
               paiv  => dumdum(1,np)
               paiy  => dumdum(33,np)
               paiset=> dumdum(34,np)
               pagset=> dumdum(35,np)
               paidum=> dumdum(36,np)
               paidum=-(Ens_mc_seed + 1000*nc)

               ! Initial values  of spectral coefficients
               ar_p(:,:,nc)=0.d0;br_p(:,:,nc)=0.d0;ai_p(:,:,nc)=0.d0;bi_p(:,:,nc)=0.d0

               do l=lmin,lmax
                  fact1=fstd*SQRT(4.*pi_8/real((2*l+1))*pspec1(l,nc))
                  br_p(lmax-l+1,1,nc)=fact1*gasdev(paiv,paiy,paiset,pagset,paidum)
                  ar_p(lmax-l+1,1,nc)=br_p(lmax-l+1,1,nc)
                  do m=2,l+1
                     br_p(lmax-l+1,m,nc)=fact1*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                     ar_p(lmax-l+1,m,nc)=br_p(lmax-l+1,m,nc)
                     bi_p(lmax-l+1,m,nc)=fact1*gasdev(paiv,paiy,paiset,pagset,paidum)/sqrt(2.)
                     ai_p(lmax-l+1,m,nc)=bi_p(lmax-l+1,m,nc)
                  end do
               end do
            end do
      end if


1000 format( &
           /,'INITIALIZE SCHEMES CONTROL PARAMETERS (S/R ENS_MARFIELD_PTP_SPP)', &
           /,'======================================================')
1005 format (/' Problem:Cannot find Markov file: ',a)
2000 format (/' MARKOV: ',a,' FILE ',a)

3000 format (/' S/R  ENS_MARFIELD_PTP : problem in opening MRKV_PARAM_PTP_SPP file: ',a)

!3000 format(' S/R  ENS_MARFIELD_PTP : problem in opening MRKV_PARAM_PTP_SPP file)', &
!              ',a,')

4000 format(' S/R  ENS_MARFIELD_PTP : problem in WRITING  MRKV_PARAM_PTP_SPP file)', &
              /,'======================================================')

return

contains

 subroutine pleg(l, m, jlat, nlat , plg )
      use, intrinsic :: iso_fortran_env
 implicit none

      integer l,m ,i,j , jlat, nlat
      real(kind=REAL64)  factor , x , plg , lat, theta
      real(kind=REAL64) , dimension(0:l+1) :: pl
      real(kind=REAL64), parameter :: ZERO=0.0D0  , ONE_8=1.0d0 , TWO_8=2.0d0
!
      if ( m < 0 .OR. m > l ) then
         print*, ' error :  m must non-negative and m <l '
         stop
      end if

      lat=(-90.D0+90.D0/DBLE(nlat)+DBLE(jlat-1)*180.D0/DBLE(nlat))*pi_8/180.D0
      theta=pi/2.D0-lat
      x=DCOS(theta)

      pl=ZERO
      if ( m <= l ) then
         pl(m) = ONE_8
         factor = ONE_8

         do i = 1, m
            pl(m) = -pl(m)*factor*sqrt(ONE_8 - x**2)/ &
                   dsqrt(dble((l+i)*(l-m+i)))
            factor = factor + TWO_8
         end do
         plg=pl(m)
      end if

      if ( m + 1 <= l ) then
         plg = x * dble ( 2 * m + 1 ) * pl(m)
         pl(m+1)=plg
      end if

      do j = m + 2, l
          pl(j) = ( x * dble (2*j-1) * pl(j-1) &
                     + dble (-j-m+1) * pl(j-2) ) &
                     / dble (j-m)
      end do

      plg=pl(l)

  end subroutine pleg

end subroutine ens_markov_spp_init 
