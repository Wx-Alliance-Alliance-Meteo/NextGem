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

!**s/r ens_skeb_apply -Apply skeb to UT1,VT1,TT1  
!
      subroutine ens_skeb_apply (F_ugwdt1,F_vgwdt1,F_difut1,F_difvt1)
      use dcst
      use ens_param
      use ens_gmm_var
      use ens_options
      use gem_options
      use HORgrid_options
      use geomh
      use gmm_vt1
      use tdpack
      use glb_ld
      use rmn_gmm
      use var_gmm
      use omp_lib
      implicit none
!
      real    F_ugwdt1(l_minx:l_maxx,l_miny:l_maxy,l_nk), F_vgwdt1(l_minx:l_maxx,l_miny:l_maxy,l_nk)
      real    F_difut1(l_minx:l_maxx,l_miny:l_maxy,l_nk), F_difvt1(l_minx:l_maxx,l_miny:l_maxy,l_nk)

!author
!     Lubos Spacek - rpn - apr 2005

      integer i, j, k, i0,in,j0,jn 
      integer :: local_np, HLT_start, HLT_end
      real    dummy
      real  , dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: dummy_tab
      real(kind=REAL64), parameter :: half= 0.5

!     ---------------------------------------------------------------

      call HLT_split (1, l_nk, local_np, HLT_start, HLT_end)

      call gem_xch_halo ( ut1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      call gem_xch_halo ( vt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      call gem_xch_halo ( F_difut1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      call gem_xch_halo ( F_difvt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      call gem_xch_halo ( F_ugwdt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      call gem_xch_halo ( F_vgwdt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )

!  Calculate kinetic energy of diffusion tendency

!  Diffusion backscatter

      if(Ens_stat)then
!!$omp single
         call glbstat(ut1,'URT1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(vt1,'VRT1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(F_difut1,'DUT1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(F_difvt1,'DVT1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(F_ugwdt1,'UGW1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(F_vgwdt1,'VGW1','AT BEG',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if
           
      if(Ens_skeb_dif)then
!!$omp do collapse(2) 
         do k= 1, l_nk
            do j= l_miny, l_maxy
               do i= l_minx,l_maxx
            	  F_difut1(i,j,k) = ut1(i,j,k)*F_difut1(i,j,k)
            	  F_difvt1(i,j,k) = vt1(i,j,k)*F_difvt1(i,j,k)
               end do
            end do
         end do
!!$omp enddo 

!!$omp do collapse(2) 
         do k= 1, l_nk
            do j= 0, l_nj
               do i= 0, l_ni
            	  dsp_dif(i,j,k) = half*sqrt(( F_difut1(i,j,  k) + F_difut1(i,j-1,  k))**2 &
                                          +  ( F_difvt1(i,j-1,k) + F_difvt1(i-1,j-1,k))**2 )
               end do
            end do
         end do
!!$omp enddo  nowait
      end if

      if(Ens_stat)then
!!$omp single
         call glbstat(F_difut1,'DUT1','DIFF',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat(F_difvt1,'DVT1','DIFF',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

!  Gravity wave drag backscatter

      if(Ens_skeb_gwd)then
!!$omp do 
         do k= 1, l_nk
	    do j= l_miny, l_maxy
               do i= l_minx,l_maxx
                  F_difut1(i,j,k)=ut1(i,j,k)*F_ugwdt1(i,j,k)
                  F_difvt1(i,j,k)=vt1(i,j,k)*F_vgwdt1(i,j,k)
               end do
            end do
         end do
!!$omp enddo

!!$omp do 
         do k= 1, l_nk
	    do j= 0, l_nj
               do i= 0, l_ni
           	  dsp_gwd(i,j,k) = half*abs( (F_difut1(i,j  ,k)+F_difut1(i,j-1  ,k)) &
                                         + (  F_difvt1(i,j-1,k)+F_difvt1(i-1,j-1,k)) )
               end do
            end do
         end do
!!$omp enddo
      end if

      if(Ens_stat)then
!!$omp single
         call glbstat (F_difut1,'DUT1','GWD',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (F_difvt1,'DVT1','GWD',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

!!$omp do collapse(2) 
         do k= 1, l_nk
	    do j= 0, l_nj
               do i= 0, l_ni
                  dsp_local(i,j,k)=dsp_dif(i,j,k)+dsp_gwd(i,j,k)
               end do
            end do
         end do
!!$omp enddo

      if(Ens_stat)then
!!$omp single
         call glbstat (dsp_dif,'DSP','DIF',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (dsp_gwd,'DSP','GWD',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (dsp_local,'DSP','TOT',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

!     Apply 2D Gaussian filter

      call ens_filter_gauss(dsp_local)

      if(Ens_stat)then
!!$omp single
         call glbstat (dsp_local,'DSP','FLTTOT',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

      if(Ens_skeb_alpt/=0.0)then
!!$omp do 
         do k=1,l_nk
            do j=1,l_nj
               do i=1,l_ni
                  tt1(i,j,k)=tt1(i,j,k)+Ens_skeb_alpt*cpdi*dsp_local(i,j,k)*markov2(i,j)
               enddo
            enddo
         enddo
!!$omp enddo nowait
      end if

!!$omp do 
         do k=1,l_nk
            do j=1,l_nj
               do i=1,l_ni
                  dsp_local(i,j,k)=sqrt(dsp_local(i,j,k))* markov2(i,j)
               enddo
            enddo
         enddo
!!$omp enddo

      if(Ens_stat)then
!!$omp single
         call glbstat (dsp_local,'DSP','FLTTO2',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

      call gem_xch_halo ( dsp_local(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )


!  Compute Curl of filtered field

!!$omp do 
      do k=1,l_nk
         do j= 1, l_njv
            do i= i0, l_ni
               F_difut1(i,j,k) = (dsp_local(i,j,k)-dsp_local(i-1,j,k))*geomh_invDXv_8(j)
               F_difut1(i,j,k) = Ens_skeb_alph*deltax*F_difut1(i,j,k)
               vt1(i,j,k)    = vt1(i,j,k)+F_difut1(i,j,k)
            end do
         end do
      end do
!!$omp enddo nowait

!!$omp do 
      do k=1,l_nk
         do j= 1, l_nj !-pil_n
            do i= 1, l_niu
               F_difvt1(i,j,k) = (dsp_local(i,j,k) - dsp_local(i,j-1,k)) * geomh_invDY_8
               F_difvt1(i,j,k)= -Ens_skeb_alph*deltax*F_difvt1(i,j,k)
               ut1(i,j,k) = ut1(i,j,k)+F_difvt1(i,j,k)
            end do
         end do
      end do
!!$omp enddo 

      if(Ens_stat)then
!!$omp single
         call glbstat (ut1,'URT1','AT END',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (vt1,'VRT1','AT END',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (F_difut1,'DUT1','AT END',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         call glbstat (F_difvt1,'DVT1','AT END',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
!!$omp end single
      end if

      if(Ens_skeb_div)then
         call gem_xch_halo ( F_difut1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
         call gem_xch_halo ( F_difvt1(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
!!$omp single

         call cal_div ( ensdiv, F_difut1, F_difvt1, 0, dummy,&
                       l_minx,l_maxx,l_miny,l_maxy, l_nk )
         call cal_vor ( ensvor, dummy_tab, F_difut1, F_difvt1, 0, dummy, .false., &
                      l_minx,l_maxx,l_miny,l_maxy, l_nk )

       ! call gem_error(-1,'ens_filter','Must check calls to cal_div and cal_vor')

         if(Ens_stat)then
            call glbstat (ensdiv,'DVRG','FRCING',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
            call glbstat (ensvor,'VORT','FRCING',l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
         end if
!!$omp end single
      end if

      end subroutine ens_skeb_apply
