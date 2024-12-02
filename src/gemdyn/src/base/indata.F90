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

!**s/r indata - Read and process the input data at
!               beginning of integration

      subroutine indata()
      use gem_options
      use gmm_contiguous
      use dynkernel_options
      use adz_mem
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use inp_options
      use lun
      use metric
      use geomh
      use ver
      use step_options
      use theo_options
      use tr3d
      use svro_mod
      use outp
      use omp_timing
      use mem_tstp
      use mem_tracers
      use cstv
      implicit none

      logical :: synthetic_data_L
      integer :: i,j,k,dimens,km1,km2,km3,kp1,kp2,kp3
      integer :: HLT_start, HLT_end, local_np
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: p0
      real retval,v1,v2,v3,v4,v5,v6
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      synthetic_data_L = (Ctrl_theoc_L .or. Ctrl_testcases_L)

      if (synthetic_data_L) then
         call synthetic_data ()
      else
         call gtmg_start ( 71, 'INITIAL_input', 2)

         call get_topo ()

         Inp_src_GZ_L= .false. ; Inp_gtmg= (/82,71/)
         call inp_data (pw_uu_plus,pw_vv_plus,wt1(l_minx,l_miny,1),&
             pw_tt_plus,qt1(l_minx,l_miny,1),zdt1(l_minx,l_miny,1),&
             p0,trt1,fis0,orols,.false.,Step_runstrt_S,&
             l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr )
         do j=1-G_haloy+1,l_nj+G_haloy-1
            do i=1-G_halox+1,l_ni+G_halox-1
               fis0u (i,j) = (fis0(i+1,j)+fis0(i,j))*0.5
               fis0v (i,j) = (fis0(i,j+1)+fis0(i,j))*0.5
               orolsu(i,j) = (orols(i+1,j)+orols(i,j))*0.5
               orolsv(i,j) = (orols(i,j+1)+orols(i,j))*0.5
            end do
         end do
      
         dimens=(l_maxx-l_minx+1)*(l_maxy-l_miny+1)*G_nk
         call bitflip ( pw_uu_plus, pw_vv_plus, pw_tt_plus, &
                        perturb_nbits, perturb_npts, dimens )
         call gtmg_stop  ( 71 )
      end if

      call gemtime ( Lun_out, 'AFTER INITIAL INPUT', .false. )

!!$omp parallel private (local_np, HLT_start, HLT_end)
      if (Grd_yinyang_L) then
         if (.not.synthetic_data_L) call yyg_xchng_hlt (p0, l_minx,l_maxx,l_miny,l_maxy, &
                                                 l_ni,l_nj, 1, .false., 'CUBIC', .true.)
         call yyg_xchng_all()
      else
         call HLT_split (1, 1, local_np, HLT_start, HLT_end)
      end if
      
      if(.not.Dynamics_sw_L) call vertical_metric ()

!!$omp do
      do j=1-G_haloy+1,l_nj+G_halox-1
         do i=1-G_halox+1,l_ni+G_halox-1
            me_full(i,j) = fis0(i,j) 
            me_large(i,j)= orols(i,j) 
         end do
      end do
!!$omp end do nowait
!!$omp do
      do i= 1, ubound(trt0,1)
         trt0(i) = trt1(i)
         trt2(i) = trt1(i)
      end do
!!$omp end do nowait

      if (.not. synthetic_data_L) then
!!$omp single
         call tt2tvirt (tt1(l_minx,l_miny,1), pw_tt_plus, -G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy)
!!$omp end sin
         call HLT_split (1, G_nk, local_np, HLT_start, HLT_end)
         call gem_xch_halo (pw_uu_plus(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
         call gem_xch_halo (pw_vv_plus(l_minx,l_miny,HLT_start),l_minx,l_maxx,l_miny,l_maxy,local_np,-1)

         call hwnd_stag2 ( ut1(l_minx,l_miny,1),vt1(l_minx,l_miny,1), pw_uu_plus,pw_vv_plus   ,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk   ,&
                         1-G_halox*west ,l_niu+G_halox*east ,&
                         1-G_haloy*south,l_njv+G_haloy*north, .true. )
         call derivate_data (zdt1(l_minx,l_miny,1),wt1(l_minx,l_miny,1), &
                              ut1(l_minx,l_miny,1),vt1(l_minx,l_miny,1), &
                           tt1(l_minx,l_miny,1),p0,qt1(l_minx,l_miny,1), &
                         fis0,orols,GVM%zmom_8,GVM%ztht_8,l_minx,l_maxx, &
          l_miny,l_maxy, G_nk,.not.Inp_zd_L, .not.Inp_w_L, .not.Inp_qt_L )
      endif
      call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,&
                      pw_log_pt, pw_pm_plus_8,pw_p0_plus_8      ,&
                      l_minx,l_maxx,l_miny,l_maxy,l_nk,1 )
      call pw_update_GW ()
      call psadj_init()

!!$omp do collapse(2)
      do k=1, G_nk
         do j=1-G_haloy,l_nj+G_halox
            do i=1-G_halox,l_ni+G_halox
               pw_uu_moins(i,j,k) = pw_uu_plus(i,j,k) 
               pw_vv_moins(i,j,k) = pw_vv_plus(i,j,k) 
               pw_tt_moins(i,j,k) = pw_tt_plus(i,j,k) 
               pw_gz_moins(i,j,k) = pw_gz_plus(i,j,k) 
               pw_pm_moins(i,j,k) = pw_pm_plus(i,j,k) 
               pw_pt_moins(i,j,k) = pw_pt_plus(i,j,k) 
            end do
         end do
      end do
!!$omp end do nowait
!!$omp do
      do j=1-G_haloy,l_nj+G_halox
         do i=1-G_halox,l_ni+G_halox
            pw_pm_moins(i,j,G_nk+1) = pw_pm_plus(i,j,G_nk+1) 
            pw_pt_moins(i,j,G_nk+1) = pw_pt_plus(i,j,G_nk+1)
            pw_me_moins(i,j) = pw_me_plus(i,j) 
            pw_p0_moins(i,j) = pw_p0_plus(i,j) 
            udiag(i,j) = pw_uu_plus(i,j,G_nk)
            vdiag(i,j) = pw_vv_plus(i,j,G_nk)
            tdiag(i,j) = pw_tt_plus(i,j,G_nk)
            qdiag(i,j) = tracers_P(Tr3d_hu)%pntr(i,j,G_nk)
         end do
      end do
!!$omp end do
!!$omp end parallel

      call out_outdir()

      call iau_apply (0)

      if (.not. Grd_yinyang_L) call nest_init()
      
!!$omp parallel
      call itf_phy_step (0,Lctl_step)
!!$omp end parallel
      call OUTs_end (.false.)
      
      dynt0= dynt1
      dynt2= dynt1
      allocate (var_init(ubound(dynt1,1)))
      var_init= dynt1

!     call nest_glbstat ((/'now','deb','fin'/),3)
!
!     ---------------------------------------------------------------
!
 1000 format(/,'TREATING INITIAL CONDITIONS  (S/R INDATA)',/,41('='))

      return
      end subroutine indata

      subroutine bitflip (u,v,t,nbits,npts,n)
      implicit none

      integer, intent(in) :: n,nbits,npts
      real, dimension(n), intent(inout) :: u, v, t

      integer stride,i
      integer, dimension(n) :: u_bits, v_bits, t_bits
!
! ---------------------------------------------------------------------
!
      u_bits  = transfer(u, u_bits)
      v_bits  = transfer(v, v_bits)
      t_bits  = transfer(t, t_bits)

      if (nbits < 1) return

      stride = min(max(1,npts),n)

      do i=1,n,stride
         u_bits(i) = xor(u_bits(i), nbits)
         v_bits(i) = xor(v_bits(i), nbits)
         t_bits(i) = xor(t_bits(i), nbits)
      end do

      u = transfer(u_bits, u)
      v = transfer(v_bits, v)
      t = transfer(t_bits, t)
!
! ---------------------------------------------------------------------
!
      return
      end subroutine bitflip
