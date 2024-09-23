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

!**s/r dynstep - Advance the dynamics by one timestep

      subroutine dynstep ()
      use ISO_C_BINDING
      use step_options
      use theo_options
      use gmm_pw
      use gmm_contiguous
      use mem_tstp
      use adz_mem
      use mem_tracers
      use tr3d
      implicit none

      logical :: flag=.false.
      integer :: n
!
!     ---------------------------------------------------------------
!
      ds_i0= 1   +pil_w
      ds_in= l_ni-pil_e
      ds_j0= 1   +pil_s
      ds_jn= l_nj-pil_n
      ds_k0= 1
      ds_kn= l_nk
      
      call pw_switch    ()

      if ((Step_kount.eq.1).and.Euler_step_one) then

         call tstpdyn(0.75d0*Cstv_dt_8) ! EULER

         call t02t1 ()
         call pw_update_UV ()
         call pw_switch ()

         call t02t2 () ! tracers have been modified by physics at t=0

         call tstpdyn(0.5d0*Cstv_dt_8 ) ! BDF2

         !for first step tracer advection only
         call adz_traject (Cstv_dt_8,1,.false.,.false.)

      else
            
         call t02t2 ()
         if ((Step_kount.eq.2).and.Euler_step_one)  dynt2= var_init
         call tstpdyn(Cstv_dt_8) ! BDF2

      endif

      if (.not.Grd_yinyang_L) then
         if (Ctrl_theoc_L) then
            call height_spongeH() ! mountain cases
            call mirror () ! bubble case
         endif
      endif

      call adz_tracers_interp ()

      call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,&
                      pw_log_pt,pw_pm_plus_8,pw_p0_plus_8       ,&
                      l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )
      
      if ( Schm_psadj > 0 ) call psadj ( Step_kount )

      call adz_tracers_massfixing ()
      
      if (Grd_yinyang_L) then
         do n=1, Tr3d_ntr
            call yyg_xchng_hlt (tracers_M(n)%pntr, l_minx,l_maxx,&
               l_miny,l_maxy,l_ni,l_nj,G_nk, .true., 'CUBIC', .false.)
         end do
      else
         call nest_HOR_gwa ()
         call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,&
                         pw_log_pt, pw_pm_plus_8,pw_p0_plus_8      ,&
                         l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )
      endif

      call t02t1 ()
!     ---------------------------------------------------------------
!
      return
      end subroutine dynstep
