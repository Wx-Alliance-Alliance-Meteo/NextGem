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

!**s/r proc_split_step - Advance one timestep using process splitting method

      subroutine proc_split_step ()
      use gem_options
      use step_options
      use theo_options
      use dcmip_options
      use dynkernel_options
      use mem_nest
      use omp_timing
      use gmm_vt1
      use gmm_pw
      implicit none

      logical, save :: done=.false., psadj_L=.false.
!!$omp threadprivate(psadj_L,done)
!
!     ---------------------------------------------------------------
!
      call gtmg_start (10, 'DYNSTEP', 1)
      
      if (.not. psadj_L) call psadj_p0init ( Step_kount, psadj_L )

!$OMP BARRIER

      if (.not. done) then
!!$omp single
         call set_sol ()
!!$omp end single
      endif
      done =.true.
      
!!$omp single
      Nest_current_dayfrac_8= dble(Step_kount) * Cstv_dt_8
      Vtopo_L= (Vtopo_L.and.(Nest_current_dayfrac_8<=Vtopo_start_8+Vtopo_ndt_8))
!!$omp end single

      if (Dynamics_sw_L) then
         call SW_dynstep ()
      else
         call dynstep ()
      endif
      
      call gtmg_stop (10)

      if (Ctrl_testcases_L) call canonical_cases ("VRD")

      call hzd_main ()

      if (Grd_yinyang_L) call yyg_blend () !ut1,vt1,zdt1

      call pw_update_GW ()
      call pw_update_UV ()
      call pw_update_T  ()

      call out_dyn (.false., .true.) ! output for cascade
      
      call itf_phy_step (Step_kount, Lctl_step)

      call canonical_cases ("PHY,ERR")
!!$omp single
      call iau_apply (Step_kount)
!!$omp end single

      if ( Ctrl_phyms_L.or.Dcmip_physics_L ) then
         call tt2tvirt (tt1(l_minx,l_miny,1), pw_tt_plus, 1,l_ni, 1,l_nj)
         call itf_phy_UVupdate()
      end if
      call pw_update_GW ()
      
      if (Grd_yinyang_L) call yyg_xchng_all()
!
!     ---------------------------------------------------------------
!
      return
      end subroutine proc_split_step
