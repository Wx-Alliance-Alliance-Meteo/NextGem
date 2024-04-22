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

!**s/r gem_run - Compute all timesteps in sequence

      subroutine gem_run (F_rstrt_L)
      use cstv
      use dynkernel_options
      use lun
      use rstr
      use svro_mod
      use gem_options
      use step_options
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(OUT) :: F_rstrt_L

!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_rstrt_L     O         Is a restart required
!----------------------------------------------------------------

      logical, external :: gem_muststop
      integer, external :: model_timeout_alarm
      character(len=16) :: datev
      integer :: stepf,last_step, seconds_since, err
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: sec_in_day = 86400.0d0
!
!     ---------------------------------------------------------------
!
      dayfrac = dble(Step_kount) * Cstv_dt_8 / sec_in_day
      call incdatsd (datev,Step_runstrt_S,dayfrac)

      if (Lun_out > 0) write (Lun_out,900) datev

      call blocstat (.true.)

      call gemtime ( Lun_out, 'STARTING TIME LOOP', .false. )

      stepf= Step_total
      last_step = Step_total + Step_initial

      F_rstrt_L = .false.
      if ( .not. Rstri_rstn_L ) then
         call out_outdir()
         call out_dyn (.true., .true.)
         call OUTs_end (.false.)
         if (gem_muststop (stepf)) then
            seconds_since = model_timeout_alarm(Step_alarm)
            if (Lun_out > 0) write(Lun_out,4000) Lctl_step
            return
         end if
      end if

      call canonical_cases ("ERR")

      do while (Step_kount < stepf)

         seconds_since = model_timeout_alarm(Step_alarm)

         Lctl_step= Lctl_step + 1  ;  Step_kount= Step_kount + 1
         call MPI_barrier (COMM_GRID, err)

         if (Lun_out > 0) write (Lun_out,1002) Lctl_step, last_step

         call out_outdir()

!!$omp parallel
         call proc_split_step ()
!!$omp end parallel

         call out_dyn (.true., .false.) ! regular output

         call blocstat (.false.)

         if (Lun_out > 0) write(Lun_out,3000) Lctl_step

         call save_restart()

         F_rstrt_L= gem_muststop (stepf)

         if (F_rstrt_L) exit

      end do

      seconds_since = model_timeout_alarm(Step_alarm)

      if (Lun_out > 0) write(Lun_out,4000) Lctl_step

 900  format (/'STARTING THE INTEGRATION WITH THE FOLLOWING DATA: VALID ',a)
 1002 format(/,'FISL-H: PERFORMING TIMESTEP #',I9,' OUT OF ',I9, &
             /,'=========================================================')
 1003 format(/,'PICSL: PERFORMING TIMESTEP #',I9,' OUT OF ',I9, &
             /,'=========================================================')
 3000 format(/,'THE TIME STEP ',I8,' IS COMPLETED')
 4000 format(/,'GEM_RUN: END OF THE TIME LOOP AT TIMESTEP',I8, &
             /,'===================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
