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

!**s/r stop_world_view - Update status file and stop MPI

      subroutine stop_world_view
      use phy_itf
      use step_options
      use gem_options
      use sol_options
      use ctrl
      use lun
      use path
      use clib_itf_mod
      use adz_mem
      use inp_mod
      use svri_mod
      use svro_mod
      use ptopo
      use version
      use omp_timing
      use omp_timing
      implicit none

#include <rmnlib_basics.hf>
      include 'mpif.h'
      character(len=256) :: postjob_S
      character(len=16 ) :: dateV_S
      logical continue_L
      integer i,err,send(5),tag,ireq
      integer :: status(mpi_status_size,Ptopo_numproc)
!
!-------------------------------------------------------------------
!
      if (Ctrl_phyms_L) err = phy_terminate()

      continue_L= (Step_kount < Step_total)

      if (Lun_out > 0) then
         write (postjob_S,34) Lctl_step
         postjob_S='_endstep='//trim(postjob_S)
         call write_status_file3 (postjob_S)
         if (continue_L) then
            call write_status_file3 ('_status=RS')
         else
            call write_status_file3 ('_status=ED')
         end if
         call close_status_file3 ()
      end if

      if (Lun_out > 0) call out_stat ()

      if (.not. continue_L) then
         err = clib_remove('gem_restart')
         err = clib_remove('gmm_restart')
      end if

      call gtmg_stop ( 1 )
      call gtmg_terminate( Ptopo_myproc, 'GEMDM' )
!     Free all one-sided MPI windows
      call MPI_Win_free  (Inp_window,err)

      if (OUTs_server_L) then
         call gemtime ( Lun_out, 'ENDING ServerO', .false. )
         ireq= 3
         call MPI_waitall (ireq,OUTs_request,status,err)
         if (OUTs_comm_L) then
            send= -1 ; tag= 601
            call MPI_barrier (OUTs_COMM, err)
            call MPI_send ( send,size(send),MPI_INTEGER,&
                            OUTs_1on1,tag,OUTs_COMM,err )
         endif
         call MPI_barrier (COMM_multigrid,err)
         call gemtime ( Lun_out, 'ENDING ServerO DONE', .false. )
      endif

      if (INs_server_L) then
         call gemtime ( Lun_out, 'ENDING ServerI', .false. )
         if (INs_comm_L) then
            tag= 801
            ireq= -1
            call MPI_send (ireq,1,MPI_INTEGER,&
                    INs_1on1,tag,INs_COMM,err )
         endif
      endif
      
      if (ADZ_OD_L) then
          call MPI_Win_free  (Adz_Win_list,err)
          call MPI_Win_free  (Adz_Win_pos,err)
          call MPI_Win_free  (Adz_expq%winreqs,err)
          call MPI_Win_free  (Adz_expu%winreqs,err)
          call MPI_Win_free  (Adz_expv%winreqs,err)
          call MPI_Win_free  (Adz_expt%winreqs,err)
          call MPI_Win_free  (Adz_offs_win,err)
          call MPI_Win_free  (Adz_expq%wintraj,err)
          call MPI_Win_free  (Adz_expu%wintraj,err)
          call MPI_Win_free  (Adz_expv%wintraj,err)
          call MPI_Win_free  (Adz_expt%wintraj,err)
          call MPI_Win_free  (Adz_wincor,err)
      endif

      if (Lun_out > 0) then
         if (ADZ_OD_L) then
            do i=1, size(nexports)
               print*, 'EXPORTS: ',i,nexports(i)/Ptopo_numproc
            enddo
         endif
      end if

      call memusage ( Lun_out )
      call gemtime ( Lun_out, 'ENDING MY GEMDM DOMAIN', .false. )
      if (Lun_out > 0) then
         err= exfin ('ECCC NEXTGEM NWP MODEL - Development version',&
                     '','OK')
      end if
      
 34   format (i10.10)
!
!-------------------------------------------------------------------
!
      return
      end
