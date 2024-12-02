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

!**s/r init_component

      subroutine init_component (F_COMMs, F_nc)
      use iso_c_binding
      use gempi
      use clib_itf_mod
      use app
      use dcst
      use glb_ld
      use domains
      use HORgrid_options
      use dynkernel_options
      use svri_mod
      use svro_mod
      use path
      use ptopo
      use step_options
      use tdpack
      use wb_itf_mod
      use omp_timing
      use version
      implicit none

      integer, intent(IN) :: F_nc
      integer, intent(IN) :: F_COMMs(F_nc)
      
      include 'mpif.h'
      include 'gemdyn_version.inc'
#include <rmnlib_basics.hf>

      integer, external :: model_timeout_alarm, OMP_get_max_threads

      character(len=256) :: my_dir
      logical :: alongY_L
      integer :: i,ierr,mydomain
      integer :: npe_per
      integer :: GROUP_WORLD,gridpeer_group
      integer, dimension(:), allocatable :: process_ranks
      integer, parameter :: nargs=11, npos=0
      character(len=16) listec(nargs)
      character(len=2048) def(nargs), val(nargs)
!     
!--------------------------------------------------------------------
!
      COMM_world= F_COMMs(1)
      OUTs_COMM=F_COMMs(2) ; INs_COMM=F_COMMs(3)
      
      ierr = model_timeout_alarm (Step_alarm)
      
      ! List of non-positional optional arguments to the binary maingem
      ! with coresponding first and second default values (val & def)
      ! In MPI the arguments are provided using export CCARD_ARGS=
      listec = [ CHARACTER(LEN=16) :: &
                    'npex.', 'npey.', 'smtdyn.', 'smtphy.', 'ngrids.',&
                    'dom_start.', 'dom_end.', 'dom_last.', 'along_Y.',&
                    'input.', 'output.' ]
                    
      val= [ '1  ', '1  ', '0  ', '0  ', '1  ', '0  ', '0  ', '0  ', '.t.', '   ', '   ' ]
      def= [ '1  ', '1  ', '0  ', '0  ', '1  ', '0  ', '0  ', '0  ', '.f.', '   ', '   ' ]

      ! Obtain values of calling arguments (listec)
      call ccard (listec,def,val,nargs,npos)
      read (val(1),*) Ptopo_npex
      read (val(2),*) Ptopo_npey
      read (val(3),*) Ptopo_nthreads_dyn
      read (val(4),*) Ptopo_nthreads_phy
      read (val(5),*) Domains_ngrids
      read (val(6),*) Domains_deb
      read (val(7),*) Domains_fin
      read (val(8),*) Domains_last
      read (val(9),*) alongY_L
      Path_input_S = trim(val(10))
      Path_output_S= trim(val(11))

      ierr = clib_getenv ('PWD', Path_work_S)

      Domains_num= Domains_fin - Domains_deb + 1
      Ptopo_last_domain_L = (Domains_fin == Domains_last)
      Ptopo_alongY_L      = alongY_L .or. (Ptopo_npey == 1)

      Grd_yinyang_L = .false. ; Grd_yinyang_S = ''
      if ((Domains_ngrids < 1) .or. (Domains_ngrids > 2)) then
         write(6,'(/a,i4)') 'Unknown multigrid configuration: ABORT in init_component with Domains_ngrids=',Domains_ngrids
         stop -1
      else
         if (Domains_ngrids == 2) Grd_yinyang_L = .true.
      endif

      ierr = wb_put( 'model/Hgrid/is_yinyang',Grd_yinyang_L,&
                       WB_REWRITE_NONE+WB_IS_LOCAL )
      ierr = Lib_LogLevelNo(APP_LIBVGRID,APP_FATAL)  

      call MPI_COMM_size (COMM_world,gMPI_Wnumproc,ierr)
      call MPI_COMM_rank (COMM_world,gMPI_Wmyproc ,ierr)

      if (mod(gMPI_Wnumproc,Domains_num) /= 0) then
         if(gMPI_Wmyproc==0) write(6,'(/a,2i7)') 'MPI rank does NOT match requirement: RANK, Domains_num ',gMPI_Wnumproc,Domains_num
         call MPI_finalize(ierr)
         stop         
      endif

      npe_per = gMPI_Wnumproc / Domains_num
      mydomain= gMPI_Wmyproc  / npe_per + Domains_deb
      call MPI_Comm_split (COMM_world, mydomain, gMPI_Wmyproc, COMM_multigrid, ierr)
      call MPI_COMM_size (COMM_multigrid,gMPI_MGnumproc,ierr)
      call MPI_COMM_rank (COMM_multigrid,gMPI_MGmyproc ,ierr)
      npe_per = gMPI_MGnumproc / Domains_ngrids
      Ptopo_couleur= gMPI_MGmyproc  / npe_per
      call MPI_Comm_split (COMM_multigrid, Ptopo_couleur, gMPI_MGmyproc, COMM_grid, ierr)
      call MPI_COMM_size (COMM_grid,Ptopo_numproc,ierr)
      call MPI_COMM_rank (COMM_grid,Ptopo_myproc ,ierr)
      call MPI_COMM_GROUP(COMM_multigrid, GROUP_WORLD, ierr )
      if ( Domains_ngrids > 1) then
         allocate (process_ranks(0:Domains_ngrids-1))
         do i=0, Domains_ngrids-2
            process_ranks(i  )= Ptopo_myproc
            process_ranks(i+1)= process_ranks(i)+Ptopo_numproc
         end do
         call MPI_GROUP_INCL(GROUP_WORLD, Domains_ngrids, process_ranks, gridpeer_group, ierr)
         call MPI_COMM_CREATE(COMM_multigrid, gridpeer_group, COMM_gridpeers, ierr)
         deallocate (process_ranks)
      endif
      if (Ptopo_alongY_L) then
         Ptopo_mycol = Ptopo_myproc/Ptopo_npey
         Ptopo_myrow = mod(Ptopo_myproc, Ptopo_npey)
      else
         Ptopo_myrow = Ptopo_myproc/Ptopo_npex
         Ptopo_mycol = mod(Ptopo_myproc, Ptopo_npex)
      endif
      call MPI_Comm_split (COMM_grid, Ptopo_myrow, gMPI_Wmyproc, COMM_row, ierr)
      call MPI_Comm_split (COMM_grid, Ptopo_mycol, gMPI_Wmyproc, COMM_col, ierr)
      call MPI_COMM_size (COMM_row,num_in_row,ierr)
      call MPI_COMM_rank (COMM_row,me_in_row ,ierr)
      call MPI_COMM_size (COMM_col,num_in_col,ierr)
      call MPI_COMM_rank (COMM_col,me_in_col ,ierr)
      
      write(my_dir,'(a,i4.4)') 'cfg_',mydomain

      Path_input_S  = trim(Path_input_S ) // '/' // trim(my_dir)
      Path_work_S   = trim(Path_work_S  ) // '/' // trim(my_dir)
      Path_output_S = trim(Path_output_S) // '/' // trim(my_dir)

      ! cd to work directory
      ierr = clib_chdir (trim(Path_work_S))

      Dcst_rayt_8      = rayt_8
      Dcst_inv_rayt_8  = 1.d0 / rayt_8
      Dcst_omega_8     = omega_8
                  
      lun_out     = -1
      Lun_debug_L = .false.

      if (Ptopo_myproc == 0) then
         lun_out = output_unit
         ierr= exdb('ECCC NEXTGEM NWP MODEL - Development version',&
                    '','NON')
         ierr= clib_mkdir (Path_output_S)

         call  open_status_file3 (trim(Path_output_S)//'/status_MOD.dot')
         call write_status_file3 ('_status=ABORT' )
      endif

      call gtmg_init ( Ptopo_myproc, COMM_grid, 'MOD' )
      call gtmg_start ( 1, 'GEMDM', 0)

      ! Init color associated variables
      if (Grd_yinyang_L) then
         call MPI_COMM_size (COMM_multigrid,Ptopo_world_numproc,ierr)
         call MPI_COMM_rank (COMM_multigrid,Ptopo_world_myproc ,ierr)
         if (Ptopo_couleur == 0) Grd_yinyang_S = 'YIN'
         if (Ptopo_couleur == 1) Grd_yinyang_S = 'YAN'
         if (Ptopo_myproc  == 0) ierr= clib_mkdir (trim(Grd_yinyang_S))
         call MPI_barrier (COMM_grid, ierr)
         ierr= clib_chdir(trim(Grd_yinyang_S))
         Ptopo_ncolors = 2
      else
         Ptopo_couleur = 0
         Ptopo_ncolors = 1
         Ptopo_world_numproc = Ptopo_numproc
         Ptopo_world_myproc  = Ptopo_myproc
      end if
      
      ! Initialize OpenMP
      Ptopo_npeOpenMP = OMP_get_max_threads()
      if (Ptopo_nthreads_dyn < 1) Ptopo_nthreads_dyn=Ptopo_npeOpenMP
      if (Ptopo_nthreads_phy < 1) Ptopo_nthreads_phy=Ptopo_npeOpenMP
      call set_num_threads ( Ptopo_nthreads_dyn, 0 )

      if (Lun_out > 0) then
         write (Lun_out, 8255) Ptopo_npex, Ptopo_npey, Ptopo_alongY_L,&
                   Ptopo_npeOpenMP,Ptopo_nthreads_dyn, &
                   Ptopo_nthreads_phy
         write (Lun_out, 8256) trim(Path_work_S)
      end if

      call msg_set_can_write (Ptopo_myproc == 0)

      call pe_all_topo()

      ! Initialize local sub domain boundaries flags
      G_periodx = .false.
      G_periody = .false.

      l_west  = (0 == Ptopo_mycol)
      l_east  = (Ptopo_npex-1 == Ptopo_mycol)
      l_south = (0 == Ptopo_myrow)
      l_north = (Ptopo_npey-1 == Ptopo_myrow)

      north = 0 ; south = 0 ; east = 0 ; west = 0
      if (l_north) north = 1
      if (l_south) south = 1
      if (l_east ) east  = 1
      if (l_west ) west  = 1

 8255 format (/," MPI CONFIG (npex x npey): ",i4,' x ',i4,&
              8x,'Ptopo_alongY_L= ',l1/, &
              " OMP CONFIG (npeOpenMP x nthreads_dyn x nthreads_phy): ",&
              i4,' x ',i3,' x ',i3)
 8256 format (/," WORKING DIRECTORY:"/a/)
!
!--------------------------------------------------------------------
!
      return
      end subroutine init_component
