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
!---------------------------------- LICENCE END ----------------------

!*s/r spn_init - initialize spectral nudging profile, filter

      subroutine spn_init()
      use HORgrid_options
      use spn_options
      use glb_ld
      use glb_pil
      use cstv
      use lun
      use ver
      use ldnh
      use step_options
      use tdpack
      use ptopo
      use, intrinsic :: iso_c_binding
      use mpi_f08
      use gem_fft
      implicit none

      logical, external :: decomp
      integer, parameter :: lowest = 2
      integer :: minx, maxx, n, npartiel, n0
      integer :: k, err1, err2, tmdt
      real    :: t_turn, b_turn
      real(kind=real64) :: fft_norm_factor
      type(c_ptr) :: spn_c_ptr ! c_ptr used to allocate the spn_wrk array via MPI_Alloc_mem
      type(MPI_Comm) :: grid_comm_object
!
!----------------------------------------------------------------------
!

      ! Modification: May 2023, Christopher Subich, for standalone transpose module
      !
      ! The overall goal of the Spectral Nudging code is to take a field from the current simulation
      ! (a regional model) plus a field from an outer nesting simulation (typically a lower-resolution
      ! model, often global), and nudge the former towards the latter at the coarse scales.  The 
      ! underlying assumption is that this regional simulation has skill in the fine details, but the
      ! outer simulation is more likely to be correct over very large distances.

      ! Mathematically, this involves computing the difference between the two fields, applying a low-pass
      ! filter via FFT (cosine transform), then adding some of that difference back, with a level-dependent
      ! weight.

      ! Numerically, we want to apply this filter only inside the valid region of the grid (inside the
      ! piloting region).  We also begin with a grid decomposed over several MPI processes in an x-y grid,
      ! so we must also perform MPI-parallel array transposes in order to have full colums in X and Y for
      ! the FFT (DCT).

      ! If this is a global model or the nudging frequency is non-positive, we do not perform spectral
      ! nudging and can return
!      if ( Grd_yinyang_L .or. (Spn_freq<=0) ) return
      if ( Grd_yinyang_L .or. (.not.Spn_ON_L) ) return

      err1= 0 ; err2= 0

      ! Canonicalize the spectral nudging profile transition shape
      call low2up( Spn_trans_shape_S, Spn_trans_shape_S )

      if (Lun_out > 0) write(Lun_out,1000)
      
      ! if (Lun_out > 0) write(Lun_out,1002) ' Transpose 1===>2 for SPN:', &
      !            ' G_nk distributed on Ptopo_npex PEs', G_nk,Ptopo_npex

      ! As a reminder of GEM grid logic, the local tracer grid has a shape of
      ! arr(l_minx:l_maxx, l_miny:l_maxy, G_nk=l_nk).  The {min,max}{xy} values
      ! exist to accomodate the grid halo region, and the valid interior grid data
      ! exists from 1:l_ni and 1:l_nj in x and y respectively; this divides the global
      ! grid of G_ni x G_nj x G_nk points.

      ! However, this global grid also includes the piloting region, and we do not want
      ! to include this in the spectral nudging.  We want to find and remove this region
      ! for our purposes:

      ! Find the local extent of any piloting region.  This also gives the offset between
      ! the coordinates of the Spn_wrk array (grid-internal only) and the input array,
      ! which includes the piloting region
      if (l_west) then ! Is this process touching the west boundary?
         Spn_pil_w = pil_w
      else
         Spn_pil_w = 0
      end if
      if (l_east) then ! Is this process touching the east boundary?
         Spn_pil_e = pil_e
      else
         Spn_pil_e = 0
      end if
      if (l_south) then ! South boundary?
         Spn_pil_s = pil_s
      else
         Spn_pil_s = 0
      end if
      if (l_north) then ! North boundary?
         Spn_pil_n = pil_n
      else
         Spn_pil_n = 0
      end if

      ! The initial working grid is the z-global grid, which is the subset of the tracer grid that
      ! excludes the piloting region.  This gives local extents of:
      Spn_zgrid_lnx = max(l_ni - Spn_pil_w - Spn_pil_e,0)
      Spn_zgrid_lny = max(l_nj - Spn_pil_s - Spn_pil_n,0)
      Spn_zgrid_lnz = G_nk

      ! Allocate the working array via MPI_Alloc_mem.  This might offer slightly better performance with
      ! the MPI-RMA routines used for the transpose

      call MPI_Alloc_mem(Spn_zgrid_lnx * Spn_zgrid_lny * Spn_zgrid_lnz * c_sizeof(real(0,kind=REAL64)), & ! Allocation size
                         MPI_INFO_NULL, & ! MPI_Info parameter (null)
                         Spn_c_ptr)
      call c_f_pointer(Spn_c_ptr,Spn_wrk,[Spn_zgrid_lnx, Spn_zgrid_lny, Spn_zgrid_lnz])

      ! Use the transpose module to define the two transposes

      ! Set the MPI_VAL member of grid_comm_object, wrapping the integer handle in the type(MPI_Comm)
      ! opaque object.  See section 2.5.1 of the MPI Standard for details.
      grid_comm_object%MPI_VAL = COMM_grid

      ! zx_transpose is the first transpose, which takes a z-contiguous array and gives an x-contiguous array
      call transpose_create(descriptor = zx_transpose, & ! Transpose descriptor object
                            rank_i = Ptopo_mycol, & ! Process rank along the i (first) dimension
                            rank_j = Ptopo_myrow, & ! Process rank along the j (second) dimension
                            source = Spn_wrk,     & ! Source array for the transpose
                            comm = grid_comm_object)! MPI communicator (this grid)

      ! The destination of the first transpose is the x-contiguous grid
      Spn_xgrid => zx_transpose%dst_array

      ! xy_transpose is the second transpose, which takes an x-contiguous array and gives a y-contiguous array.
      ! Note that the transpose operator labels the dimensions as (i,j,k) from the perspective of the source
      ! array, so for this transpose i=y, j=x, and k=z.  Thus, the meaning of rank_i and rank_j switch from
      ! zx_transpose.
      call transpose_create(descriptor = xy_transpose, rank_i = Ptopo_myrow, &
                            rank_j = Ptopo_mycol, source = Spn_xgrid, comm=grid_comm_object)
      
      ! The transpose descriptor objects contain bounds information for the local grids, including their places
      ! within the global array.  We will need this information to compute the filter coefficients (in wavenumber
      ! space) and to apply the filter.
      !
      ! In the descriptor object, local bounds and extents are stored as arrays of size (numproc), so we must select
      ! the proper index based on Ptopo_mycol and Ptopo_myrow.  The latter MPI ranks are zero-based, so add one.

      ! for zx_transpose, i->x, j->y, k->z; for xy_transpose, i->y, j->z, k->x
      ! Spn_zgrid has order (x,y,z)
      Spn_zgrid_llbx = zx_transpose%llbi(Ptopo_mycol+1)
      ! Spn_zgrid_lnx  = zx_transpose%lni(Ptopo_mycol+1) ! Matches value calculated above
      ! The y dimension is invariant for the zx_transpose, so that descriptor does not store the full set of
      ! global information.  We can obtain it from the xy_transpose descriptor, however
      Spn_zgrid_llby = xy_transpose%llbi(Ptopo_myrow+1) 
      ! Spn_zgrid_lny  = xy_transpose%lni(Ptopo_myrow+1) ! Matches value calculated above
      Spn_zgrid_llbz = 1
      ! Spn_zgrid_lnz  = zx_transpose%gnk ! Matches Gnk, used above

      ! Spn_xgrid has order (y,z,x)
      Spn_xgrid_llby = xy_transpose%llbi(Ptopo_myrow+1) ! Matches the y-dimension of Spn_zgrid
      Spn_xgrid_lny  = xy_transpose%lni(Ptopo_myrow+1)
      Spn_xgrid_llbz = zx_transpose%llbk(Ptopo_mycol+1)
      Spn_xgrid_lnz  = zx_transpose%lnk(Ptopo_mycol+1)
      Spn_xgrid_llbx = 1
      Spn_xgrid_lnx  = zx_transpose%gni

      ! Spn_ygrid has order (z,x,y)
      Spn_ygrid_llbz = zx_transpose%llbk(Ptopo_mycol+1) ! Matches the z-dimension of Spn_xgrid
      Spn_ygrid_lnz  = zx_transpose%lnk(Ptopo_mycol+1)
      Spn_ygrid_llbx = xy_transpose%llbk(Ptopo_myrow+1)
      Spn_ygrid_lnx  = xy_transpose%lnk(Ptopo_myrow+1)
      Spn_ygrid_llby = 1
      Spn_ygrid_lny  = xy_transpose%gni
      
      ! ! Compute the first parallel transpose, which takes the split-x, split-y, global-z
      ! ! local grid and gives a split-y, split-z, global-x grid.  
      ! if (.not. decomp (G_nk, minx, maxx, n, npartiel, 0, n0, &
      !           .true., .true., Ptopo_npex, -1, .false., 3 )) err1 = -1
      
      ! Set legacy grid values: Spn_12s{etc} are values relating to the split of
      ! G_nk vertical levels across the Ptopo_npex columns (indexed by Ptopo_mycol)
      ! In the context of the new transforms, these are the bounds and extents on
      ! the z-dimension (second) of the x-contiguous array
      ! Spn_12smin = 1             ! Index of the local lower bound
      ! Spn_12smax = Spn_xgrid_lnz ! Index of the local upper bound
      ! Spn_12sn   = Spn_xgrid_lnz ! Extent
      ! Spn_12smin = minx
      ! Spn_12smax = maxx
      ! Spn_12sn   = n
      
      ! if (Lun_out > 0) write(Lun_out,1002) ' Transpose 2===>2 for SPN:', &
      !            ' G_ni distributed on Ptopo_npey PEs', G_ni,Ptopo_npey
      if (Lun_out > 0) then 
         write(Lun_out,'(" Transpose 1 (zx) for SPN: ",I0," (",I0,") G_nk points distributed over ",I0," (Ptopo_npex) PEs")') &
                           Spn_zgrid_lnz, G_nk, Ptopo_npex
         write(Lun_out,'(" Transpose 2 (xy) for SPN: ",I0," (",I0,"-",I0,") G_ni points distributed over ",I0," (Ptopo_npey) PEs")') &
                           Spn_xgrid_lnx, G_ni, 2*Grd_extension, Ptopo_npey
      end if

      ! if (.not. decomp (G_ni, minx, maxx, n, npartiel, 0, n0, &
      !           .false., .true., Ptopo_npex, lowest, .false., 0 )) err1 = -1

      ! Set legacy grid values: Spn_22{etc} values relate to the split of G_ni
      ! horizontal levels across the Ptopo_npey rows (indexed by Ptopo_myrow)
      ! For the new transform, these are the bounds and extents on the x-dimension
      ! (second) of the y-contiguous array
      ! Spn_22min = 1              ! Local lower bound
      ! Spn_22max = Spn_ygrid_lnx  ! Local upper bound
      ! Spn_22n   = Spn_ygrid_lnx  ! Extent
      ! Spn_22n0  = Spn_ygrid_llbx ! Index of the lower bound along the global array
      ! Spn_22min = minx
      ! Spn_22max = maxx
      ! Spn_22n   = n
      ! Spn_22n0  = n0

      ! The pilot region is now removed prior to the transposes, so the 'pilot region'
      ! has extent zero.
      ! Spn_22pil_w= 0 ;  Spn_22pil_e= 0
      ! if (Spn_22n0==1)              Spn_22pil_w= Grd_extension
      ! if (Spn_22n0+Spn_22n-1==G_ni) Spn_22pil_e= Grd_extension

      ! if (err1<0) goto 999
      
      allocate ( prof(G_nk), & ! Vertical weight profile
                 Spn_flt(Spn_ygrid_lnx,Spn_ygrid_lny)) ! 2D spectral weights, valid on y-contiguous grid

      ! Initialize arrays
      prof = 0
      Spn_xgrid = 0
      Spn_ygrid = 0
      Spn_wrk = 0
      Spn_flt = 0
                 
      ! allocate ( prof(G_nk), &
      !    Spn_fft(ldnh_maxy ,Spn_12smax,G_ni+2+Ptopo_npex),& ! Order (y,z,x)
      !    Spn_fdg(Spn_12smax,Spn_22max ,G_nj  +Ptopo_npey),& ! Order (z,x,y)
      !    Spn_wrk(ldnh_maxx,ldnh_maxy,l_nk) ) ! Order (x,y,z)
      ! prof=0. ; Spn_fft= 0. ; Spn_fdg= 0. ; Spn_wrk= 0.
     
      ! Spn_njnh  = ldnh_maxy-ldnh_miny+1
      ! Spn_nk12  = Spn_12smax-Spn_12smin+1
      ! Spn_ni22  = Spn_22max-Spn_22min+1

      tmdt      = int(Cstv_dt_8) ! Inverse delta-t
      Spn_interval = Spn_freq/tmdt ! Interval (#steps) between applications of spectral nudging
      Spn_interval = max(1,Spn_interval) ! Bounds check: ensure at the target interval isn't 0
      Spn_ws = Step_nesdt/tmdt ! Ratio between delta t and delta t of nested field
      Spn_weight= 1.0 ! Temporal weight (intialized here, overwritten elsewhere)
      
      ! Calculate horizontal filter coefficients in terms of spectral coefficients.
      call spn_calfiltre ()

      ! Calculate the vertical profle weights
      t_turn= max( Spn_up_const_lev,Ver_hyb%m(  1 ) ) ! Level above which the weight is a constant 1
      b_turn= min( Spn_start_lev   ,Ver_hyb%m(G_nk) ) ! Level below which the weight is a constant 0

      if (Spn_trans_shape_S == 'COS2' ) then
         ! Transition between 0 and 1 with a cos^2 shape in eta
         do k=1,G_nk
            if (Ver_hyb%m(k) <= b_turn .and. Ver_hyb%m(k) >= t_turn) then
               prof(k) = cos(pi_8-pi_8*(b_turn-Ver_hyb%m(k))/(b_turn-t_turn))
            elseif (Ver_hyb%m(k) < t_turn) then
               prof(k)=1.
            else
               prof(k)=0.
            end if
            prof(k) = prof(k)*prof(k)
         end do

      elseif (Spn_trans_shape_S == 'LINEAR' ) then
         ! Transition linearly between 0 and 1, as a function of eta
         do k=1,G_nk
            if (Ver_hyb%m(k) <= b_turn .and. Ver_hyb%m(k) >= t_turn) then
               prof(k) =  (b_turn-Ver_hyb%m(k))/(b_turn-t_turn)
            elseif (Ver_hyb%m(k) < t_turn) then
               prof(k)=1.
            else
               prof(k)=0.
            end if
         end do

      else
         ! Invalid shape specified
         err2 = -1
      end if

      999  call gem_error ( min(err1,err2),'spn_init',&
               'Wrong choice for Spn_trans_shape_S or transpose problems')

      ! Compute FFT plans
      call make_r2r_dft_plan(fft_x_forward, & ! Plan variable
                             Spn_xgrid, Spn_xgrid, & ! Source and destination arrays, matching for in-place transform
                             3, 'QCOS', 1) ! Logical dimension of transform, transform type, direction (forward)
      call make_r2r_dft_plan(fft_x_reverse, Spn_xgrid, Spn_xgrid, 3, 'QCOS', 0)
      call make_r2r_dft_plan(fft_y_forward, Spn_ygrid, Spn_ygrid, 3, 'QCOS', 1)
      call make_r2r_dft_plan(fft_y_reverse, Spn_ygrid, Spn_ygrid, 3, 'QCOS', 0)

      ! The gem_fft_mod driver does not automatically normalize the transforms:
      ! f = norm_factor * T^(-1)(T(f)).
      fft_norm_factor = get_dft_norm_factor(Spn_xgrid_lnx,'QCOS') * get_dft_norm_factor(Spn_ygrid_lny,'QCOS')
      
      ! Scale the vertical level weights by Spn_relax_hours, setting the time constant
      ! of (unreduced) nudging to that value, and also apply the FFT normalization factor.
      do k=1,G_nk
         prof(k) = prof(k) * Cstv_dt_8/(Spn_relax_hours*3600.) * fft_norm_factor
      end do

!      Spn_ON_L= .true.

 1000 format (/' SPN_INIT: Initialization of spectral nudging'/)
 1002 format (a/a45,i6,' /',i5)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_init
