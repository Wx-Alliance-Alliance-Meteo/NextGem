!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2023 - Division de Recherche en Prevision Numerique
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


! Created -- April 2023, Christopher Subich, based on 
! https://gitlab.science.gc.ca/csu001/transpose_mod
! commit a73fa0ce

! This module implements a generic parallel transposition of a three-dimensional
! array, divided approximately evenly among a logically-2D grid of processors.

! The typical 3D array in GEM is divided between (pi*pj) processors along the x (i) 
! and y (j) dimension, but it is contiguous in z (k).  For a global array of 
! (gnx, gny, gnz) elements, each local processor holds approximately 
! (gnx/pi, gny/pj, gnz) elements.

! The typical use case of this transpose is to rearrange the array for a Fourier
! transform along the i or j dimension (or both, in sequence).  For these transforms,
! each process needs to own a set of complete i- or j-rows, which requires the parallel
! transpose.

! The logical model of this process is:

! 0: (pi, pj) processors each owning arrays of size (gnx/pi, gny/pj, gnz)
! 1: (pi, pj) processors each owning arrays of size (gnx, gny/pj, gnz/pi)
!    (Fourier transform along i)
! 2: (pi, pj) processors each owning arrays of size (gnx/pj, gny, gnz/pi)
!    (Fourier transform along j)
! Process the resulting array / invert Fourier transforms and transposes

! Notice that in both of these transposes, one of the "divided" dimensions
! remains constant.  In the first transpose, the second logical dimension 
! (size gny/pj) is unchanged, and in the second transpose the third logical
! dimensin (size gnz/ni) is unchanged. 

! This limits the scope of communication required.  For the first transform,
! the full communicator of size pi*pj splits into pj separate communicators,
! each of size pi.  Likewise, the second transform splits the communicator
! into pi communicators, of size pj

! The symmetry of the transposes becomes even more clear if we allow ourselves
! to permute the arrays in memory during the parallel transpose:

! 0: (pi, pj) processors owning arrays of size (gnx/pi, gny/pj, gnz) [0,1,2]
! 1: (pj, pi) processors owning arrays of size (gny/pj, gnz/pi, gnx) [1,2,0]
! 2: (pi, pj) processors owning arrays of size (gny/pj, gnx/pj, gny) [2,0,1]

! That is, each transpose shifts the array layout by one dimension.  The global
! dimension starts as the third dimension and becomes the second (being split), 
! the preserved dimension starts as the second and becomes the first (unchanged), 
! and the split dimension starts as the first and becomes the third (becoming 
! contiguous).  

! With this restriction on input data ordering, we only need one transpose operator.

! To simplify the implementation, we will also make one other restrictive assumption:
! transposes happen from and to fixed memory locations.  This restriction lets us
! use one-sided MPI communication (MPI_Put) and associate both the source and
! destination with an MPI window.  This is not a major restriction for the transpose-
! destination arrays, since FFTW also likes to have fixed operating locations, but it
! is something to keep note of for the source arrays. 

! This "universal transpose" formulation has one nice conceptual advantage from the MPI
! perspective: processes will send and receive contiguous chunks of data to each other.
! Each local chunk being sent needs a local transposition, however, and for performance
! reasons we will handle that with a set of workspace arrays.  In theory the MPI derived
! datatype system can construct a transposed type, but in practice the performance is
! abysmal.  

! Notation note: the indices i, j, and k will refer to the first, second, and third
! indices from the perspective of the forward transpose.  That is, the initial array
! takes the form source(lni,lnj,gnk) and transposes to dest(lnj,lnk,gni).  The l prefix
! denotes a local extent, and the g prefix denotes a global extent.  lnj is invariant,
! and the distinction between lni and lnk should be clear from the forward/reverse context.

module transpose
    use iso_c_binding
    use, intrinsic :: iso_fortran_env
    use mpi_f08 ! MPI module, required per modern MPI standard
    implicit none

    ! Derived type to collect the information for a single transpose
    type transpose_descriptor
        ! Pointers to the source and destination arrays of the transpose
        ! src_array has ordering (i,j,k), and dst_array has ordering (j,k,i)
        real(kind=real64), dimension(:,:,:), contiguous, pointer :: src_array, dst_array

        ! Array sizes
        integer :: gni, gnk ! Global extents along i and k 
        integer :: lnj      ! Local extent along j (invariant during transpose)
        ! Local extents along i and k, separate per processor in row
        integer, dimension(:), allocatable :: lni, lnk

        ! Lower bounds along i and k, for index calculation
        integer, dimension(:), allocatable :: llbi, llbk

        ! Segregate the private MPI objects, so that code simply using this module does not
        ! accidentally break internal consistency, e.g. by starting a communications epoch
        ! on an MPI window

        ! Communicator for the transpose
        type(MPI_Comm), private :: communicator = MPI_COMM_NULL

        ! MPI_Group, for managing PSCW synchronization
        type(MPI_Group), private :: group

        ! Communicator information, for later reference
        integer, private :: numproc, rank

        ! MPI windows for the source and destination arrays.  The destination window is
        ! used for the MPI_Put operations corresponding to a forward transpose, and the
        ! source window is used for the reverse transpose
        type(MPI_Win), private :: src_window = MPI_WIN_NULL, dst_window = MPI_WIN_NULL

        ! Index into the module-level open_windows 
        integer, private :: window_index = 0

        ! MPI datatype for the transpose
        type(MPI_Datatype), private :: mpitype ! Scalar with respect to the source array

        ! Array for in-flight MPI requests.  This is held in the descriptor object rather
        ! than local to the transpose routines because the descriptor object is shared among
        ! OpenMP threads, whereas stack-local variables in the transpose subroutines are
        ! thread-local
        type(MPI_Request), private, dimension(:), allocatable :: mpirequests

        ! Workspace variables.  src_workspace matches the ordering of src_array, and
        ! dst_workspace matches dst_array.  These are used to perform the in-memory
        ! tranapose of data in flight, before sending it via MPI_Put directly to
        ! the destination process
        type(workspace), dimension(:), allocatable, private :: src_workspace, dst_workspace

        ! Initialized flag, used for synchronization in an OpenMP context to prevent multiple
        ! creation or destruction of a transpose descriptor
        logical, private:: initialized=.false.
    end type

    ! Derived type to wrap a workspace pointer.  In Fortran, it is impossible to directly
    ! refer to arrays of pointers.  However, it is easy to refer to arrays of objects that
    ! contain just a single pointer
    type workspace
        real(kind=real64), dimension(:,:,:), contiguous, pointer :: arr
    end type

    !! Control over MPI synchronization styles for the parallel transpose
    !! In testing with 2400 processes (30 nodes), distributed 75x32 for an 1801x601x83 grid,
    !! SYNC_LOCK very slightly outperformed SYNC_FENCE, which outperformed SYNC_PSCW.

    !! Note: as described in https://gitlab.science.gc.ca/hpc/support/issues/742, PSCW 
    !! seems to have some problems on robert/underhill, so it should be used with caution
    !! and considerable testing.

    integer, parameter :: SYNC_LOCK = 1,  & !! Synchronize via MPI_Win_lock and MPI_Win_unlock (passive-target)
                          SYNC_FENCE = 2, & !! Synchronize via MPI_Win_fence
                          SYNC_PSCW = 3,  & !! Synchronize via MPI_Win_(post/start/complete/wait)
                          SYNC_FLUSH = 4    !! Synchronize with long-lived locks and MPI_Win_flush (also passive-target)
    integer, parameter :: SYNC_STYLE = SYNC_LOCK

    !! Parameters relating to the cleanup of MPI_Win objects during the MPI_Finalize call
    logical, private, target :: callbacks_initialized = .false.
    integer, private :: comm_keyval

    ! Record of open MPI_Win objects, for cleanup at program exit.  This array is
    ! append-only, and it is extended as necessary by transpose_create.  The implied
    ! performance assumption is that transpose descriptor objects will be created and
    ! destroyed rarely
    type(MPI_Win), dimension(:), pointer, private :: open_windows
    integer, private :: window_count = 0

    !! Set implementation details as private
    private :: workspace !! Workspace pointer-wrapper
    private :: parallel_equal, transpose_create_impl, sync_open, sync_close, transpose_clean_windows !! Internal subroutines

contains

    function parallel_equal(val, comm)
        !! Helper function to check whether a specified nonnegative value is equal on all processes in a
        !! communicator
        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        integer, intent(in) :: val !! value to check
        type(MPI_Comm), intent(in) :: comm !! communicator
        logical parallel_equal 
        integer, dimension(2) :: val_array
        integer :: ierr ! MPI error flag

        val_array(1) = val
        val_array(2) = -val

        call MPI_Allreduce(MPI_IN_PLACE,val_array, & ! source (unused), dest
                           2, MPI_INT, MPI_MAX, comm,ierr) ! size, type, op, comm
        if (val_array(1) == -val_array(2)) then
            parallel_equal = .true.
        else
            parallel_equal = .false.
        endif

        return
    end function

    subroutine transpose_create(descriptor, rank_i, rank_j, source, dest, comm, lni, lnj, lnk)
        !! Create a transpose descriptor, given processor rank information and either a source
        !! and destination array or information about the local array size.  (Wrapper for OpenMP
        !! consistency)        
        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        !use omp_lib
        implicit none

        ! Note: descriptor is intent(inout) for OpenMP reasons.  With just intent(out), threads might
        ! separately initialize the shared descriptor when this function is called, clobbering the
        ! actual initialization.
        type(transpose_descriptor), intent(inout) :: descriptor !! Generated array transpose descriptor
        integer, intent(in) :: rank_i, rank_j 
            !! Process location within the implied global 2D process grid.  The source array of
            !! process (rank_i, rank_j) covers approximately 
            !! ((rank_i)*gni/Pi : (rank_i+1)*gni/Pi, (rank_j)*gnj/Pj : (rank_j+1)*gnj/Pj, :)
            !! @warning Base 0 addressing of process ranks, as in the MPI library

        REAL(kind=real64), dimension(:,:,:), pointer, contiguous, optional, intent(in) :: source, dest
            !! Source and destination arrays.  Optional arguments.  If specified, these arrays
            !! will be used for the source and destination arrays of the transpose (respectively),
            !! and the arrays will be exposed to one-sided communication through a window creation.
            !! If one or both are not specified, the missing arrays will be allocated dynamically.

        type(MPI_Comm), optional, intent(in) :: comm
            !! Global MPI communicator for the transpose; rank_i and rank_j are with respect to this
            !! communicator.  If not specified, default to MPI_COMM_WORLD.

        integer, optional, intent(in) :: lni, lnj, lnk
            !! If neither source nor dest are specified, (lni, lnj, lnk) specifies the local sizes of
            !! the source and destination arrays to be allocated.  The source array will have shape
            !! (lni, lnj, gnk), and the destination array will have shape (lnj, lnk, gni), with the
            !! global i- and k-sizes implicitly computed by summing the lni and lnk values over the
            !! respective communicators.

        ! This function can be called in an OpenMP context, inside a parallel region.  In turn,
        ! the transpose-description creation calls MPI routines.  We need to ensure that:
        ! *) A transpose is created only once, and
        ! *) MPI thread limits are respected.

        ! In turn, there are two possibilities for how this function is being called:
        ! *) We could be creating many transposes, if each thread or parallel team has its
        !    own transpose descriptor
        ! *) We could be creating one single transpose, if each thread is referring to a shared
        !    (global, module-level, or explicitly shared by the parallel region) descriptor.

        ! It is okay to serialize transpose creation, but we must not break MPI's rules or double-create
        ! a transpose.

        logical :: exit_loop 
        integer :: ierr, rank, mpi_thread_type
        logical :: mainthread 
        call MPI_Query_thread(mpi_thread_type)

        call MPI_Is_thread_main(mainthread,ierr)

        call MPI_Comm_rank(MPI_Comm_world,rank)

        if (mpi_thread_type == MPI_THREAD_FUNNELED .and. .not. mainthread) then
            ! This thread is not allowed to make MPI calls.  The descriptor _must_ be shared
            ! with the master thread, but we cannot check that.  We must wait until the descriptor
            ! is created by the master thread before returning.  If the master thread is not creating
            ! this descriptor, this loop will hang.
            exit_loop = .false.
            do while (.not. exit_loop)
                !$omp critical (transpose_mod)
                ! Note that the critical section acts as a memory barrier, so descriptor%initialized
                ! should be fetched from memory with each iteration.  Otherwise, the compiler could
                ! load the value into a register and never update it.
                exit_loop = descriptor%initialized
                !$omp end critical (transpose_mod)
            end do
        else 
            ! This thread may make MPI calls, so it is permitted to create the descriptor.
            ! We do not want to double-create the descriptor if it is already initialized,
            ! otherwise we can call the implementation function to create the descriptor.

            ! The first thread to reach this critical section will "win" and create the transpose
            ! descriptor.  If many threads are trying to create one single transpose, then
            ! subsequent threads will see that the descriptor has been initialized and exit.

            ! If many threads are trying to create many transposes (one each), then this critical
            ! section will protect shared data.  However, it does not ensure synchronization between
            ! processes (i.e. process 1 thread 1 might talk to process 2 thread 5); that must be
            ! provided externally.  I recommend not letting this situation happen, it would be very
            ! confusing.
            !$omp critical (transpose_mod)
            if (descriptor%initialized .eqv. .false.) then
                ! Call the implementing subroutine
                !write(0,*) 'rank', rank, 'creating descriptor on thread', omp_get_thread_num()
                call transpose_create_impl(descriptor, rank_i, rank_j, source, dest, comm, lni, lnj, lnk)
            end if
            ! Else descriptor%initialized == .true., so nothing to do.
            !$omp end critical (transpose_mod)
        end if
    end subroutine


    subroutine transpose_create_impl(descriptor, rank_i, rank_j, source, dest, comm, lni, lnj, lnk)
        !! Create a transpose descriptor, given processor rank information and either a source
        !! and destination array or information about the local array size.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        
        type(transpose_descriptor), intent(out) :: descriptor !! Generated array transpose descriptor
        integer, intent(in) :: rank_i, rank_j 
            !! Process location within the implied global 2D process grid.  The source array of
            !! process (rank_i, rank_j) covers approximately 
            !! ((rank_i)*gni/Pi : (rank_i+1)*gni/Pi, (rank_j)*gnj/Pj : (rank_j+1)*gnj/Pj, :)
            !! @warning Base 0 addressing of process ranks, as in the MPI library

        REAL(kind=real64), dimension(:,:,:), pointer, contiguous, optional, intent(in) :: source, dest
            !! Source and destination arrays.  Optional arguments.  If specified, these arrays
            !! will be used for the source and destination arrays of the transpose (respectively),
            !! and the arrays will be exposed to one-sided communication through a window creation.
            !! If one or both are not specified, the missing arrays will be allocated dynamically.

        type(MPI_Comm), optional, intent(in) :: comm
            !! Global MPI communicator for the transpose; rank_i and rank_j are with respect to this
            !! communicator.  If not specified, default to MPI_COMM_WORLD.

        integer, optional, intent(in) :: lni, lnj, lnk
            !! If neither source nor dest are specified, (lni, lnj, lnk) specifies the local sizes of
            !! the source and destination arrays to be allocated.  The source array will have shape
            !! (lni, lnj, gnk), and the destination array will have shape (lnj, lnk, gni), with the
            !! global i- and k-sizes implicitly computed by summing the lni and lnk values over the
            !! respective communicators.

        !! local variables, some of which will be copied into the transpose descriptor
        integer :: nproc_i ! number of processors along i
        integer :: myrank_i ! My local rank in i, after the communicator split
        integer :: gni, gnk ! Global sizes along i, k

        integer :: ii ! loop index variable
        integer :: ierr ! MPI error flag

        type(MPI_Comm) :: global_communicator, lcomm
        !type(MPI_Info) :: window_info

        !integer, dimension(2) :: comm_int ! Workspace for MPI operations, integer type
        integer :: temp ! temporary value
        logical :: ltemp ! temporary logical value

        type(c_ptr) :: allocated_window

        if (.not. present(comm)) then
            global_communicator = MPI_COMM_WORLD
        else
            global_communicator = comm
        endif

        ! Split the global communicator into separate row communicators.  Since the forward transpose
        ! splits z and consolidates i, the row needs to cover all processes that share the same j-rank.
        ! That forms the 'colour' argument to mpi_comm_split; rank_i is then the key to give the proper
        ! ordering inside the new communicator.
        call MPI_Comm_split(global_communicator, rank_j, rank_i, lcomm,ierr)
        descriptor%communicator = lcomm
        call MPI_Comm_group(lcomm, descriptor%group,ierr)

   
        ! In order to clean up MPI_Windows on program exit regardless of how that occurs, we want to 
        ! register a callback with the MPI implementation.  Per the MPI specification, callbacks on
        ! MPI_COMM_SELF will be called by MPI_Finalize.
        if (.not. callbacks_initialized) then
            ! comm_keyval = int(loc(callbacks_initialized),kind=MPI_ADDRESS_KIND)
            call MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, transpose_clean_windows, & 
                                        comm_keyval,int(LOC(callbacks_initialized),kind=MPI_ADDRESS_KIND), ierr)

            call mpi_comm_set_attr(MPI_COMM_SELF,comm_keyval,int(0,kind=MPI_ADDRESS_KIND))
            callbacks_initialized = .true.
        end if

        ! Get the number of processes and local rank along the row communicator
        call MPI_Comm_size(descriptor%communicator,nproc_i,ierr)
        call MPI_Comm_rank(descriptor%communicator,myrank_i,ierr)
        descriptor%numproc = nproc_i
        descriptor%rank = myrank_i
        allocate(descriptor%mpirequests(nproc_i))

        if (myrank_i /= rank_i) then
            write(0,*) 'Error in MPI split, process', rank_i, rank_j, 'erroneous rank', myrank_i
            stop 'error'
        end if

        ! Perform very basic sanity checks: make sure that either all or none of the processors
        ! in this communicator specify source and dest.  Otherwise, the grid extent computations
        ! will deadlock / give invalid results when the collective operations are not executed
        ! in sync.

        ! use temp to check flags:
        ! temp = (source_present + 2*dest_present)
        temp = 0
        if (present(source)) then
            temp = 1 
        end if
        if (present(dest)) then
            temp = temp + 2
        end if
        if (.not. parallel_equal(temp,descriptor%communicator)) then
            write(0,*) 'Inconsistent transpose_create call ', rank_i, rank_j, present(source), present(dest)
            stop 'error'
        end if

        !! From here, use the provided information to compute the local extents of the transposing arrays
        !! (source and destination), as well as the global extent along i and k (dimensions 1 and 3 of the
        !! source array).  

        !! This calculation is a bit complicated because there are several possible combinations of admissible/
        !! reasonable parameters.  Source and dest might be present or absent, and lni/lnj/lnk are all optional.


        ! Allocate the extent and lower bound arrays in the descriptor
        allocate(descriptor%lni(nproc_i))
        allocate(descriptor%lnk(nproc_i))
        allocate(descriptor%llbi(nproc_i))
        allocate(descriptor%llbk(nproc_i))

        ! Allocate the arrays of workspace pointers
        allocate(descriptor%src_workspace(nproc_i))
        allocate(descriptor%dst_workspace(nproc_i))

        ! Set signal values for array extents, indicating that they have not yet been
        ! defined by a provided array
        gni = 0
        gnk = 0
        descriptor%lni = 0
        descriptor%lnk = 0
        descriptor%llbi = 0
        descriptor%llbk = 0
        descriptor%lnj = 0

        ! Take any grid-defining parameters
        if (present(lni)) then
            descriptor%lni(myrank_i+1) = lni
        end if
        if (present(lnj)) then
            descriptor%lnj = lnj
        end if
        if (present(lnk)) then
            descriptor%lnk(myrank_i+1) = lnk
        end if

        ! Check if lnjs are equal so far
        if (.not. parallel_equal(descriptor%lnj,lcomm)) then
            write(0,*) 'inconsistent parameter lnj', rank_i, rank_j, 'has', descriptor%lnj
            stop 'error'
        end if

        ! If the source array is present, assign it to the descriptor.  The source array provides information
        ! about lnj, gnk, and this process's lni, which later can be aggregated to gni.  Check to ensure
        ! that lnj and gnk are consistent among all processes in the row and with any previous definition.
        if (present(source)) then
            descriptor%src_array => source ! Bind the src_array pointer
            gnk = size(source,3)
            descriptor%gnk = gnk
            if (.not. parallel_equal(size(source,2),lcomm)) then
                write(0,*) 'inconsistent source j-shape', rank_i, rank_j, 'has', size(source,2)
                stop 'error'
            end if
            if (descriptor%lnj == 0) then
                descriptor%lnj = size(source,2)
            elseif (descriptor%lnj /= size(source,2)) then
                write(0,*) 'inconsistent source j-shape', rank_i, rank_j, 'has', size(source,2), 'needs', descriptor%lnj
                stop 'error'
            end if

            ! Infer the lni size from the local array, and write it to the descriptor
            ! in the corrent place. 
            
            ! Check whether the local array size is consistent with lni, if previously-specified
            if (descriptor%lni(rank_i+1) == 0) then
                ! Match: no previous specification, so the array size wins
                temp = 1
            elseif (descriptor%lni(rank_i+1) == size(source,1)) then
                ! Match: previous specification matches array size
                temp = 1
            else
                ! Mismatch: specified lni inconsistent with array size
                temp = 0
            endif

            ltemp = parallel_equal(temp,lcomm) ! Ensure collective call always happens
            if (.not. ltemp .or. temp == 0) then
                write(0,*) 'inconsistent source i-shape', rank_i, rank_j, 'array has', size(source,1), 'parameter', descriptor%lni(myrank_i+1)
                stop 'error'
            end if
            ! Write the array size back to the appropriate lni index; this either matches the previous
            ! value or overwrites a 0
            descriptor%lni(rank_i+1) = size(source,1) 

            ! Since the source array has global-k extent, ensure that all processes in this row
            ! do in fact see the same k-size
            if (.not. parallel_equal(gnk,lcomm)) then
                write(0,*) 'inconsistent source k-shape', rank_i, rank_j, 'has', gnk
                stop 'error'
            end if
        end if

        ! Gather all of the lni sizes specified so far.
        call MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, & ! Source buffer(=dest), count, type (ignored)
                           descriptor%lni, 1, MPI_INT, & ! Dest buffer, count, type
                           descriptor%communicator,ierr)
        ! The global i-extent is the sum of the local extents.
        gni = sum(descriptor%lni)
        descriptor%gni = gni
        ! If lni has not yet been specified by any process, this will just re-write the signal value
        ! of 0.

        ! If the destination array is present, it provides lnk values along with gni and lnj.  Assign
        ! the destination array to the descriptor, then check this array information for consistency
        ! with both any previous setting (from parameters or source) and across processes.  As a reminder,
        ! the destination array is after the transpose and has shape (lnj, lnk, gni)
        if (present(dest)) then
            descriptor%dst_array => dest
            descriptor%lnk(myrank_i+1) = size(dest,2)

            ! Check equality of gni implied by dest array
            if (.not. parallel_equal(size(dest,3),lcomm)) then
                write(0,*) 'inconsistent dest i-shape (dim 3)', rank_i, rank_j, 'has', size(dest,3)
                stop 'error'
            end if

            if (gni == 0) then
                ! Set gni based on the destination array
                gni = size(dest,3)
                descriptor%gni = gni
            elseif (gni /= size(dest,3)) then
                ! The i-distribution has been specified, but it's inconsistent with this array
                write(0,*) 'inconsistent dest-i shape (dim 3)', rank_i, rank_j, 'has', size(dest,3), 'needs', gni
                stop 'error'
            endif

            ! Check equality of lnj implied by dest array
            if (.not. parallel_equal(size(dest,1),lcomm)) then
                write(0,*) 'inconsistent dest j-shape (dim 1)', rank_i, rank_j, 'has', size(dest,1)
                stop 'error'
            end if
            if (descriptor%lnj == 0) then
                ! Write the implied lnj to the descriptor
                descriptor%lnj = size(dest,1)
            elseif (descriptor%lnj /= size(dest,1)) then
                write(0,*) 'inconsistent dest j-shape (dim 1)', rank_i, rank_j, 'has', size(dest,1), 'needs', descriptor%lnj
                stop 'error'
            end if

            ! Check whether the destination array is compatible with lnk, if previously specified
            if (descriptor%lnk(rank_i+1) == 0) then
                ! Match: no previous specifiation
                temp = 1
            elseif (descriptor%lnk(rank_i+1) == size(dest,2)) then
                ! Match: definitions match
                temp = 1
            else
                ! No match, inconsistent definitions
                temp = 0
            end if

            ltemp = parallel_equal(temp,lcomm)
            if (.not. ltemp .or. temp == 0) then
                write(0,*) 'inconsistent dest k-shape (dim 2)', rank_i, rank_j, 'has', size(dest,2), 'needs', descriptor%lnk(rank_i+1)
                stop 'error'
            end if
            ! Write the lnk implied by the destination array back to the descriptor
            descriptor%lnk(rank_i+1) = size(dest,2)
        end if

        ! Gather all of the lnk sizes specified so far.
        call MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, & ! Source buffer(=dest), count, type (ignored)
                           descriptor%lnk, 1, MPI_INT, & ! Dest buffer, count, type
                           descriptor%communicator,ierr)
        ! Check for consistency with gnk
        if (sum(descriptor%lnk) /= 0) then
            ! no lnk has been specified, so we need things to be consistent
            if (descriptor%gnk == 0) then
                ! If descriptor%gnk is 0, the lnk values define gnk
                descriptor%gnk = sum(descriptor%lnk)
            elseif (descriptor%gnk /= sum(descriptor%lnk)) then
                write(0,*) 'Inconsistent global k-extent', rank_i, rank_j, 'has local', descriptor%lnk(rank_i+1), 'total needs', descriptor%gnk
                stop 'error'
            endif
        endif
        
        !! At this point, the following parts of the descriptor might be missing:
        ! lni (if dest is specified but not source or lni)
        ! lnk (if source is specified but not dest or lnk)
        ! lower bounds

        gni = descriptor%gni
        gnk = descriptor%gnk

        ! We can now take the opportunity to fix this up, completing the array specifiation.  To most evenly
        ! distribute points over processors, extents will be assigned through division, maximally spreading out
        ! any load imbalance.  For example, assigning 20 points over 6 processes would give a split of
        ! [3,3,4,3,3,4].

        if (sum(descriptor%lni) == 0) then        
            ! Define both lower bound and extent
            do ii=1,nproc_i
                descriptor%llbi(ii) = 1 + (gni*(ii-1))/nproc_i
                descriptor%lni(ii) = 1 + (gni*ii)/nproc_i - descriptor%llbi(ii)
            end do
        else
            ! Define just lower bound
            descriptor%llbi(1) = 1
            do ii=2,nproc_i
                descriptor%llbi(ii) = descriptor%llbi(ii-1) + descriptor%lni(ii-1)
            end do
        end if
        if (sum(descriptor%lnk) == 0) then
            ! Define both lower bound and extent
            do ii=1,nproc_i
                descriptor%llbk(ii) = 1 + (gnk*(ii-1))/nproc_i
                descriptor%lnk(ii) = 1 + (gnk*ii)/nproc_i - descriptor%llbk(ii)
            end do
        else
            ! Define just lower bound
            descriptor%llbk(1) = 1
            do ii=2,nproc_i
                descriptor%llbk(ii) = descriptor%llbk(ii-1) + descriptor%lnk(ii-1)
            end do
        end if

        !! With the logical size information complete, we can now start the MPI work.

        if (associated(descriptor%src_array))  then
            ! The source array was specified, and it is already associated in the descriptor.  We do not need
            ! to allocate memory, and we must use MPI_win_create for the MPI window.  Since we previously ensured
            ! that we don't mix supplied and allocated source arrays, this call will be correctly collective.

            ! As an implementation caution, note that the window creation implicitly assumes double-precision
            ! arrays (real(kind=real64)).

            call MPI_Win_create(descriptor%src_array, & ! Base address of the memory location
                                int(8*descriptor%lni(rank_i+1)*descriptor%lnj*descriptor%gnk,kind=MPI_ADDRESS_KIND), & ! Window size in bytes.
                                8, & ! Displacement unit, sizeof(double) allowing for index-based addressing over byte-based
                                MPI_INFO_NULL, & ! MPI_Info object -- possibly pass in no_locks
                                descriptor%communicator, & ! Communicator
                                descriptor%src_window,ierr) ! MPI shared memory window (ouptut)
        else
            ! The source array was not specified, so we should allocate it, via MPI_Win_allocate.
            ! This returns a type(c_ptr) at the base of the allocated memory, so we must also bind it
            ! to the src_array pointer in the descriptor.

            call MPI_Win_allocate(int(8*descriptor%lni(rank_i+1)*descriptor%lnj*descriptor%gnk,kind=MPI_ADDRESS_KIND), & ! Window size in bytes.
                                  8, & ! Displacement unit, sizeof(double) allowing for index-based addressing over byte-based
                                  MPI_INFO_NULL, & ! MPI_Info object -- possibly pass in no_locks
                                  descriptor%communicator, & ! Communicator
                                  allocated_window, & ! Base c_ptr (output)
                                  descriptor%src_window,ierr) ! MPI shared memory window (ouptut)

            ! Assign the allocated c_ptr to src_array in the descriptor, and specify extents
            call c_f_pointer(allocated_window, descriptor%src_array, [descriptor%lni(rank_i+1), descriptor%lnj, descriptor%gnk])
        end if

        ! The destination array is handled in essentially the same fashion
        if (associated(descriptor%dst_array))  then
            ! The destination array was alerady specified, so just create the window.
            ! Note that the proper shape of the destination array is (lnj, lnk, gni)

            call MPI_Win_create(descriptor%dst_array, & ! Base address of the memory location
                                int(8*descriptor%lnj*descriptor%lnk(rank_i+1)*descriptor%gni,kind=MPI_ADDRESS_KIND), & ! Window size in bytes.
                                8, & ! Displacement unit, sizeof(double) allowing for index-based addressing over byte-based
                                MPI_INFO_NULL, & ! MPI_Info object -- possibly pass in no_locks
                                descriptor%communicator, & ! Communicator
                                descriptor%dst_window,ierr) ! MPI shared memory window (ouptut)
        else
            ! The destination array is not specified, so allocate it.

            call MPI_Win_allocate(int(8*descriptor%lnj*descriptor%lnk(rank_i+1)*descriptor%gni,kind=MPI_ADDRESS_KIND), & ! Window size in bytes.
                                  8, & ! Displacement unit, sizeof(double) allowing for index-based addressing over byte-based
                                  MPI_INFO_NULL, & ! MPI_Info object -- possibly pass in no_locks
                                  descriptor%communicator, & ! Communicator
                                  allocated_window, & ! Base c_ptr (output)
                                  descriptor%dst_window,ierr) ! MPI shared memory window (ouptut)

            ! Assign the allocated c_ptr to src_array in the descriptor, and specify extents
            call c_f_pointer(allocated_window, descriptor%dst_array, [descriptor%lnj, descriptor%lnk(rank_i+1), descriptor%gni])
        end if

        ! If we're using passive-target synchronization, lock and unlock the windows on the self-process.
        ! This is not semantically required, but Vincent Magnoux reports at 
        ! https://gitlab.science.gc.ca/hpc/support/issues/742 that this is helpful for some odd reason in
        ! avoiding performance problems at the first nontrivial lock of the window
        if (SYNC_STYLE == SYNC_LOCK .or. SYNC_STYLE == SYNC_FLUSH) then
            call MPI_Barrier(descriptor%communicator)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank_i, MPI_MODE_NOCHECK, descriptor%src_window)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank_i, MPI_MODE_NOCHECK, descriptor%dst_window)
            call MPI_Win_unlock(rank_i, descriptor%src_window)
            call MPI_Win_unlock(rank_i, descriptor%dst_window)
            call MPI_Barrier(descriptor%communicator)
        end if

        ! If we're using SYNC_FLUSH synchronization, the window remains locked at all times and
        ! MPI_Win_flush(_all) is used to ensure ordering; the locks can be opened here.
        if (SYNC_STYLE == SYNC_FLUSH) then
            call MPI_Win_lock_all(MPI_MODE_NOCHECK, descriptor%src_window)
            call MPI_Win_lock_all(MPI_MODE_NOCHECK, descriptor%dst_window)
        end if

        ! Next, create the workspace arrays for the forward and reverse transposes.  These workspace arrays
        ! allow us to reconfigure the data-to-send to match the ordering on the other side of the MPI_Put call,
        ! turning the MPI RMA operation into a "dumb," contiguous remote memory copy.  This appears to be well-
        ! optimized by at least Intel MPI.

        ! For the "road not taken," an alternative approach would be to eschew the workspace arrays in favour
        ! of more complicated MPI datatypes.  By using MPI_type_create_hvector, it's possible to create one
        ! datatype per destination process such that the data elements are in the same logical order as the
        ! source while having a memory order that matches the destination.  

        ! See https://stackoverflow.com/questions/13093301/permuting-a-3-d-array-with-mpi-datatypes
        ! for an example of this approach in C.
        
        ! However, Intel MPI does not seem to optimize this case very much, and quick testing shows over two 
        ! orders of magnitude reduction in bandwidth for large transposes (100MByte/sec -> 19000 MByte/sec 
        ! for a dummy send that performs no in-memory transpose).

        ! The MPI type for the arrays are fixed, MPI_DOUBLE_PRECISION until and unless this code is extended
        ! to handle real(kind=real32) or noncontiguous source array slices.

        descriptor%mpitype = MPI_DOUBLE_PRECISION

        do ii = 1,nproc_i
            ! For the forward transpose, rank_i sends (lni(rank_i+1)), lnj, lnk(ii)) to the remote process.
            ! This data is in (j,k,i) order, matching the remote process's dst_array
            allocate(descriptor%dst_workspace(ii)%arr(descriptor%lnj, &
                                                      descriptor%lnk(ii),&
                                                      descriptor%lni(rank_i+1)))

            ! For the reverse transpose, rank_i sends (lni(ii), lnj, lnk(rank_i+1)) to the remote process.
            ! This data is in (i,j,k) order, matching the remote process's src_array
            allocate(descriptor%src_workspace(ii)%arr(descriptor%lni(ii), &
                                                      descriptor%lnj, &
                                                      descriptor%lnk(rank_i+1)))
        end do

        ! Finally, update the open_windows array to record the source and destination windows

        if (window_count == 0) then
            ! open_windows is not associated; allocate it with a default size
            allocate(open_windows(10))
            open_windows = MPI_WIN_NULL ! Initialize the array with a harmless value
            window_count = 2 ! Add two windows to the count, for src and dst_window
            descriptor%window_index = 1 ! Write the bookkeeping index of the first window to the descriptor
            ! Write the window objects to the open_windows bookkeeping array
            open_windows(1) = descriptor%src_window
            open_windows(2) = descriptor%dst_window
        else
            ! Check whether the open_windows array can accommidate this descriptor's windows
            if (ubound(open_windows,1) < window_count + 2) then
                block
                    ! Reallocate the open_windows array
                    type(MPI_Win), dimension(ubound(open_windows,1)) :: window_copy
                    type(MPI_Win), dimension(:), pointer :: new_open_windows
                    integer idx

                    do idx=1,ubound(open_windows,1)
                        window_copy(idx) = open_windows(idx)
                    end do
                    ! Allocate the new window bookeeping array; double the array size
                    allocate(new_open_windows((ubound(open_windows,1)*2)))
                    ! Initialize
                    new_open_windows = MPI_WIN_NULL
                    ! Deallocate the existing bookkeeping array and reassociate the pointer
                    deallocate(open_windows)
                    open_windows => new_open_windows
                end block
            end if
            ! Now, we're guaranteed to have enough space
            descriptor%window_index = window_count + 1
            open_windows(window_count + 1) = descriptor%src_window
            open_windows(window_count + 2) = descriptor%dst_window
            window_count = window_count + 2
        end if

        ! Mark the descriptor as initialized before returning
        descriptor%initialized = .true.
    end subroutine

    subroutine sync_open(window, group, err)
        !! Execute MPI RMA synchronization, opening an epoch where MPI_Puts to the provided
        !! window are valid.  The exact set of calls is determiend by the module-level SYNC_STYLE
        !! parameter.

        !! This subroutine makes MPI calls, but it does not provide OpenMP synchronization.  Calling
        !! routines must use directives (master/single/critical) as appropriate to prevent invalid
        !! or multiple calls.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(MPI_Win), intent(in) :: window !! MPI_Win to synchronize
        type(MPI_Group), intent(in) :: group !! MPI_Group of processes involved in RMA communication, for PSCW sync
        integer :: ierr !! MPI error parameter
        integer, intent(out), optional :: err !! Optional error parameter

        !! Synchronize the MPI window, opening an epoch where MPI_Puts to the window are valid.
        if (SYNC_STYLE == SYNC_LOCK) then
            !! Use passive-target synchronization; the source calls MPI_Win_lock (with a shared lock), after
            !! which MPI_Puts are valid.  The target process does not need to do anything; data just shows
            !! up.  After the source's call to MPI_Win_unlock, the MPI_Puts have completed on the remote target.
            !! This approach also requires an MPI_Barrier after synchronization to ensure that writes to this
            !! process have completed.
            call MPI_Win_lock_all(MPI_MODE_NOCHECK, & ! Assertion flag: no exclusive locks will be used
                                window,ierr) ! RMA window object
        elseif (SYNC_STYLE == SYNC_FENCE) then
            !! Use active-target syhcnronization with MPI_Win_fence.  This is collective over the full window.
            !! The first call to Win_fence opens an access/exposure epoch, and subsequent calls separate epochs.
            !! Communication is guaranteed to have completed at both the source and target after the second fence
            !! call.

            !! This first call can supply the hint MPI_MODE_NOPRECEDE, indicating that this call completes no
            !! RMA communication.
            call MPI_Win_fence(MPI_MODE_NOPRECEDE,window,ierr)
        elseif (SYNC_STYLE == SYNC_PSCW) then
            !! Use Post/Start/Complete/Wait synchronization.  MPI_Win_post announces an intention to communicate
            !! with a process group.  MPI_Win_start announces an intention to be communicated with by a process
            !! group, and it might block until all the corresponding _post calls have completed.

            !! To end the communication epoch, MPI_Win_complete ensures that all Puts issued by this process have
            !! completed on the local process, and MPI_Win_wait ensures that all writes targeting this process have
            !! completed.  After the _wait call, this process has finished the sends/receives for its part of the
            !! transpose, and no further collective synchronization should be required.

            !! As an implementation caution, there may be an unfixed synchronization bug with PSCW mode; large-nproc
            !! timing tests resulted in a randomly-reproducible model deadlock for unclear reasons.
            call MPI_Win_post(group, 0, window,ierr)
            call MPI_Win_start(group, 0, window,ierr)
        elseif (SYNC_STYLE == SYNC_FLUSH) then
            !! SYNC_FLUSH synchronization uses a long-lived lock on the MPI window, so no synchronization calls are
            !! required here.  We only need to set the ierr parameter for return.
            ierr = MPI_SUCCESS
        endif
        if (present(err)) then
            err = ierr
        end if
    end subroutine

    subroutine sync_close(window, comm, err)
        !! Execute MPI RMA synchronization, closing an epoch where MPI_Puts to the provided
        !! window are valid.  The exact set of calls is determiend by the module-level SYNC_STYLE
        !! parameter.

        !! This subroutine makes MPI calls, but it does not provide OpenMP synchronization.  Calling
        !! routines must use directives (master/single/critical) as appropriate to prevent invalid
        !! or multiple calls.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(MPI_Win), intent(in) :: window !! MPI_Win to synchronize
        type(MPI_Comm), intent(in) :: comm !! MPI Communicator responsible for the window, for local synchronization with SYNC_LOCK
        integer :: ierr !! MPI error parameter
        integer, intent(out), optional :: err !! Optional error parameter
    
        !! Synchronization at the end of the RMA communication
        if (SYNC_STYLE == SYNC_LOCK) then
            !! Passive-target synchronization calls unlock_all to ensure that this process's puts have completed,
            !! then a barrier to ensure that all proceses in the communicator have reached that point (thus, all
            !! writes destined to this process have completed)
            call MPI_Win_unlock_all(window,ierr)
            call MPI_barrier(comm,ierr)
        elseif (SYNC_STYLE == SYNC_FENCE) then
            !! Fence synchronization calls the fence, which ensures all writes have completed.  MPI_MODE_NOSUCCEED
            !! states that there will be no more RMA calls until the next fence, which is at the beginning of
            !! this function
            call MPI_Win_fence(MPI_MODE_NOSUCCEED,window,ierr)
        elseif (SYNC_STYLE == SYNC_PSCW) then
            !! MPI_Win_complete ensures that all RMA ops issued by this process have completed locally.  MPI_Win_wait
            !! ensures that all RMA ops with this process as a target have completed (with the corresponding _complete
            !! already issued.)  This implies a synchronization of all processes involved in RMA ops with this process.
            call MPI_Win_complete(window,ierr)
            call MPI_Win_wait(window,ierr)
        elseif (SYNC_STYLE == SYNC_FLUSH) then
            !! Use MPI_Win_flush_all to ensure that all communications started by this process have completed, both
            !! here and at the destination.  SYNC_FLUSH synchronization uses a long-lived lock on the window, so no
            !! unlocking is required.
            call MPI_Win_flush_all(window,ierr)
            call MPI_Barrier(comm,ierr)
        endif
        if (present(err)) then
            err = ierr
        end if
    end subroutine


    subroutine transpose_forward(descriptor)
        !! Using a transpose descriptor, execute the forward transpose.
        !! From the user perspective, this subroutine is blocking: it will not return until the transpose has
        !! completed on all processes, and descriptor%dst_array is fully available for local access.

        !! When OpenMP parallelism is used, THIS IS A COLLECTIVE CALL WITHIN A PARALLEL REGION.  This subroutine
        !! uses internal ("orphaned") work-sharing directives, and it will function incorrectly if called by some
        !! but not all threads in an active parallel region.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(transpose_descriptor), intent(inout) :: descriptor !! The descriptor object

        integer :: rank_i, numproc !! Rank and number of processes in communicator
        integer :: ii, jj, kk, pp, tosend ! Loop variable, number of elements to send
        integer(kind=MPI_ADDRESS_KIND) :: dst_disp ! Offset in the destination array for the target of MPI_Put
        integer :: ierr ! mpi error flag
        integer :: mpithread ! MPI threading support level
        logical :: mainthread 

        type(MPI_Status), dimension(descriptor%numproc) :: request_status

        rank_i = descriptor%rank
        numproc = descriptor%numproc
        call MPI_Query_thread(mpithread)
        call MPI_Is_thread_main(mainthread)

        !! Synchronize window, depending on threading level
        if (mpithread == MPI_THREAD_FUNNELED) then
            ! Only the master thread may issue MPI calls
            if (mainthread) then
                call sync_open(descriptor%dst_window, descriptor%group)
                descriptor%mpirequests = MPI_REQUEST_NULL
            end if
            !$omp barrier
        else
            ! Any thread can issue MPI calls, and synchronization must happen just once
            !$omp single
            call sync_open(descriptor%dst_window, descriptor%group)
            descriptor%mpirequests = MPI_REQUEST_NULL
            !$omp end single
        endif

        !! Loop over each process in the row.  Perform the local transpose of the data
        !! to send, and queue the MPI_Put operation.
        !!$omp do
        do pp=1,numproc
            ! Compute the number of elements to send to rank (pp)
            tosend = descriptor%lni(rank_i+1)*descriptor%lnj*descriptor%lnk(pp) 
            if (tosend == 0) then
                ! If that amount is zero, then there's no need to call MPI_Put
                cycle
            end if

            ! Locally transpose the data from the right portion of src_array to the
            ! correct element of dst_workspace.  OpenMP parallelism happens at this stage, since
            ! nthread is probably a better-balanced divisor of lni*lnk than it is of nproc
            if (pp /= rank_i+1) then
                !$omp do collapse(2)
                do ii=1,descriptor%lni(rank_i+1)
                    do kk=1,descriptor%lnk(pp)
                        ! The OMP simd directive here is required here for vectorization.  Without it, the Intel
                        ! compiler worries that there can be aliasing between the pointer variables, which is not
                        ! forbidden by the Fortran language standard.  The compiler can be reassured about these
                        ! variables' independence through an associate block, but in testing that caused the Intel
                        ! Fortran compiler to lose the information about array strides.
                        !$omp simd
                        do jj=1,descriptor%lnj
                            descriptor%dst_workspace(pp)%arr(jj,kk,ii) = descriptor%src_array(ii,jj,kk+descriptor%llbk(pp)-1)
                        end do
                    end do
                end do
            else
                !$omp do collapse(2)
                do ii=1,descriptor%lni(rank_i+1)
                    do kk=1,descriptor%lnk(pp)
                        !$omp simd
                        do jj=1,descriptor%lnj
                            descriptor%dst_array(jj,kk,ii+descriptor%llbi(rank_i+1)-1) = descriptor%src_array(ii,jj,kk+descriptor%llbk(pp)-1)
                        end do
                    end do
                end do
                cycle ! Do not queue an MPI Put, since the transpose is local.
            end if             

            ! Compute the destination displacement.  With (j,k,i) ordering, we want to write to the memory
            ! beginning at (1,1,llb_i(rank_i)), which has a one-dimensional offset given by:
            dst_disp = descriptor%lnj * descriptor%lnk(pp) * (descriptor%llbi(rank_i+1)-1)
            
            ! Enqueue the Put call.  We want to perform this only once per destination process.
            if (mpithread /= MPI_THREAD_MULTIPLE) then
                ! If we're not in a THREAD_MULTIPLE environment, then only one thread at a time may execute MPI
                ! calls.  Since this call is itself inside a loop, we would need to protect it with both $omp single
                ! and $omp critical directives, since the $omp single is independent per loop iteration.  However,
                ! we're already synchronizing at the end of the local transpose, so there's no problem in assigning
                ! the call to the master thread only.
                if (mainthread) then
                    call MPI_RPut(descriptor%dst_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                                tosend, & ! Number of elements to send
                                descriptor%mpitype, & ! Origin datatype
                                pp - 1, & ! Target rank (base 0 indexing)
                                dst_disp, & ! Offset (#elements) in destination array
                                tosend, & ! Target count (always 1)
                                descriptor%mpitype, & ! Destination datatype
                                descriptor%dst_window, & ! ierr)!, & ! Destination window
                                descriptor%mpirequests(pp), ierr) ! Request object
                end if
                ! No barrier is necessary.  Only master will execute the corresponding waitall call,
                ! and obviously it will have completed all of the required MPI_Rputs.
            else 
                ! Threads can execute MPI calls in parallel, but we want each individual remote put to happen
                ! only once.  The $omp single enforces the latter.
                !$omp single
                call MPI_RPut(descriptor%dst_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                            tosend, & ! Number of elements to send
                            descriptor%mpitype, & ! Origin datatype
                            pp - 1, & ! Target rank (base 0 indexing)
                            dst_disp, & ! Offset (#elements) in destination array
                            tosend, & ! Target count (always 1)
                            descriptor%mpitype, & ! Destination datatype
                            descriptor%dst_window, & ! ierr)!, & ! Destination window
                            descriptor%mpirequests(pp), ierr) ! Request object
                !$omp end single nowait
            end if

        end do

        ! Wait for completion of all put operations
        if (mpithread /= MPI_THREAD_MULTIPLE) then
            if (mainthread) then
                call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
                call sync_close(descriptor%dst_window, descriptor%communicator)
            end if
            !$omp barrier
        else
            !$omp barrier
            !$omp single
            call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
            call sync_close(descriptor%dst_window, descriptor%communicator)
            !$omp end single
        endif

    end subroutine

    subroutine transpose_forward_dummy(descriptor)
        !! Using a transpose descriptor, execute the forward transpose in a "dummy" fashion, without in-memory
        !! rearrangement.  This routine tests the implied performance penalty of complicated MPI types by
        !! using the same memory locations, windows, and data transfers while ignoring the "transpose" part
        !! of the process.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(transpose_descriptor), intent(inout) :: descriptor !! The descriptor object

        integer :: rank_i, numproc !! Rank and number of processes in communicator
        integer :: pp, tosend ! Loop variable, number of elements to send
        integer(kind=MPI_ADDRESS_KIND) :: dst_disp ! Offset in the destination array for the target of MPI_Put
        integer :: ierr ! Mpi error flag
        type(MPI_Status), dimension(descriptor%numproc) :: request_status
        integer :: mpithread
        logical :: mainthread 

        rank_i = descriptor%rank
        numproc = descriptor%numproc
        call MPI_Query_thread(mpithread)
        call MPI_Is_thread_main(mainthread)

        !! Synchronize window, depending on threading level
        if (mpithread == MPI_THREAD_FUNNELED) then
            ! Only the master thread may issue MPI calls, and this thread isn't it.
            if (mainthread) then
                call sync_open(descriptor%dst_window, descriptor%group)
                descriptor%mpirequests = MPI_REQUEST_NULL
            end if
            !$omp barrier
        else
            ! Any thread can issue MPI calls, and synchronization must happen just once
            !$omp single
            call sync_open(descriptor%dst_window, descriptor%group)
            descriptor%mpirequests = MPI_REQUEST_NULL
            !$omp end single
        endif

        do pp=1,numproc
            ! Compute the number of elements to send to rank (pp)
            tosend = descriptor%lni(rank_i+1)*descriptor%lnj*descriptor%lnk(pp) 
            if (tosend == 0) then
                ! If that amount is zero, then there's no need to call MPI_Put
                cycle
            end if

            !! The full forward_transpose does not execute an MPI_Put call for the "self" portion
            !! of the transpose.  However, if we outright skip it here we are excluding not just
            !! the local-transpose time, but also the entire data movement.  Thus, the "dummy" transpose
            !! here will include an MPI_Put for the pp=rank_i+1 transpose.

            dst_disp = descriptor%lnj * descriptor%lnk(pp) * (descriptor%llbi(rank_i+1)-1)
            
            if (mpithread /= MPI_THREAD_MULTIPLE) then
                if (mainthread) then
                    call MPI_RPut(descriptor%dst_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                                tosend, & ! Number of elements to send
                                descriptor%mpitype, & ! Origin datatype
                                pp - 1, & ! Target rank (base 0 indexing)
                                dst_disp, & ! Offset (#elements) in destination array
                                tosend, & ! Target count (always 1)
                                descriptor%mpitype, & ! Destination datatype
                                descriptor%dst_window, & ! ierr)!, & ! Destination window
                                descriptor%mpirequests(pp), ierr) ! Request object
                end if
            else 
                !$omp single
                call MPI_RPut(descriptor%dst_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                            tosend, & ! Number of elements to send
                            descriptor%mpitype, & ! Origin datatype
                            pp - 1, & ! Target rank (base 0 indexing)
                            dst_disp, & ! Offset (#elements) in destination array
                            tosend, & ! Target count (always 1)
                            descriptor%mpitype, & ! Destination datatype
                            descriptor%dst_window, & ! ierr)!, & ! Destination window
                            descriptor%mpirequests(pp), ierr) ! Request object
                !$omp end single nowait
            end if
        end do

        ! Wait for completion of all put operations
        if (mpithread /= MPI_THREAD_MULTIPLE) then
            if (mainthread) then
                call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
                call sync_close(descriptor%dst_window, descriptor%communicator)
            end if
            !$omp barrier
        else
            !$omp barrier
            !$omp single
            call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
            call sync_close(descriptor%dst_window, descriptor%communicator)
            !$omp end single
        endif
    end subroutine

    subroutine transpose_reverse(descriptor)
        !! Using a transpose descriptor, execute the reverse transpose.
        !! From the user perspective, this subroutine is blocking: it will not return until the transpose has
        !! completed on all processes, and descriptor%src_array is fully available for local access.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(transpose_descriptor), intent(inout) :: descriptor !! The descriptor object

        integer :: rank_i, numproc !! Rank and number of processes in communicator
        integer :: ii, jj, kk, pp, tosend ! Loop variable, number of elements to send
        integer(kind=MPI_ADDRESS_KIND) :: src_disp ! Offset in the source array for the target of MPI_Put
        integer :: ierr ! mpi error
        type(MPI_Status), dimension(descriptor%numproc) :: request_status
        integer :: mpithread
        logical :: mainthread

        rank_i = descriptor%rank
        numproc = descriptor%numproc
        call MPI_Query_thread(mpithread)
        call MPI_Is_thread_main(mainthread)

        !! For the reverse transpose, we take data initially in descriptor%dst_array and write it to
        !! the correct part of descriptor%src_array on the target process.  This uses descriptor%src_window.
        !! Synchronize window, depending on threading level
        if (mpithread == MPI_THREAD_FUNNELED) then
            if (mainthread) then
                call sync_open(descriptor%src_window,descriptor%group)
                descriptor%mpirequests = MPI_REQUEST_NULL
            end if
            !$omp barrier
        else
            !$omp single
            call sync_open(descriptor%src_window,descriptor%group)
            descriptor%mpirequests = MPI_REQUEST_NULL
            !$omp end single
        endif

        ! Loop over each process in the row, MPI_Putting the correct data
        do pp=1,numproc
            ! Compute the number of elements to send to rank (pp)
            tosend = descriptor%lni(pp)*descriptor%lnj*descriptor%lnk(rank_i+1) 
            if (tosend == 0) then
                ! If that amount is zero, then there's no need to call MPI_Put
                cycle
            end if

            
            if (pp /= rank_i+1) then
                ! Locally transpose the data from the right portion of dst_array to the
                ! correct element of src_workspace
                !$omp do collapse(2)
                do kk=1,descriptor%lnk(rank_i+1)
                    do jj=1,descriptor%lnj
                        !$omp simd
                        do ii=1,descriptor%lni(pp)
                            descriptor%src_workspace(pp)%arr(ii,jj,kk) = descriptor%dst_array(jj,kk,ii+descriptor%llbi(pp)-1)
                        end do
                    end do
                end do
            else
                ! Locally transpose the data directly to src_array, avoiding the MPI Put
                !$omp do collapse(2)
                do kk=1,descriptor%lnk(rank_i+1)
                    do jj=1,descriptor%lnj
                        !$omp simd
                        do ii=1,descriptor%lni(pp)
                            descriptor%src_array(ii,jj,kk+descriptor%llbk(rank_i+1)-1) = descriptor%dst_array(jj,kk,ii+descriptor%llbi(pp)-1)
                        end do
                    end do
                end do
                cycle ! No MPI call necessary
            end if

            ! Compute the destination displacement.  With (i,j,k) ordering, we want to write to the memory
            ! beginning at (1,1,llb_k(rank_i)), which has a one-dimensional offset given by:
            src_disp = descriptor%lnj * descriptor%lni(pp) * (descriptor%llbk(rank_i+1)-1)
            
            ! Call MPI_RPut, with proper threading guards
            if (mpithread /= MPI_THREAD_MULTIPLE) then
                if (mainthread) then
                    call MPI_RPut(descriptor%src_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                                tosend, & ! Number of elements to send
                                descriptor%mpitype, & ! Origin datatype
                                pp - 1, & ! Target rank (base 0 indexing)
                                src_disp, & ! Offset (#elements) in destination array
                                tosend, & ! Target count (always 1)
                                descriptor%mpitype, & ! Destination datatype
                                descriptor%src_window, & ! Destination window
                                descriptor%mpirequests(pp), ierr) ! Request slot
                end if
            else
                !$omp single
                call MPI_RPut(descriptor%src_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                            tosend, & ! Number of elements to send
                            descriptor%mpitype, & ! Origin datatype
                            pp - 1, & ! Target rank (base 0 indexing)
                            src_disp, & ! Offset (#elements) in destination array
                            tosend, & ! Target count (always 1)
                            descriptor%mpitype, & ! Destination datatype
                            descriptor%src_window, & ! Destination window
                            descriptor%mpirequests(pp), ierr) ! Request slot
                !$omp end single nowait
            endif
        end do

        ! Wait for completion of all put operations
        if (mpithread /= MPI_THREAD_MULTIPLE) then
            if (mainthread) then
                call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
                call sync_close(descriptor%src_window, descriptor%communicator)
            end if
            !$omp barrier
        else
            !$omp barrier
            !$omp single
            call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
            call sync_close(descriptor%src_window, descriptor%communicator)
            !$omp end single
        endif
        
    end subroutine

    subroutine transpose_reverse_dummy(descriptor)
        !! As transpose_reverse, but performing no in-memory transpose.  This tests the performance limit
        !! of _just_ the MPI calls.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(transpose_descriptor), intent(inout) :: descriptor !! The descriptor object

        integer :: rank_i, numproc !! Rank and number of processes in communicator
        integer :: pp, tosend ! Loop variable, number of elements to send
        integer(kind=MPI_ADDRESS_KIND) :: src_disp ! Offset in the source array for the target of MPI_Put
        integer :: ierr ! mpi error
        type(MPI_Status), dimension(descriptor%numproc) :: request_status
        integer :: mpithread
        logical :: mainthread 

        rank_i = descriptor%rank
        numproc = descriptor%numproc
        call MPI_Query_thread(mpithread)
        call MPI_Is_thread_main(mainthread)

        if (mpithread == MPI_THREAD_FUNNELED) then
            ! Only the master thread may issue MPI calls, and this thread isn't it.
            if (mainthread) then
                call sync_open(descriptor%src_window,descriptor%group)
                descriptor%mpirequests = MPI_REQUEST_NULL
            end if
            !$omp barrier
        else
            ! Any thread can issue MPI calls, and synchronization must happen just once
            !$omp single
            call sync_open(descriptor%src_window,descriptor%group)
            descriptor%mpirequests = MPI_REQUEST_NULL
            !$omp end single
        endif

        do pp=1,numproc
            tosend = descriptor%lni(pp)*descriptor%lnj*descriptor%lnk(rank_i+1) 
            if (tosend == 0) then
                ! If that amount is zero, then there's no need to call MPI_Put
                cycle
            end if
            src_disp = descriptor%lnj * descriptor%lni(pp) * (descriptor%llbk(rank_i+1)-1)            
            if (mpithread /= MPI_THREAD_MULTIPLE) then
                if (mainthread) then
                    call MPI_RPut(descriptor%src_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                                tosend, & ! Number of elements to send
                                descriptor%mpitype, & ! Origin datatype
                                pp - 1, & ! Target rank (base 0 indexing)
                                src_disp, & ! Offset (#elements) in destination array
                                tosend, & ! Target count (always 1)
                                descriptor%mpitype, & ! Destination datatype
                                descriptor%src_window, & ! Destination window
                                descriptor%mpirequests(pp), ierr) ! Request slot
                end if
            else
                !$omp single
                call MPI_RPut(descriptor%src_workspace(pp)%arr, & ! Address, in dst_workspace, of data to send
                            tosend, & ! Number of elements to send
                            descriptor%mpitype, & ! Origin datatype
                            pp - 1, & ! Target rank (base 0 indexing)
                            src_disp, & ! Offset (#elements) in destination array
                            tosend, & ! Target count (always 1)
                            descriptor%mpitype, & ! Destination datatype
                            descriptor%src_window, & ! Destination window
                            descriptor%mpirequests(pp), ierr) ! Request slot
                !$omp end single nowait
            endif
        end do

        if (mpithread /= MPI_THREAD_MULTIPLE) then
            if (mainthread) then
                call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
                call sync_close(descriptor%src_window, descriptor%communicator)
            end if
            !$omp barrier
        else
            !$omp barrier
            !$omp single
            call MPI_Waitall(numproc,descriptor%mpirequests,request_status)
            call sync_close(descriptor%src_window, descriptor%communicator)
            !$omp end single
        endif
    end subroutine

    subroutine transpose_destroy(descriptor)
        !! Deallocates the MPI objects associated with a transpose descriptor, allowing a successful program exit.
        !! In particular, Intel MPI does not appear to like it if MPI_Window objects persist through MPI_Finalize,
        !! so calling MPI_Win_free avoids the last-minute program crash.

        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        type(transpose_descriptor), intent(inout) :: descriptor
        integer :: ierr ! mpi error flag
        integer i
        integer :: mpithread !! Supported MPI threading level
        logical :: mainthread

        call MPI_Query_thread(mpithread)
        call MPI_Is_thread_main(mainthread)

        if (mpithread == MPI_THREAD_FUNNELED .and. .not. mainthread) then
            ! This thread is not allowed to execute MPI calls.  Block until the master
            ! thread deconstructs this descriptor.
            block
                logical :: exit_loop 
                exit_loop = .false.
                do while (.not. exit_loop)
                    !$omp critical(transpose_mod)
                    exit_loop = .not. descriptor%initialized
                    !$omp end critical(transpose_mod)
                end do
            end block
            return
        end if

        !$omp critical(transpose_mod)
        if (descriptor%initialized) then
            ! Free MPI one-sided windows.  This will also free the memory associated with an allocated window;
            ! it will not free an array that was passed in to transpose_create

            if (SYNC_STYLE == SYNC_FLUSH) then
                ! Unlock the windows
                call MPI_Win_unlock_all(descriptor%src_window)
                call MPI_Win_unlock_all(descriptor%dst_window)
            end if
            call MPI_Win_free(descriptor%src_window,ierr)
            call MPI_Win_free(descriptor%dst_window,ierr)

            ! Write tombstone values to the open window array
            open_windows(descriptor%window_index) = MPI_WIN_NULL
            open_windows(descriptor%window_index + 1) = MPI_WIN_NULL

            deallocate(descriptor%mpirequests)

            do i=1,descriptor%numproc
                deallocate(descriptor%src_workspace(i)%arr)
                deallocate(descriptor%dst_workspace(i)%arr)
            end do
            deallocate(descriptor%src_workspace)
            deallocate(descriptor%dst_workspace)
            deallocate(descriptor%lni, descriptor%lnk)
            deallocate(descriptor%llbi, descriptor%llbk)

            ! Finally, free the row communicator
            call MPI_Comm_free(descriptor%communicator,ierr)

            descriptor%initialized = .false.
        end if
        !$omp end critical(transpose_mod)
    end subroutine

    ! MPI Callback routine to free created MPI windows.  This function is attached to the MPI_COMM_SELF
    ! communicator, which causes it to be called when the self communicator is deleted as part of
    ! MPI_Finalize.  This in turn affords us the chance to free the created MPI_Win objects, so we can
    ! finally exit without causing errors in the Intel MPI library.
    subroutine transpose_clean_windows(comm, comm_keyval, attribute_val, extra_state, ierror) bind(C)
        ! Implement the MPI_Comm_copy_attr_function interface to create a callback function for the comm_free
        ! operation
        use iso_c_binding
        use, intrinsic :: iso_fortran_env
        use mpi_f08
        implicit none

        ! Because of the singleton nature of work being done -- freeing MPI Windows at program exit --
        ! then we don't need to check the input parameters.
        type(MPI_Comm) :: comm
        integer :: comm_keyval, ierror
        integer(kind=MPI_ADDRESS_KIND) :: attribute_val, extra_state
        integer :: idx, ierr

        ! Free all currently open MPI windows
        do idx = 1, window_count
            if (open_windows(idx) /= MPI_WIN_NULL) then
                if (sync_style == SYNC_FLUSH) then
                    ! With SYNC_FLUSH, the windows are kept in a locked state, and it must
                    ! be unlocked before the window can be freed
                    call MPI_Win_unlock_all(open_windows(idx),ierr)
                end if
                call MPI_Win_free(open_windows(idx),ierr)
            end if
        end do
        ierror = MPI_SUCCESS
    end subroutine

end module
