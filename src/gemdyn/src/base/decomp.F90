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
!----------------------------------LICENCE END ---------------------------------
      
!! Defines grid decomposition parameters, dividing a dimension of global size F_npts over
!! the X or Y axis of the processor grid.  The resulting grid distribution optionally
!! includes a halo region, specified by the F_halo parameter.

!! After calling decomp, this process should define a decomposed array as
!! arr(F_min:F_max) (supposing a single dimension).  arr(1:F_lni) contains the valid
!! interior data for the decomposed array, with F_halo points 
      logical function decomp ( F_npts, F_min, F_max, F_lni, F_npartiel, F_halo, F_start,  &
                                F_alongx_L, F_fill_L, F_npe, F_lowestsize, F_checkparti_L, &
                                F_relax)
      use, intrinsic :: iso_fortran_env
      use ptopo               
      implicit none

      logical, intent(in) :: F_alongx_L, &  !! 1 (true) if the process distribution is along X, otherwise along Y
                             F_fill_L  , &  !! 1 (true) if the allocated array should be padded and "optimized for current machine"
                             F_checkparti_L !! 1 (true) if the decomposition should be checked for validity, but not computed
      integer, intent(inout) :: F_min, F_max, & !! Minimum and maximum bounds on the allocation of the decomposed array
                                F_Npartiel  , & !! Extent of the largest chunk given to any process
                                F_start         !! Index in the global array corresponding to arr(1)
      integer, intent(in) :: F_npts, & !! Global extent of the dimension to be split
                             F_halo, & !! Size of the halo region along the split dimension
                             F_npe     !! Number of processes along the dimension to split
      integer, intent(inout):: F_lni(F_npe)  !! Number of points owned by each process along the decomposition
      integer, intent(in) :: F_lowestsize, & !! Minimum number of points per process for a valid decomposition
                             F_relax !! Relaxation parameter:
                                     !! 0 -- constant points per process, except for the last process (invalid if last process gets 0 points)
                                     !! 1 -- F_lni(1) or F_lni(1)-1 points per process, with larger arrays at the beginning.  Must have more than 0 points per process.
                                     !! 2 -- As 1, but allows 0 points per process
                                     !! 3 -- F_lni(1) points for N processes, <F_lni(1) points for 1 process, then 0 points for remainder
                                     !! 4 -- Even point distribution; might have 0-point processes distributed throughout

      integer, external :: block_part
      integer gmin, gmax, ierr
      integer, dimension(F_npe) :: depl
!
!-------------------------------------------------------------------
!
      decomp = .false. ; F_lni= 0

      if (F_checkparti_L ) then

         ierr= block_part (0, F_npe, 1, F_npts, F_min, F_max, &
                                  F_lni, depl, F_relax)
         if (ierr >= 0) decomp = (  minval(F_lni) >= F_lowestsize )
                                  
      else

         if (F_alongx_L) then   ! distribute along X axis
          ierr = block_part(Ptopo_mycol, Ptopo_npex, 1, F_npts, gmin, gmax, F_lni, depl,F_relax)
        else                    ! distribute along Y axis
          ierr = block_part(Ptopo_myrow, Ptopo_npey, 1, F_npts, gmin, gmax, F_lni, depl,F_relax)
        end if
        decomp = (ierr == 0)
        if (decomp) then
           F_start = gmin       ! "tile" covers from nx0 to nxl in global space
           F_min = 1 - F_halo   ! "tile" dimension along this axis should be minx:maxx
           F_npartiel = F_lni(1)
           F_lni(1) = gmax - gmin + 1
           F_max = F_npartiel + F_halo
           if(F_fill_L) F_max = F_max + 1 - mod(F_npartiel,2) ! this needs to be revised for cache based machin
           if (F_lowestsize > 0) decomp = (F_lni(1) >= F_lowestsize)
        endif
      end if

!      if (.not.decomp) write (output_unit,1001) F_npts, F_npe
! 1001 format(/' DECOMP: illegal partitionning ====> ',i7,' / ',i7)
!
!-------------------------------------------------------------------
!
      return
      end

