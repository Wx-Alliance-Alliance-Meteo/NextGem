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

!**s/r inp_open - Open+link all fst input files valid at F_datev
!                 and determine vertical structure (kind)

      subroutine inp_open ( F_datev )
      use dynkernel_options
      use vGrid_Descriptors
      use path
      use clib_itf_mod
      use inp_base
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=*), intent(IN) :: F_datev

      include 'mpif.h'
      character(len=2048) fn,root,mesg
      integer i,err,err_code,unf,n123(3)
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
!
!-----------------------------------------------------------------------
!
      if (Inp_src_GZ_L) return
      
      Inp_handle= -1 ; Inp_nfiles= 0 ; err_code= 0 ; i= 0

      Inp_datev= F_datev
      call datp2f ( Inp_cmcdate, F_datev )

      if (Inp_iope >= 0) then
         root=trim(Path_input_S)//'/MODEL_INREP/VALID_'//trim(F_datev)
         err= clib_fileexist (trim(root)//'/content')
         if (err < 0) root=trim(Path_input_S)//&
                             '/MODEL_ANALYSIS/VALID_'//trim(F_datev)
         fn = trim(root)//'/content'
         unf= 0
         if (fnom( unf,trim(fn),'SEQ+FMT+OLD',0 ) /= 0) unf= 0
         if (unf == 0) goto 33
         read (unf,*,end=33) Inp_nfiles
         if (Inp_nfiles == 0) goto 33
         allocate (Inp_list_unf(Inp_nfiles))
         Inp_list_unf= 0

 55      read (unf,'(a)',end=33) fn
         i = i+1
         fn= trim(root)//'/'//trim(fn)

         if (fnom (Inp_list_unf(i),trim(fn),'RND+OLD+R/O',0) == 0) then
            if (fstouv(Inp_list_unf(i) ,'RND') < 0) err_code= -1
         else
            err_code= -1
         end if
         if (err_code == 0) goto 55
 33      if ((Inp_nfiles == 0).or.(i /= Inp_nfiles)) err_code= -1
         if (unf > 0) err= fclos(unf)
         if (err_code == 0) then
            err= fstlnk ( Inp_list_unf, Inp_nfiles )
            Inp_handle = Inp_list_unf(1)
         endif
      endif
      
      call gem_error ( err_code, 'inp_open', &
                      'Problems opening input files' )

      Inp_src_GZ_L= inp_src_vert3d ('GEOPOTENTIAL',GZ3d)

      nullify (vtbl_8) ; n123= -1
      Inp_kind= GZ3d%kind ; Inp_version= -1 ; Inp_src_hauteur_L= .true.
      if ((Inp_src_GZ_L).and.(Inp_kind == 21)) return
      
      if (Inp_iope >= 0) then
            err= vgd_new ( Inp_vgd_src, unit=Inp_handle, &
                           format='fst', ip1=-1, ip2=-1 )
            if (err == 0) then
               err= vgd_get ( Inp_vgd_src, 'VTBL', vtbl_8, quiet=.true.)
               n123(1:3) = ubound(vtbl_8)
            end if
      endif

      call MPI_bcast (n123, 3, MPI_INTEGER, Inp_bcast, COMM_grid, err)

      if (n123(1) > 0) then
         if (Inp_iope /= 0) allocate(vtbl_8(n123(1),n123(2),n123(3)))
         call MPI_bcast ( vtbl_8,size(vtbl_8), &
             MPI_DOUBLE_PRECISION, Inp_bcast, COMM_grid, err )
         if (Inp_iope /= 0) err= vgd_new ( Inp_vgd_src, vtbl_8 )
         deallocate (vtbl_8)
         
         err = vgd_get ( Inp_vgd_src, key='KIND',value=Inp_kind    )
         err = vgd_get ( Inp_vgd_src, key='VERS',value=Inp_version )

         Inp_src_hauteur_L = .not. (Inp_kind == 2)

         if ( (Inp_kind == 21) .and. (.not.Inp_src_GZ_L)) &
             Inp_src_GZ_L= inp_src_vert3d (Inp_vgd_src, GZ3d)
         
         if ( (Inp_kind == 5) .or. &
             ((Inp_kind == 1).and.(Inp_version == 3)) ) then
            err= vgd_get ( Inp_vgd_src, key='PREF',value=Inp_pref_a_8 )
         end if
      end if
      
      if ((.not.Inp_src_GZ_L).and.Inp_src_hauteur_L) err_code= -1
      call gem_error ( err_code, 'inp_open', &
         'Unable to determine input file vertical structure')
!
!-----------------------------------------------------------------------
!         
      return
      end
