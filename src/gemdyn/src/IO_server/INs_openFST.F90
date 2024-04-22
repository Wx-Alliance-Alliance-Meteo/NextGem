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

!**s/r inp_openFST - Open+link all fst input files valid at F_datev
!                    and determine vertical structure (kind)
      subroutine INs_openFST ( F_datev, F_err )
      use, intrinsic :: iso_fortran_env
      use clib_itf_mod
      use INs_base
      implicit none

      character(len=*), intent(IN ) :: F_datev
      integer         , intent(OUT) :: F_err

      logical, external :: INs_shared_mem, INs_gz3d
      character(len=2048) :: fn,root,mesg
      logical :: shared_mem_L
      integer :: i,j,err,err_code,unf
!
!-----------------------------------------------------------------------
!
      if (INs_1o1_L) err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
      fst(:)%fnt= -1 ; fst(:)%type= '0' ; SRC_GZ_L= .false.

      unf= 0
      fn= trim(Path_input_S)//'/ANALYSIS'
      if (fnom (unf,trim(fn),'RND+OLD+R/O',0) == 0) then
         if (fstouv(unf ,'RND') >= 0) then
            fst(1)%fnt = unf
            fst(1)%type= 'A'
            if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i4)') trim(fn),unf
         end if
      endif
      unf= 0
      fn= trim(Path_input_S)//'/GEOPHY/Gem_geophy.fst'
      if (fnom (unf,trim(fn),'RND+OLD+R/O',0) == 0) then
         if (fstouv(unf ,'RND') >= 0) then
            fst(2)%fnt = unf
            fst(2)%type= 'G'
            if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i4)') trim(fn),unf
         end if
      endif
      unf= 0
      fn= trim(Path_input_S)//'/CLIMATO'
      ! Climato is a directory: don't know what to do
      unf= 0
      fn= trim(Path_input_S)//'/MODEL_INPUT'
      ! Don't know what to do with MODEL_INPUT

      F_err = 0
      Inp_handle= -1 ; Inp_nfiles= 0 ; i= 0

      if ( any(SRL(:)%src == 'R') ) then
      Inp_datev= F_datev
      call datp2f ( Inp_cmcdate, F_datev )

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

 55   read (unf,'(a)',end=33) fn
      i = i+1
      fn= trim(root)//'/'//trim(fn)
      if (fnom (Inp_list_unf(i),trim(fn),'RND+OLD+R/O',0) == 0) then
         if (fstouv(Inp_list_unf(i) ,'RND') < 0) F_err= -1
      else
         F_err = -1
      end if
      if (F_err==0) then
         goto 55
      else
         Inp_list_unf(i)= 0
      endif
  
 33   if ((Inp_nfiles == 0).or.(i /= Inp_nfiles)) F_err= -1
      if (unf > 0) err= fclos(unf)
      
      if (F_err == 0) then
         err= fstlnk ( Inp_list_unf, Inp_nfiles )
         Inp_handle = Inp_list_unf(1)
         fst(3)%fnt = Inp_handle
         fst(3)%type= 'R'
         if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i4)') trim(root),Inp_handle
      endif
      if (F_err<0) return
      endif
      
      do i=1,INs_nreq
         do j=1,size(fst%type)
            if (fst(j)%type /= '0') then
               if (SRL(i)%src == fst(j)%type) then
                  SRL(i)%unf= fst(j)%fnt
                  cycle
               endif
            endif
         end do
      end do

! Establishing memory on first occurence of
! TT in an 'A' or 'R' input file
      shared_mem_L= .false.
      do i=1,INs_nreq
         if (shared_mem_L) exit
         if ( (SRL(i)%src == 'A').or.(SRL(i)%src == 'R')) then
            shared_mem_L= INs_shared_mem (SRL(i)%unf, 'TT')
            if (shared_mem_L) unf=SRL(i)%unf
         endif
      end do
      
      SRC_GZ_L= .false.
      if (shared_mem_L) SRC_GZ_L= INs_gz3d () !Obtain GZ if possible
      if (.not.SRC_GZ_L) then
         GZcaract(:)= 0
         INs_nreq=0
         F_err= -1
      endif
!
!-----------------------------------------------------------------------
!         
      return
      end subroutine INS_openFST
