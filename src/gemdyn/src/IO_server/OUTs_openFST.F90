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
      subroutine OUTs_openFST ()
      use iso_c_binding
      use OUTs
      implicit none
      
#include <rmnlib_basics.hf>

      integer :: i,indx,err
!     
!--------------------------------------------------------------------
!
!      err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)

      out(:)%fnt = 0
      if (OUTs_1o1_L) print*, nout_files,'FST files for Out_npas= ',Out_npas

      do i=1,nout_files
         out(i)%name= trim(Out_dirname_S)//'/'//trim(Out_filenames_S(i))
         indx = index(Out_filenames_S(i),"_")
         if (indx>5) then
            out(i)%type= Out_filenames_S(i)(1:2)
         else
            out(i)%type= Out_filenames_S(i)(1:indx-1)
         endif
         err= fnom  ( out(i)%fnt, trim(out(i)%name), 'STD+RND', 0 )
         err= fstouv( out(i)%fnt, 'RND' )
         if (OUTs_1o1_L) print*, 'FST: ',out(i)%fnt,trim(out(i)%name)
      end do
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_openFST
