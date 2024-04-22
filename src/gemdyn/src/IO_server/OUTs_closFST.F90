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
      subroutine OUTs_closFST (F_signal_L)
      use iso_c_binding
      use OUTs
      use clib_itf_mod
      implicit none

      logical, intent(IN) :: F_signal_L
      
#include <rmnlib_basics.hf>

      character(len=2048) :: filen, link
      integer :: i,err
!     
!--------------------------------------------------------------------
!
      if (OUTs_1o1_L) print*, 'Closing FST files unit', &
                             (out(i)%fnt,i=1,nout_files)
      do i=1,nout_files
         err = fstfrm(out(i)%fnt)
         err = fclos(out(i)%fnt)
         out(i)%fnt = 0
         out(i)%name= ''
         out(i)%type= ''
      end do

      if (OUTs_1o1_L.and.F_signal_L) then
         filen= trim(Out_dirname_S)//'/../../'//'output_ready_MASTER'
         link = trim(Out_dirname_S)//'/../../'//'output_ready'
         err = clib_symlink ( trim(filen), trim(link) )
      endif
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_closFST
