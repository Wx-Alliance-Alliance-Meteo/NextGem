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

!**s/r INs_closeFST - Close all fst input files already opened

      subroutine INs_closeFST (F_err)
      use IOs
      use INs
      implicit none

      integer, intent(OUT) :: F_err

#include <rmnlib_basics.hf>

      integer i, j, err_code(size(fst%fnt)),err
!
!--------------------------------------------------------------------
!
      F_err= 0 ; err_code=0
      do i=1,size(fst%fnt)
         if (fst(i)%fnt>0) then
            if (fst(i)%type/='R') then
               if (fstfrm (fst(i)%fnt) == 0) then
                  err_code(i)= fclos (fst(i)%fnt)
                  if ((err_code(i)==0).and.(lun_out>0)) &
                  write (6,'("Fortran unit:",i4,"is closed")') fst(i)%fnt

               else
                  err_code(i)= -1
               end if
            else
               err_code(i)= 0
               do j= 1, Inp_nfiles
                  if (fstfrm (Inp_list_unf(j)) == 0) then
                     err= fclos (Inp_list_unf(j))
                  else
                     err= -1
                  endif
                  if (err<0) then
                     err_code(i)= -1
                  else
                     if (lun_out>0) write (6,'("Fortran unit:",i4," is closed")') Inp_list_unf(j)
                  endif
               end do
            endif
         endif
      end do
      F_err= minval(err_code)
      fst(:)%fnt= -1 ; fst(:)%type= '0'

      if (associated(Inp_list_unf)) then
         deallocate (Inp_list_unf)
         nullify (Inp_list_unf)
      endif
      if (associated(GZ)) then
         call MPI_Win_free (GZ_win  ,err)
         nullify(GZ) ; deallocate (GZIP1)
      endif
      if (associated(ND)) then
         call MPI_Win_free (NEST_win,err)
         nullify(ND) ; deallocate (DIP1)
      endif
      INs_nia= -1 ; INs_nja= -1 ; INs_nka= -1
      if (INs_vgd_L) err= vgd_free(Inp_vgd_src)
!
!--------------------------------------------------------------------
!
      return
      end subroutine INs_closeFST
