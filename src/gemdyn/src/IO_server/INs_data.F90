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

!**s/r INs_data  - Reads FST input files

      subroutine INs_data ()
      use INs_base
      implicit none

      character(len=9) :: mesg
      integer :: i,deb1,deb2,err
!
!-----------------------------------------------------------------------
!
      if (INs_1o1_L) err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
      SRL(:)%vname(2) = ''
      
      INs_nplans=0
      do i=1,INs_nreq
         if (INs_1o1_L) print*, 'TREATING: ',SRL(i)%vname(1), SRL(i)%stag, SRL(i)%unf
         mesg='FOUND'
         Inp_handle= SRL(i)%unf
         if (trim(SRL(i)%vname(1)) == 'UV') then
            err= INs_hwnd (SRL(i)%vname(2),SRL(i)%nk,SRL(i)%deb)
            if (err<0) mesg='NOT FOUND'
         else
            deb1= INs_nplans+1
            deb2= INs_nplans*INs_hord+1
            SRL(i)%deb= deb1
            err= INs_read_mt (SRL(i)%vname(1),SRL(i)%stag,&
                 SRL(i)%vname(2),ND(deb2:),len(trim(SRL(i)%stag)),&
                 DIP1(deb1),SRL(i)%nk)
            INs_nplans= INs_nplans+max(0,SRL(i)%nk)
            if (err<0) mesg='NOT FOUND'
         endif
      end do
      if (INs_1o1_L) err = fstopc ('MSGLVL','WARNIN',RMN_OPT_SET)
!
!-----------------------------------------------------------------------
!
      return
      end subroutine INs_data
