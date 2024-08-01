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

!The arguments:
!
!	1.  F_stk:  stack arrays, i.e the src and dst in adz_main.F90. The values
!               used for the interpolation and where the result is stored
!   2.  F_nptr: how many stack arrays, e.g for the thermo variables this would
!			    normally be 3 (w,t,z), but for continuity this should be 1
!   3.  F_xyz:  points to interpolate at; the result from the SL advection
!   4.  F_geom: ???
!	5.  F_num : Adz_num_q = how many q points there are?
!	6.  F_i0:   start i index of q, i.e Adz_i0
!	7.  F_in:   end i index of q, i.e Adz_in
!	8.  F_j0:   start j index of q, i.e Adz_j0
!	9.  F_jn:   end j index of q, i.e Adz_j 
!	10. F_k0:   start index of q, i.e 1
!	11. F_kn:   end index of q, i.e. l_nk+1
!	12. F_Quint_l: flag for order 5 interpolation
!


module SL_interp_mod
  use ISO_C_BINDING
  use adz_mem
  use adz_options
  use HORgrid_options
  implicit none
  public

contains
      subroutine SL_interp ( F_stk,F_nptr,F_xyz,F_geom,F_num,&
                             F_i0,F_in,F_j0,F_jn,F_k0,F_kn,F_Quint_L)
      implicit none

      logical, intent(IN), optional :: F_Quint_L
      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0,F_kn
      real, dimension(*), intent(in ) :: F_xyz
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      logical quint_L
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nk,nij,n,n1,n2,np,slc,dimext
      integer :: slice, HLT_np, HLT_start, HLT_end
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
      real, dimension(:,:,:,:), pointer :: extended
      real, dimension(:,:), pointer :: xchg
!
!---------------------------------------------------------------------
!
      ni = F_in-F_i0+1
      nj = F_jn-F_j0+1
      nk = F_kn-F_k0+1
      nij= ni*nj
      if (nij<=0) return
      dimext= Adz_nij*nk !total dimension; all d.o.f 
      quint_L = .false.
      if (present(F_Quint_L)) quint_L= F_Quint_L
      stkpntr(1)= c_loc(F_stk(1)%src) !grab the source array into stkpntr?
    !  call C_F_POINTER ( stkpntr(1), extended, [Adz_lmaxx-Adz_lminx+1,Adz_lmaxy-Adz_lminy+1,nk,F_nptr] )
      call C_F_POINTER ( stkpntr(1), extended, [ni,nj,nk,F_nptr] ) !save stkptr into extended array/pointer?
      
      !if there were more stack arrays, then stkpntr would grab the source array
      do n=1,F_nptr
         stkpntr(n)= c_loc(F_stk(n)%src(1,1,1))
      end do

      slice = ni*nj !one vertical slize 

      do slc= 1, F_num/slice + min(1,mod(F_num,slice))
         n1 = (slc-1)*slice + 1
         n2= min(n1+slice-1,F_num)
         np= n2-n1+1
         k1= n1/nij + 1
         k2= n2/nij + min(1,mod(n2,nij))
         ij= 0
         

         if (quint_L) then

            call SL_quint ( wrkc,lin,mi,ma,extended,&
              F_xyz((n1-1)*3+1),np, Adz_lminx,Adz_lmaxx,Adz_lminy,&
              Adz_lmaxy,dimext,1,F_k0,F_kn,'m' )

           else
              !call tricublin_zyx1_m_n ( wrkc,stkpntr,F_xyz((n1-1)*3+1),&
              !                          F_geom ,np, F_nptr )

            call SL_cubic ( wrkc,lin,mi,ma,extended,&
              F_xyz((n1-1)*3+1),np, Adz_lminx,Adz_lmaxx,Adz_lminy,&
              Adz_lmaxy,dimext,1,F_k0,F_kn,'m' )

         endif

           !save the interpolated array into the destination
           do n=1,F_nptr
              do k=k1,k2
                 kk= (k-1)*nj
                 j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
                 j2=min(nj,n2/ni  -kk) + F_j0 - 1
                 do j=j1,j2
                    do i=F_i0,F_in
                       ij= ij+1
                       F_stk(n)%dst(i,j,k)= wrkc(ij)
                    end do
                 end do
              end do
           end do

      end do
!---------------------------------------------------------------------
!
      return
      end subroutine SL_interp
      
end module SL_interp_mod
