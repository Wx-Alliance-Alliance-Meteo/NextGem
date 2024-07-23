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

!**s/r out_stkecr

      subroutine out_stkecr ( fa,lminx,lmaxx,lminy,lmaxy, &
                              meta,nplans, g_id,g_if,g_jd,g_jf )
      use iso_c_binding
      use out_collector, only: block_collect_fullp, Bloc_me
      use HORgrid_options
      use out_options
      use glb_ld
      use out_mod
      use out_meta
      use out3
      use ptopo
      use step_options
      use omp_timing
      implicit none

      integer, intent(IN) :: lminx,lmaxx,lminy,lmaxy,nplans
      integer, intent(IN) :: g_id,g_if,g_jd,g_jf
      real   , intent(IN) :: fa(lminx:lmaxx,lminy:lmaxy,nplans)
      type (meta_fstecr), dimension(:), pointer, intent(IN) :: meta

#include <rmnlib_basics.hf>
      include 'mpif.h'

      integer  nz, err, nis, njs, k, kk, wk_njs, tag, stat
      integer, dimension (:)    , pointer     :: zlist
      real   , dimension (:,:,:), pointer     :: wk, wk_glb
      real   , dimension (:,:)  , pointer     :: vec1,vec2
!
!----------------------------------------------------------------------
!
      nis = g_if - g_id + 1
      njs = g_jf - g_jd + 1
      wk_njs = -1
      if ( (nis < 1) .or. (njs < 1) ) return

      nullify (wk_glb, zlist)
      call block_collect_fullp ( fa, l_minx,l_maxx,l_miny,l_maxy, &
                                 nplans, wk_glb, nz, zlist )

      if (Bloc_me == 0) then
         if (nz>0) then
            if ((Grd_yinyang_L) .and. (Ptopo_couleur == 0)) then
               wk_njs = njs*2
            else
               wk_njs = njs
            end if
            allocate (wk(nis,wk_njs,nz))

      !      allocate (wk_glb(G_ni,G_nj,nz)) ; wk_glb=0. !to short-circuit the collector
            wk(1:nis,1:njs,1:nz) = wk_glb(g_id:g_if,g_jd:g_jf,1:nz)
            deallocate (wk_glb)
         else
            allocate (wk(1,1,1))
         end if

         allocate (vec1(nis*njs,1),vec2(nis*njs,2))

         do k= nz, 1, -1

            if (zlist(k) > 0) then
               kk= zlist(k)

               if ( (Grd_yinyang_L) .and. (.not.Out_reduc_l) ) then
                  !Merge from Yang (couleur 1) to Yin (couleur 0)
                  !reshape must have equal elements in src and destination

                  tag=401
                  if (Ptopo_couleur == 0) then
                     vec2 = reshape(wk(:,:,k), (/nis*njs,2/))
                     call MPI_recv ( vec2(1,2), nis*njs, MPI_REAL, 1, &
                                     tag, COMM_gridpeers, stat, err )
                     wk(:,:,k) = reshape(vec2, (/nis, wk_njs/))
                     err = fstecr ( wk(:,:,k),wk,-meta(kk)%nbits,Out_unf   ,&
                                    Out_dateo,Out_deet,Out_npas,nis,2*njs,1,&
                                    meta(kk)%ip1,meta(kk)%ip2,meta(kk)%ip3 ,&
                                    Out_typvar_S,meta(kk)%nv,Out_etik_S,'U',&
                                    meta(kk)%ig1,meta(kk)%ig2,meta(kk)%ig3 ,&
                                    Out_ig4,meta(kk)%dtyp,Out_rewrit_L )
                  else
                     vec1 = reshape(wk(:,:,k), (/nis*njs,1/))
                     call MPI_send ( vec1     , nis*njs, MPI_REAL, 0, &
                                     tag, COMM_gridpeers,         err )
                  end if
               else
                  err = fstecr ( wk(:,:,k),wk,-meta(kk)%nbits,Out_unf,&
                              Out_dateo,Out_deet,out_npas,nis,njs,1  ,&
                              meta(kk)%ip1,meta(kk)%ip2,meta(kk)%ip3 ,&
                              Out_typvar_S,meta(kk)%nv,Out_etik_S,'Z',&
                              meta(kk)%ig1,meta(kk)%ig2,meta(kk)%ig3 ,&
                              Out_ig4,meta(kk)%dtyp,Out_rewrit_L )
               end if
            end if

         end do

         if (associated(wk))    deallocate (wk)
         if (associated(vec1))  deallocate (vec1)
         if (associated(vec2))  deallocate (vec2)
         if (associated(zlist)) deallocate (zlist)

      end if
!
!--------------------------------------------------------------------
!
      return
      end
