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

!**s/r out_dq - Compute and output divergence and vorticity

      subroutine out_dq (levset,set)
      use vertical_interpolation
      use gmm_vt1
      use gmm_pw
      use gem_options
      use out_options
      use glb_ld
      use out_meta
      use out3
      use levels
      use svro_mod
      use outd
      use ver
      use rmn_gmm
      use outgrid
      implicit none

      integer, intent(IN) :: levset,set

      logical write_diag_lev
      integer i,istat,kind,nko,pndd,pnqq,pnqr,pnxx,gridset
      integer :: HLT_start, HLT_end, local_np, dim
      integer, dimension(:), allocatable :: indo
      real   , dimension(:), allocatable :: rf
      real, dimension(:), pointer, contiguous :: vorticity
      real, dimension(:,:,:), pointer :: uu_pres,vv_pres,cible, &
                                         div,vor,qr
!
!----------------------------------------------------------------------
!
      pndd=0 ; pnqq=0 ; pnqr=0 ; write_diag_lev = .false.

      do i=1,Outd_var_max(set)
        if (Outd_var_S(i,set) == 'DD') pndd=i
        if (Outd_var_S(i,set) == 'QQ') pnqq=i
        if (Outd_var_S(i,set) == 'QR') pnqr=i
      enddo

      if (pndd+pnqq+pnqr == 0) return

      Out_stag_S(3:3)= ' '
      if (Level_typ_S(levset) == 'M') then

         Out_stag_S= 'MM '
         kind= Level_kind_ip1
         allocate (indo( min(Level_max(levset),Level_momentum) ))
         call out_slev ( Level(1,levset), Level_max(levset), &
                          Level_momentum , indo,nko,write_diag_lev )

!!$omp parallel private (HLT_start, HLT_end, local_np)
         call HLT_split (1, 2*(G_nk+3), local_np, HLT_start, HLT_end)
         call gem_xch_halo ( ut1(l_minx,l_miny,HLT_start),&
                   l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
!!$omp end parallel

         if (pndd > 0) then
            allocate ( div(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
            call cal_div ( div, ut1, vt1 , Outd_filtpass(pndd,set),&
                           Outd_filtcoef(pndd,set)                ,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk )
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'Mass_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            call out_fstecr( div, l_minx,l_maxx,l_miny,l_maxy        ,&
                              Ver_hyb%m,'DD  ',Outd_convmult(pndd,set),&
                              Outd_convadd(pndd,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pndd,set),.false.)
            deallocate (div)
         endif

         if ((pnqq > 0).or.(pnqr > 0)) then

            Out_stag_S= 'FM '
            dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*l_nk
            allocate (vorticity(2*dim))
            vor(l_minx:l_maxx,l_miny:l_maxy,1:G_nk)=>vorticity(    1:)
            qr (l_minx:l_maxx,l_miny:l_maxy,1:G_nk)=>vorticity(dim+1:)
            if(pnqq > 0)then
               pnxx=pnqq
            else
               pnxx=pnqr
            endif
            call cal_vor ( qr, vor, ut1, vt1 , Outd_filtpass(pnxx,set),&
                           Outd_filtcoef(pnxx,set),(pnqq > 0)        ,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk )
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'F_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            if (pnqq > 0) &
            call out_fstecr( vor, l_minx,l_maxx,l_miny,l_maxy        ,&
                              Ver_hyb%m,'QQ  ',Outd_convmult(pnqq,set),&
                              Outd_convadd(pnqq,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pnqq,set),.false.)
            if (pnqr > 0) &
            call out_fstecr( qr, l_minx,l_maxx,l_miny,l_maxy         ,&
                              Ver_hyb%m,'QR  ',Outd_convmult(pnqr,set),&
                              Outd_convadd(pnqr,set),kind,-1          ,&
                              G_nk, indo,nko, Outd_nbit(pnqr,set),.false.)
            deallocate (vorticity)

          endif

      else

         istat= gmm_get(gmmk_pw_log_pm_s, pw_log_pm)
         kind= 2
         nko = Level_max(levset)
         dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*nko
         allocate (vorticity(2*dim))
         allocate ( indo(nko), rf(nko), &
                     cible(l_minx:l_maxx,l_miny:l_maxy,nko) )
         uu_pres(l_minx:l_maxx,l_miny:l_maxy,1:nko)=>vorticity(    1:)
         vv_pres(l_minx:l_maxx,l_miny:l_maxy,1:nko)=>vorticity(dim+1:)

         do i = 1, nko
            indo(i)= i
            rf  (i)= Level(i,levset)
            cible(:,:,i)= log(rf(i) * 100.0)
         enddo

         call vertint2 ( uu_pres,cible,nko, ut1,pw_log_pm,G_nk      ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_niu,1,l_nj,&
                         inttype=Out3_vinterp_type_S )
         call vertint2 ( vv_pres,cible,nko, vt1,pw_log_pm,G_nk      ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_njv,&
                         inttype=Out3_vinterp_type_S )

!!$omp parallel private (HLT_start, HLT_end, local_np)
         call HLT_split (1, 2*nko, local_np, HLT_start, HLT_end)
         call gem_xch_halo ( uu_pres(l_minx,l_miny,HLT_start),&
                         l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
!!$omp end parallel

         if (pndd > 0) then
            allocate ( div(l_minx:l_maxx,l_miny:l_maxy,nko) )
            call cal_div ( div, uu_pres, vv_pres  ,&
                           Outd_filtpass(pndd,set),&
                           Outd_filtcoef(pndd,set),&
                           l_minx,l_maxx,l_miny,l_maxy, nko )
            Out_stag_S= 'MP '
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'Mass_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            call out_fstecr( div, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'DD  ',Outd_convmult(pndd,set),&
                              Outd_convadd(pndd,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pndd,set),.false.)
            deallocate (div)
         endif

         if ((pnqq > 0).or.(pnqr > 0)) then

            allocate ( vor(l_minx:l_maxx,l_miny:l_maxy,nko),&
                        qr(l_minx:l_maxx,l_miny:l_maxy,nko) )
            if(pnqq > 0)then
               pnxx=pnqq
            else
               pnxx=pnqr
            endif
            call cal_vor ( qr, vor, uu_pres, vv_pres          ,&
                           Outd_filtpass(pnxx,set)            ,&
                           Outd_filtcoef(pnxx,set),(pnqq > 0),&
                           l_minx,l_maxx,l_miny,l_maxy, nko )
            Out_stag_S= 'FP '
            if ( .not. OUTs_server_L) then
            gridset = Outd_grid(set)
            call out_href ( 'F_point', &
                    OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                    OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )
            endif
            if (pnqq > 0) &
            call out_fstecr( vor, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'QQ  ',Outd_convmult(pnqq,set),&
                              Outd_convadd(pnqq,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pnqq,set),.false.)
            if (pnqr > 0) &
            call out_fstecr ( qr, l_minx,l_maxx,l_miny,l_maxy, &
                              rf,'QR  ',Outd_convmult(pnqr,set),&
                              Outd_convadd(pnqr,set), kind,-1  ,&
                              nko, indo, nko, Outd_nbit(pnqr,set),.false.)
            deallocate (vor, qr)
          endif

          deallocate (rf,cible,vorticity)

      endif

      deallocate (indo)
!
!----------------------------------------------------------------------
!
      return
      end
