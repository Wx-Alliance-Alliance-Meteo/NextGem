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

!**s/r dynstep - Advance the dynamics by one timestep

      subroutine dynstep ()
      use ISO_C_BINDING
      use step_options
      use theo_options
      use gmm_vt1
      use gmm_vt2
      use gmm_pw
      use gmm_contiguous
      use mem_tstp
      use adz_mem
      implicit none

      integer i,j,k,icn,np
      integer :: HLT_start, HLT_end, local_np
!
!     ---------------------------------------------------------------
!
      call pw_switch ()

!!$omp do collapse(2)
      do k=1, l_nk
         do j=1, l_nj
            do i= 1, l_ni
               Adz_uu_ext(i,j,k) = pw_uu_moins(i,j,k)
               Adz_vv_ext(i,j,k) = pw_vv_moins(i,j,k)
               Adz_ww_ext(i,j,k) =        zdt1(i,j,k)
               Adz_uu_dep_ext(i,j,k) = pw_uu_plus(i,j,k)
               Adz_vv_dep_ext(i,j,k) = pw_vv_plus(i,j,k)
               Adz_ww_dep_ext(i,j,k) =       zdt2(i,j,k)
            end do
         end do
      end do
!!$omp enddo

      call HLT_split (1, 1, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( st1, l_minx,l_maxx,l_miny,l_maxy,local_np,-1 )
      np= 4*G_nk+1
      call HLT_split (1, np, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( ut1(l_minx,l_miny,HLT_start),&
                l_minx,l_maxx,l_miny,l_maxy,local_np,-1)
      call HLT_split (1, 6*l_nk, local_np, HLT_start, HLT_end)
      call gem_xch_halo ( Adz_uu_ext(Adz_lminx,Adz_lminy,HLT_start),&
                Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, local_np,-1) 
                
      if (Dynamics_sw_L) then
                   
         if (Step_kount.eq.1) then
            call SW_tstpdyn(0.75d0*Cstv_dt_8) !Shallow-Water (1st substep = EULER)
            call t02t1 ()
            call SW_tstpdyn(0.5d0 *Cstv_dt_8) !Shallow-Water (2nd substep = BDF2 )
         else
            call SW_tstpdyn(Cstv_dt_8) !Shallow-Water (all other timesteps= BDF2 )
         endif

      else

         if (Step_kount.eq.1) then
            call tstpdyn(0.75d0*Cstv_dt_8) !Fully-Compressible Euler (1st substep = EULER)
            
            call t02t1 ()
            call t02t2 ()
!!$omp do collapse(2)
            do k=1, l_nk
               do j=1, l_nj
                  do i= 1, l_ni
                     Adz_uu_ext(i,j,k) = ut1(i,j,k)
                     Adz_vv_ext(i,j,k) = vt1(i,j,k)
                     Adz_ww_ext(i,j,k) = zdt1(i,j,k)
                     Adz_uu_dep_ext(i,j,k) = ut2(i,j,k)
                     Adz_vv_dep_ext(i,j,k) = vt2(i,j,k)
                     Adz_ww_dep_ext(i,j,k) = zdt2(i,j,k)
                  end do
               end do
            end do
!!$omp enddo

            call HLT_split (1, 6*l_nk, local_np, HLT_start, HLT_end)
            call gem_xch_halo ( Adz_uu_ext(Adz_lminx,Adz_lminy,HLT_start),&
                      Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, local_np,-1)

            call tstpdyn(0.5d0*Cstv_dt_8 ) !Fully-Compressible Euler (2nd substep = BDF2 )

!     special for tracer advection only
            call adz_traject(Cstv_dt_8)

         else  
            if (Step_kount.eq.2) then
               dynt2= var_init
            else
               call t02t2 ()
            endif
            call tstpdyn(Cstv_dt_8) !Fully-Compressible Euler (all other timesteps= BDF2)

         endif

      endif

      if (Ctrl_theoc_L .and. .not.Grd_yinyang_L) call theo_bndry ()

      call adz_tracers_interp ()

      call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,&
                      pw_log_pt,pw_pm_plus_8,pw_p0_plus_8       ,&
                      l_minx,l_maxx,l_miny,l_maxy,l_nk,0 )
      
      if ( Schm_psadj > 0 ) call psadj ( Step_kount )

      call adz_tracers_massfixing ()
 
      call t02t1 ()
      call pw_update_GW ()      
      call pw_update_UV ()
      call pw_update_T  ()

      if (Grd_yinyang_L) &
      call yyg_xchng_vec_q2q (pw_uu_plus,pw_vv_plus,&
                    l_minx,l_maxx,l_miny,l_maxy,G_nk)

      call HOR_bndry ()
!
!     ---------------------------------------------------------------
!
      return
      end subroutine dynstep
