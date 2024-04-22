      subroutine OUTs_wrtref ( F_stag_S )
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use mpi_f08
      use IOs
      use OUTs
      use vGrid_Descriptors
      implicit none

      character(len=2), intent(IN) :: F_stag_S

#include <rmnlib_basics.hf>
      
      character(len=1) :: familly_uencode_S
      integer :: i0,in,j0,jn,vesion_uencode,niyy,sindx
      integer :: k,err,dimx,dimy,nis,njs,c1,c2
      integer, dimension(size(Level_allpres)) :: ip1s
      real(kind=REAL64), dimension(size(Level_allpres)) :: zero
      real :: wk
      real, dimension(:), allocatable :: yy
      real, dimension(:), pointer :: posx,posy
      type(vgrid_descriptor) :: vgd
!     
!--------------------------------------------------------------------
!
      i0 = max( 1   , Out_i0)
      in = min( G_ni, Out_in)
      j0 = max( 1   , Out_j0)
      jn = min( G_nj, Out_jn)
      nis = in - i0 + 1  ;  njs = jn - j0 + 1
      dimx= G_ni+2*G_halox ; dimy= G_nj+2*G_haloy
      
      if (F_stag_S(1:1) =='M') then
         posx => geomh_longs
         posy => geomh_latgs
         Out_ig3  = 1
      end if
      if (F_stag_S(1:1) =='U') then
         posx => geomh_longu
         posy => geomh_latgs
         Out_ig3  = 2
         in = min( G_ni-1, Out_in)
      end if
      if (F_stag_S(1:1) =='V') then
         posx => geomh_longs
         posy => geomh_latgv
         Out_ig3  = 3
         jn = min( G_nj-1, Out_jn)
      end if
      if (F_stag_S(1:1) =='F') then
         posx => geomh_longu
         posy => geomh_latgv
         Out_ig3  = 4
         in = min( G_ni-1, Out_in)
         jn = min( G_nj-1, Out_jn)
      end if
      Out_ig4 = 0

      call OUTs_igs ( Out_ig1, Out_ig2, posx(1), posy(1), G_ni,G_nj,&
                      Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4           ,&
                      i0,in,1, j0,jn,1 )

      if (OUTs_1o1_L) then
         if ( (Grd_yinyang_L) .and. (.not.Out_reduc_l) ) then
            vesion_uencode    = 1
            familly_uencode_S = 'F'

            niyy=5+2*(10+nis+njs)
            allocate (yy(niyy))

            yy(1 ) = iachar(familly_uencode_S)
            yy(2 ) = vesion_uencode
            yy(3 ) = 2          ! 2 grids (Yin & Yang)
            yy(4 ) = 1          ! the 2 grids have same resolution
            yy(5 ) = 1          ! the 2 grids have same area extension
!YIN
            sindx  = 6
            yy(sindx  ) = nis
            yy(sindx+1) = njs
            yy(sindx+2) = posx(i0)
            yy(sindx+3) = posx(i0+nis-1)
            yy(sindx+4) = posy(j0)
            yy(sindx+5) = posy(j0+njs-1)
            yy(sindx+6) = Grd_xlat1
            yy(sindx+7) = Grd_xlon1
            yy(sindx+8) = Grd_xlat2
            yy(sindx+9) = Grd_xlon2
            yy(sindx+10    :sindx+9+nis    )= &
            posx(i0:i0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(j0:j0+njs-1)
!YAN
            sindx  = sindx+10+nis+njs
            yy(sindx  ) = nis
            yy(sindx+1) = njs
            yy(sindx+2) = posx(i0)
            yy(sindx+3) = posx(i0+nis-1)
            yy(sindx+4) = posy(j0)
            yy(sindx+5) = posy(j0+njs-1)
            yy(sindx+6) = Grd_xlat1Y
            yy(sindx+7) = Grd_xlon1Y
            yy(sindx+8) = Grd_xlat2Y
            yy(sindx+9) = Grd_xlon2Y
            yy(sindx+10    :sindx+9+nis    )= &
            posx(i0:i0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(j0:j0+njs-1)

            err= fstecr(yy,yy, -32, Out_unf,Out_dateo,0,0,niyy,1,1  ,&
                        Out_ig1,Out_ig2,Out_ig3,'X','^>',Out_etik_S ,&
                        familly_uencode_S,vesion_uencode,0,0,0      ,&
                        5, .true.)
            deallocate (yy, STAT = err)
                         
         else
         err=fstecr(posx(i0),wk,-32,Out_unf,Out_dateo   ,&
                    0,0, nis,1,1, Out_ig1,Out_ig2,Out_ig3,'X', '>>'    ,&
                    Out_etik_S,'E',&
                    Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, 5, .true.)
         err=fstecr(posy(j0),wk,-32,Out_unf,Out_dateo   ,&
                    0,0, 1,njs,1,Out_ig1,Out_ig2,Out_ig3,'X', '^^'    ,&
                    Out_etik_S,'E',&
                    Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, 5, .true.)
         endif
                    
         if ( F_stag_S(2:2) == 'S' ) return
!if (( F_stag_S(2:2) == 'M' ) .or. ( F_stag_S(2:2) == 'T' )) then
         if ( F_stag_S(2:2) == 'M' ) then
            err = vgd_put  (gem_vgd,'IP_1 - record ip1',Out_ig1)
            err = vgd_put  (gem_vgd,'IP_2 - record ip2',Out_ig2)
            err = vgd_write(gem_vgd,unit=Out_unf,format='fst')
         else
            if ( F_stag_S(2:2) == 'P' ) then
               zero = 0.d0 ; c1=2 ; c2=1
               do k=1,size(Level_allpres)
                  call convip(ip1s(k),Level_allpres(k),c1,c2,'',.false.)
               end do
               err = vgd_new(vgd,            &
               kind     = 2,              &
               version  = 1,              &
               nk       = size(Level_allpres),     &
               ip1      = Out_ig1,        &
               ip2      = Out_ig2,        &
               a_m_8    = dble(Level_allpres*100.),&
               b_m_8    = zero,           &
               ip1_m    = ip1s)
               err = vgd_write(vgd,unit=Out_unf,format='fst')
               err = vgd_free(vgd)
            endif
         endif
      endif
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_wrtref
