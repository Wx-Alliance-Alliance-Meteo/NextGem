      integer :: i,j
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            F_dst(i,j, 0) = F_cte*F_dst(i,j,1)
            F_dst(i,j,-1) = F_cte*F_dst(i,j,2)
            F_dst(i,j,-2) = F_cte*F_dst(i,j,3)
            F_dst(i,j,G_nk+1) = F_cte*F_dst(i,j,G_nk  )
            F_dst(i,j,G_nk+2) = F_cte*F_dst(i,j,G_nk-1)
            F_dst(i,j,G_nk+3) = F_cte*F_dst(i,j,G_nk-2)
         end do
      end do