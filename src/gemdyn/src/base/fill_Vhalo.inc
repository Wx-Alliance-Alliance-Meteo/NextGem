      integer :: i,j,nl,n
      nl=1-Minx
      if (present(Nlvls)) nl=Nlvls
      if (nl>1-Minx) call gem_stop
      do n=1, nl
      do j=1-G_haloy+1,l_nj+G_haloy-1
         do i=1-G_halox+1,l_ni+G_halox-1
            F_dst(i,j, 1-n) = F_cte*F_dst(i,j,n)
            F_dst(i,j,G_nk+n) = F_cte*F_dst(i,j,G_nk-n+1  )
         end do
      end do
      end do