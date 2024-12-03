      do n= 0,F_k0,-1
         F_dst(:,:, n) = F_cte*F_dst(:,:,-(n-1))
      end do
      do n= G_nk+1,F_kn
         F_dst(:,:,n) = F_cte*F_dst(:,:,2*G_nk-n+1)
      end do
