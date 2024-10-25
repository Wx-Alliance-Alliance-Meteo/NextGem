      real(kind=REAL64), intent(in) :: odx3rd
      real(kind=REAL64), parameter :: l1_d3rd = -1.d0/24.d0 , l2_d3rd = 9.d0/8.d0
      real(kind=REAL64), parameter :: l3_d3rd = -9.d0/8.d0 , l4_d3rd =  1.d0/24.d0
      real(kind=REAL64) :: val_deriv_3rd
      val_deriv_3rd = -odx3rd* (v1_3rd*l1_d3rd + v2_3rd*l2_d3rd + v3_3rd*l3_d3rd + v4_3rd*l4_d3rd)

