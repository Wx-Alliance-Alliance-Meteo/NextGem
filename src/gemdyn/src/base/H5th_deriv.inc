      real(kind=REAL64), intent(in) :: h
      real(kind=REAL64), parameter :: c1= 1.d0/64.d0, c2= 0.1d0, c3= 1.d0/6.d0
      real(kind=REAL64) :: val_deriv
      val_deriv = (h*c1) * ( 3.d0  * (v6 - v1) *c2 &
                           + 25.d0 * (v2 - v5) *c3 &
                           + 75.d0 * (v4 - v3) )
