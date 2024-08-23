contains
      real(kind=REAL64) function Hderiv(v1,v2,v3,v4,v5,v6,h)
      real, intent(in) :: v1,v2,v3,v4,v5,v6
      include 'H5th_deriv.inc'
      Hderiv = val_deriv
      end function Hderiv

      real(kind=REAL64) function Hderiv8(v1,v2,v3,v4,v5,v6,h)
      real(kind=REAL64), intent(in) :: v1,v2,v3,v4,v5,v6
      include 'H5th_deriv.inc'
      Hderiv8 = val_deriv
      end function Hderiv8

      real(kind=REAL64) function Hstag(v1,v2,v3,v4,v5,v6)
      real             , intent(in) :: v1,v2,v3,v4,v5,v6
      include 'H5th_stag.inc'
      Hstag = val_stag
      end function Hstag

      real(kind=REAL64) function Hstag8(v1,v2,v3,v4,v5,v6)
      real(kind=REAL64), intent(in) :: v1,v2,v3,v4,v5,v6
      include 'H5th_stag.inc'
      Hstag8 = val_stag
      end function Hstag8
