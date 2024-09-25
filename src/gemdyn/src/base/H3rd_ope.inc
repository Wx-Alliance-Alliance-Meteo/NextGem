contains
      real(kind=REAL64) function Hderiv(v1_3rd,v2_3rd,v3_3rd,v4_3rd,odx3rd)
      real, intent(in) :: v1_3rd,v2_3rd,v3_3rd,v4_3rd
      include 'H3rd_deriv.inc'
      Hderiv = val_deriv_3rd
      end function Hderiv

      real(kind=REAL64) function Hderiv8(v1_3rd,v2_3rd,v3_3rd,v4_3rd,odx3rd)
      real(kind=REAL64), intent(in) :: v1_3rd,v2_3rd,v3_3rd,v4_3rd
      include 'H3rd_deriv.inc'
      Hderiv8 = val_deriv_3rd
      end function Hderiv8

      real(kind=REAL64) function Hstag(v1_3rd,v2_3rd,v3_3rd,v4_3rd)
      real             , intent(in) :: v1_3rd,v2_3rd,v3_3rd,v4_3rd
      include 'H3rd_stag.inc'
      Hstag = val_stag_3rd
      end function Hstag

      real(kind=REAL64) function Hstag8(v1_3rd,v2_3rd,v3_3rd,v4_3rd)
      real(kind=REAL64), intent(in) :: v1_3rd,v2_3rd,v3_3rd,v4_3rd
      include 'H3rd_stag.inc'
      Hstag8 = val_stag_3rd
      end function Hstag8
