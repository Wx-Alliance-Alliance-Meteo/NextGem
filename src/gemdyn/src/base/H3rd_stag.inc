      real(kind=REAL64), parameter :: int3rdw1 = -1.d0/16.d0
      real(kind=REAL64), parameter :: int3rdw2 =  9.d0/16.d0
      real(kind=REAL64) :: val_stag_3rd
      val_stag_3rd = (v1_3rd+v4_3rd)*int3rdw1 + (v2_3rd+v3_3rd)*int3rdw2
