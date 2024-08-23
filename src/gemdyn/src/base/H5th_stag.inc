      real(kind=REAL64), parameter :: intw1 = 3.d0/256.d0,  &
                                      intw2 = -25.d0/256.d0,&
                                      intw3 = 75.d0/128.d0
      real(kind=REAL64) :: val_stag
      val_stag = (v1+v6)*intw1 + (v2+v5)*intw2 + (v3+v4)*intw3