! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The Reaction Rates File
!
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
!
! ADOM2.eqn - Rx 36 and 43 are merged with Rx 35 and 42, respectively
! { 1.} NO2 + hv =  NO + O; RCONST(1) = JCOR  ! JNO2 in sec-1; photo_type1
! { 2.} O   + O2  + M  = O3; RCONST(2) = (3.00e-28/(TEMP**2.3)); type 6
! { 3.} O + NO2 = NO; RCONST(3) = (ARR2(6.5E-12,120.0)); type1
! { 4.} O + NO2   = NO3; RCONST(4) = (TYPE5(0.60,8.10e-27,-2.0,2.20e-11,0.0)); type5
! { 5.} NO + O3 = NO2; RCONST(5) = (ARR2(1.8E-12,-1370.0)); type1
! { 6.} NO2 + O3 =  NO3; RCONST(6) = (ARR2(1.2E-13,-2450.0)); type1
! { 7.} NO3 + NO = 2 NO2; RCONST(7) = (ARR2(1.7E-11,150.0)); type1
! { 8.} 2 NO  + O2  =  2 NO2; RCONST(8) = (ARR2(3.30E-39,529.0)); type1
! { 9.} NO3 + NO2   =  N2O5; RCONST(9) = (TYPE5(0.60,9.86e-20,-4.3,2.60e-11,-0.5)); type5
! {10.} N2O5   =  NO3 + NO2; RCONST(10) = (RCONST(9)*ARR2(9.09e26,-11200.)); type 6
! {11.} N2O5 + H2O = 2 HNO3; RCONST(11) = 1e-21; type0
! {12.} NO3 + NO2 = NO + NO2; RCONST(12) = (ARR2(2.5E-14,-1229.0)); type1
! {13.} NO3 + hv =  NO; RCONST(13) = 2.08 *RCONST(1); photo_type4
! {14.} NO3 + hv =  NO2 + O; RCONST(14) = 18.9 *RCONST(1); photo_type4
! {15.} O3 + hv = O; RCONST(15) = 5.2E-2*RCONST(1); photo_type4
! {16.} O3 + hv = O1D; RCONST(16) = 3.8E-3*RCONST(1); photo_type1
! {17.} O1D + H2O = 2OH; RCONST(17) = 2.2e-10; type0
! {18.} O1D + M = O; RCONST(18) = 2.9e-11; type0
! {19.} OH + NO   =  HONO; RCONST(19) = (TYPE5(0.60,1.93e-24,-2.6,2.60e-10,-0.5)); type5
! {20.} HONO + hv =  OH + NO; RCONST(20) = 0.181*RCONST(1); photo_type4
! {21.} NO2 + H2O = HONO + HNO3 - NO2; RCONST(21) = 1e-24; type0
! {22.} OH + NO2   =  HNO3; RCONST(22) = (TYPE5(0.60,2.20e-22,-3.2,4.00e-8,-1.3)); type5
! {23.} HNO3 + hv  =  NO2 + OH; RCONST(23) = 6.4e-5*RCONST(1); photo_type4
! {24.} OH + HNO3   =  NO3; RCONST(24) = (ARR2(9.4E-15,778.0)); type1
! {25.} OH + CO  = HO2; RCONST(25) = 2.4e-13; type 6
! {26.} O3 + OH = HO2; RCONST(26) = (ARR2(1.6E-12,-942.0)); type1
! {27.} HO2 + NO = OH + NO2; RCONST(27) = (ARR2(3.7E-12,240.0)); type1
! {28.} NO2 + HO2 = HNO4; RCONST(28) = (TYPE5(0.60,1.52e-23,-3.2,1.38e-8,-1.4)); type5
! {29.} HNO4 = NO2 + HO2; RCONST(29) = (RCONST(28)*ARR2(4.76e26,-10940.0)); type 6
! {30.} HNO4 + hv = NO2 + HO2; RCONST(30) = 9.0e-4*RCONST(1); photo_type4
! {31.} HNO4 + OH = NO2; RCONST(31) = (ARR2(1.30e-12,380.0)); type1
! {32.} O3 + HO2 = OH; RCONST(32) = (ARR2(1.1E-14,-502.0)); type1
! {33.} 2 HO2 = H2O2; RCONST(33) = (ARR2(2.2E-13,619.0)); type1
! {34.} 2 HO2 + M = H2O2; RCONST(34) = (ARR2(1.9E-33,982.0)); type1
! {35.+ 36.} 2 HO2 + H2O = H2O2 + {36.} H2O2 + DUMMY; RCONST(35) = (ARR2(3.1E-34,2818.0)+ARR2(2.700e-54,3137.)*CFACTOR*1.0D6); ! type1
! {37.} H2O2 + hv = 2 OH; RCONST(36) = 8.4E-4*RCONST(1); photo_type4
! {38.} OH + H2O2 = HO2; RCONST(37) = (ARR2(3.3E-12,-200.0)); type1
! {39.} NO3 + HO2 = HNO3; RCONST(38) = (ARR2(2.27e-13,771.0)); type1
! {40.} NO3 + HO2 + M = HNO3; RCONST(39) = (ARR2(1.9e-33,982.0)); type1
! {41.} NO3 + HO2 + H2O = HNO3; RCONST(40) = (ARR2(3.1e-34,2818.0)+ARR2(2.700e-54,3137.)*CFACTOR*1.0D6); type1
! {42.+ 43.} NO3 + HO2 + H2O = HNO3 + DUMMY + {43}  SO4 + HO2; RCONST(41) = (TYPE5(0.60,4.48e-23,-3.3,1.50e-12,0.0)); type 6
! {44.} RO2 + NO  =  NO; RCONST(42) = (ARR2(4.20e-12,180.0)); type1
! {45.} RO2 + HO2 =  HO2; RCONST(43) = (ARR2(1.75e-13,1000.0)); type1
! {46.} 2 RO2 = PROD; RCONST(44) = 1e-15; type0
! {47.} RO2 + MCO3 =  MCO3; RCONST(45) = 3e-12; type0
! {48.} ROOH + hv =  HO2 + OH; RCONST(46) = 8.2e-4*RCONST(1); photo_type4
! {49.} HCHO + hv {+ 2 O2} = 2 HO2 + CO; RCONST(47) = 3.3E-3*RCONST(1); photo_type4
! {50.} HCHO +  hv = CO; RCONST(48) = 4.8E-3*RCONST(1); photo_type4
! {51.} HCHO + OH  =  HO2 + CO; RCONST(49) = (ARR2(1.60e-11,-110.0)); type1
! {52.} HCHO + NO3  = HNO3 + HO2 + CO; RCONST(50) = (ARR2(2.80e-12,-2518.0)); type1
! {53.} HCHO + HO2 =  RO2R + RO2; RCONST(51) = (1.1E-13*(1.-20./(20.+ARR2(4.2E-18,180.0)*C(ind_NO)/CFACTOR*1.E6_dp*CFACTOR))); type 6
! {54.} ALD2 + OH = MCO3; RCONST(52) = (ARR2(5.60e-12,311.0)); type1
! {55.} ALD2 + hv {+ 2 O2} = HCHO + RO2 + RO2R  + CO + HO2; RCONST(53) = 6.6E-4*RCONST(1); photo_type4
! {56.} ALD2 + NO3  = MCO3 + HNO3; RCONST(54) = (ARR2(1.40e-12,-1867.0)); type1
! {57.} MCO3 + NO  = HCHO + RO2 + RO2R + NO2; RCONST(55) = (ARR2(4.20e-12,180.0)); type1
! {58.} MCO3 + NO2 = PAN; RCONST(56) = (TYPE5(0.19,6.29e-19,-4.1,4.92e-3,-3.6)); type5
! {59.} MCO3 + HO2 =  ROOH + HCHO; RCONST(57) = (ARR2(1.75e-13,1000.0)); type1
! {60.} 2 MCO3 = 2 HCHO + 2 HO2; RCONST(58) = 5.3e-12; type0
! {61.} PAN = MCO3 + NO2; RCONST(59) = (ARR2(2.00e+16,-13542.0)); type1
! {62.} MEK + hv = ALD2 + MCO3 + RO2R + RO2; RCONST(60) = 1.76e-4*RCONST(1); photo_type4
! {63.} MEK + OH = 0.5 ALD2 + 0.5 HCHO + 1.5 RO2R + 1.5 RO2 + MCO3; RCONST(61) = (ARR2(1.20e-11,-745.0)); type1
! {64.} MGLY + hv = MCO3 + CO + HO2; RCONST(62) = 1.67e-2*RCONST(1); photo_type4
! {65.} OH + MGLY = MCO3 + CO; RCONST(63) = 1.7e-11; type0
! {66.} NO3 + MGLY = CO + HNO3 + MCO3; RCONST(64) = (ARR2(3.00e-13,-1427.0)); type1
! {67.} CH4 + OH =  HCHO + RO2R + RO2; RCONST(65) = (ARR2(2.40e-12,-1710.0)); type1
! {68.} C2H6 + OH = ALD2 + RO2R + RO2; RCONST(66) = (ARR2(1.70e-11,-1232.0)); type1
! {69.} C3H8 + OH = 0.3 ALD2 + 0.5 MEK + RO2R + RO2; RCONST(67) = (TEMP*TEMP*ARR2(1.27e-17,14.0)); type 6
! {70.} ALKA + OH = 0.111 HCHO + 0.530 ALD2  + 0.640 MEK + 0.131 RO2N + 0.869 RO2R + 0.7 R2O2 + 1.7 RO2; RCONST(68) =
!     (ARR2(1.017e-11,-354.0)*0.517+ARR2(2.312e-11,-289.0)*(1.-0.517)); type 6
! {71.} RNO3 + OH = NO2 + 0.14 MEK + 1.52 ALD2 + 0.16 HCHO + 1.39 R2O2 + 1.39 RO2; RCONST(69) = (ARR2(2.19e-11,-709.0)); type1
! {72.} RO2N + NO = RNO3; RCONST(70) = (ARR2(4.20e-12,180.0)); type1
! {73.} RO2N + HO2 = ROOH + MEK; RCONST(71) = (ARR2(1.75e-13,1000.0)); type1
! {74.} RO2N + RO2 = RO2 + 0.5 HO2 + MEK; RCONST(72) = 1e-15; type0
! {75.} RO2N + MCO3 = HCHO + HO2 + MEK; RCONST(73) = 3e-12; type0
! {76.} R2O2 + NO = NO2; RCONST(74) = (ARR2(4.20e-12,180.0)); type1
! {77.} R2O2 + HO2 = ROOH; RCONST(75) = (ARR2(1.75e-13,1000.0)); type1
! {78.} R2O2 + RO2 = RO2; RCONST(76) = 1e-15; type0
! {79.} R2O2 + MCO3 = HCHO + HO2; RCONST(77) = 3e-12; type0
! {80.} RO2R + NO = NO2  + HO2; RCONST(78) = (ARR2(4.20e-12,180.0)); type1
! {81.} RO2R + HO2 = ROOH; RCONST(79) = (ARR2(1.75e-13,1000.0)); type1
! {82.} RO2R + RO2 = 0.5 HO2 + RO2; RCONST(80) = 1e-13; type0
! {83.} RO2R + MCO3 = HCHO + HO2; RCONST(81) = 3e-12; type0
! {84.} OH + ETHE = RO2 + 1.56 HCHO  + RO2R + 0.22 ALD2; RCONST(82) = (ARR2(2.15e-12,411.0)); type1
! {85.} O3 + ETHE = HCHO + 0.42 CO + 0.12 HO2 + 0.4 CRG1; RCONST(83) = (ARR2(1.20e-14,-2634.0)); type1
! {86.} O + ETHE = HCHO + HO2 + CO + RO2R + RO2; RCONST(84) = (ARR2(1.04e-11,-792.0)); type1
! {87.} NO3 + ETHE = 2 HCHO + NO2 + R2O2 + RO2; RCONST(85) = (ARR2(3.70e-12,-2925.0)); type1
! {88.} OH + ALKE = 0.667 HCHO + 1.334 ALD2 + RO2R + RO2; RCONST(86) =
!     (ARR2(5.323e-12,504.0)*0.667+ARR2(1.074e-11,549.0)*(1.-0.667)); type 6
! {89.} O3 + ALKE = 0.667 ALD2 + 0.427 HCHO + 0.187 CO + 0.183 HO2 + 0.177 RO2 + 0.177 RO2R+ 0.080 OH + 0.133 CRG1 + 0.133 CRG2;
!     RCONST(87) = (ARR2(1.323e-14,-2105.0)*0.667+ARR2(7.333e-15,-1137.0)*(1.-0.667)); type 6
! {90.} O + ALKE = 0.133 ALD2 + 0.267 HO2 + 0.400 RO2 + 0.267 CO + 0.267 HCHO + 0.400 RO2R + 0.333 MEK; RCONST(88) =
!     (ARR2(1.18e-11,-324.0)*0.667+ARR2(2.26e-11,10.0)*(1.-0.667)); type 6
! {91.} NO3 + ALKE =  NO2 + 0.667 HCHO + 1.334 ALD2 + R2O2 + RO2; RCONST(89) =
!     (ARR2(1.143e-11,-1935.0)*0.667+ARR2(3.23e-11,-975.0)*(1.-0.667)); type 6
! {92.} SO2 + CRG1 = SO4 + HCHO; RCONST(90) = 1e-13; type0
! {93.} SO2 + CRG2 = SO4 + ALD2; RCONST(91) = 1e-13; type0
! {94.} CRG1 + H2O = PROD; RCONST(92) = 2.3e-17; type0
! {95.} CRG2 + H2O = PROD; RCONST(93) = 2.3e-17; type0
! {96.} CRG1 + HCHO = PROD; RCONST(94) = 2.5e-14; type0
! {97.} CRG2 + HCHO = PROD; RCONST(95) = 2.5e-14; type0
! {98.} CRG1 + ALD2 = PROD; RCONST(96) = 2.5e-14; type0
! {99.} CRG2 + ALD2 = PROD; RCONST(97) = 2.5e-14; type0
! {100.} TOLU + OH = 0.16 CRES + 0.16 HO2 + 0.84 RO2R + 0.4 DIAL + 0.84 RO2 + 0.13 MGLY + 0.11 HCHO + 0.11 CO;
!        RCONST(98) = (ARR2(2.10e-12,322.0)); type1
! {101.}  AROM + OH = 0.17 CRES + 0.17 HO2 + 0.83 RO2R + 0.83 RO2 + 0.590 DIAL + 0.518 MGLY + 0.0597 HCHO + 0.0597 CO;
!        RCONST(99) = (ARR2(1.407e-11,116.0)*0.628+4.77e-11*(1.-0.628)); type 6
! {102.} DIAL + OH = MCO3; RCONST(100) = 3e-11; type0
! {103.} DIAL + hv = HO2 + CO + MCO3; RCONST(101) = 5.89e-2*RCONST(1); photo_type4
! {104.} OH + CRES = 0.2 MGLY + 0.14 RO2N + 0.85 RO2R + RO2 + 0.08 CRES; RCONST(102) = 4e-11; type0
! {105.} NO3 + CRES = HNO3 + BZO + 0.5 CRES; RCONST(103) = 2.2e-11; type0
! {106.} NO2 + BZO = RNO3; RCONST(104) = 1.5e-11; type0
! {107.} HO2 + BZO = PROD; RCONST(105) = (ARR2(1.75e-13,1000.0)); type1
! {108.} BZO = PROD; RCONST(106) = 0.001; type0
! {109.} OH + ISOP = HCHO + RO2 + 0.7 HO2 + 0.27 MGLY + 0.2 MCO3 + ETHE + 0.2 ALD2 + 0.1 RO2N + 0.9 R2O2 + 0.320 MVK + 0.230 MACR;
!        RCONST(107) = (ARR2(1.50e-11,500.0)); type1
! {110.} O3 + ISOP = HCHO + 0.4 ALD2 + 0.5 ETHE + 0.2 MGLY + 0.2 CRG2 + 0.4 HO2 + 0.1 OH + 0.160 MVK + 0.39 MACR;
!        RCONST(108) = (ARR2(7.00e-15,-1900.0)); type1
! {111.} O + ISOP =  0.6 HO2 + ALD2 + 0.5 RO2 + 0.5 R2O2 + ETHE; RCONST(109) = 1.8e-11; type0
! {112.} NO3 + ISOP =  NO2 + HCHO + ALD2 + R2O2 + RO2; RCONST(110) = 3.2e-13; type0
! {113.} OH + HO2 =  PROD; RCONST(111) = (ARR2(4.60e-11,230.0)); type1
! {114.} OH + ROOH =  0.5 RO2R + 0.5 OH; RCONST(112) = (ARR2(4.00e-12,180.0)); type1
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!! Modified for ADOM box model (convert to subroutine)
!! removed unused functions
! return gas-phase only reaction rate for input temperature
! modified to output 2D rate constants

!! Note by Kenjiro Toyota (2017-09-11)
!! Double precision bugs have been corrected manually for ARR2
!
!!if_on
subroutine mach_kpp_adom2_rates(rct, temp2d, p2d, cfactor2d, cno2d, gni, gnk)
  use mach_pkg_adom2_mod,   only: nreact, nprcf, nreact_kpp
!!if_off
  use chm_nml_mod,          only: nk_start
  use mach_adom2_rates_mod, only: mach_adom2_bcs_coeffs, xx, xc, y, yc, z, zc
  implicit none
!
!!if_on
  integer(kind=4), intent(in) :: gni, gnk
  real(kind=8),    intent(out) :: rct(gni * gnk, nreact)
  real(kind=4),    intent(in)  :: temp2d(gni, gnk) ! Deg. Kelvin
  real(kind=4),    intent(in)  :: cfactor2d(gni, gnk) ! conversion from ppm to molecules/cm3
  real(kind=4),    intent(in)  :: cno2d(gni, gnk)  ! NO conc in molec/cm3
  real(kind=4),    intent(in)  :: p2d(gni, gnk)    ! Pressure profile in Pa
!!if_off
!
  integer(kind=4) :: i, k, jj, nn, nr, npt
  real(kind=8)    :: cfactor, temp, cno
  real(kind=8)    :: bcs(gni * gnk, nprcf)

  rct = tiny(1.0d0)

  nn = (nk_start - 1) * gni
  do k = nk_start, gnk
  do i = 1, gni
     nn = nn + 1
     temp = dble(temp2d(i, k))
     cfactor = dble(cfactor2d(i, k))
     cno = dble(cno2d(i, k))

   ! RCT(nn,1) = ! Photolysis J reaction  1: no2 + hv ---> no
     RCT(nn,2) = (3.00d-28/(TEMP**2.3d0))
     RCT(nn,3) = (ARR2(6.5D-12,120.d0,TEMP))
     RCT(nn,4) = (TYPE5(0.60d0,8.10d-27,-2.0d0,2.20d-11,0.0d0,TEMP,CFACTOR))
     RCT(nn,5) = (ARR2(1.8D-12,-1370.d0,TEMP))
     RCT(nn,6) = (ARR2(1.2D-13,-2450.d0,TEMP))
     RCT(nn,7) = (ARR2(1.7D-11,150.d0,TEMP))
     RCT(nn,8) = (ARR2(3.30D-39,529.d0,TEMP))
     RCT(nn,9) = (TYPE5(0.60d0,9.86d-20,-4.3d0,2.60d-11,-0.5d0,TEMP,CFACTOR))
     RCT(nn,10) = (RCT(nn,9)*ARR2(9.09D+26,-11200.d0,TEMP))
   ! RCT(nn,11) = constant rate coefficient
     RCT(nn,12) = (ARR2(2.5D-14,-1229.d0,TEMP))
   ! RCT(nn,13) = ! Photolysis J reaction 13: NO3 + hv =  NO
   ! RCT(nn,14) = ! Photolysis J reaction 14: NO3 + hv =  NO2 + O
   ! RCT(nn,15) = ! Photolysis J reaction 15: O3 + hv = O
   ! RCT(nn,16) = ! Photolysis J reaction 16: O3 + hv = O1D
   ! RCT(nn,17) = constant rate coefficient
   ! RCT(nn,18) = constant rate coefficient
     RCT(nn,19) = (TYPE5(0.60d0,1.93d-24,-2.6d0,2.60d-10,-0.5d0,TEMP,CFACTOR))
   ! RCT(nn,20) = ! Photolysis J reaction 20: HONO + hv =  OH + NO
   ! RCT(nn,21) = constant rate coefficient
     RCT(nn,22) = (TYPE5(0.60d0,2.20d-22,-3.2d0,4.00d-8,-1.3d0,TEMP,CFACTOR))
   ! RCT(nn,23) = ! Photolysis J reaction 23: HNO3 + hv  =  NO2 + OH
     RCT(nn,24) = (ARR2(9.4D-15,778.d0,TEMP))
   ! RCT(nn,25) = constant rate coefficient
     RCT(nn,26) = (ARR2(1.6D-12,-942.d0,TEMP))
     RCT(nn,27) = (ARR2(3.7D-12,240.d0,TEMP))
     RCT(nn,28) = (TYPE5(0.60d0,1.52d-23,-3.2d0,1.38d-8,-1.4d0,TEMP,CFACTOR))
     RCT(nn,29) = (RCT(nn,28)*ARR2(4.76D+26,-10940.d0,TEMP))
   ! RCT(nn,30) = ! Photolysis J reaction 30: HNO4 + hv = NO2 + HO2
     RCT(nn,31) = (ARR2(1.30D-12,380.d0,TEMP))
     RCT(nn,32) = (ARR2(1.1D-14,-502.d0,TEMP))
     RCT(nn,33) = (ARR2(2.2D-13,619.d0,TEMP))
     RCT(nn,34) = (ARR2(1.9D-33,982.d0,TEMP))
     RCT(nn,35) = (ARR2(3.1D-34,2818.d0,TEMP)+ARR2(2.700D-54,3137.d0,TEMP)*CFACTOR*1.0D6)
   ! RCT(nn,36) = ! Photolysis J reaction 37: H2O2 + hv = 2 OH
     RCT(nn,37) = (ARR2(3.3D-12,-200.d0,TEMP))
     RCT(nn,38) = (ARR2(2.27D-13,771.d0,TEMP))
     RCT(nn,39) = (ARR2(1.9D-33,982.d0,TEMP))
     RCT(nn,40) = (ARR2(3.1D-34,2818.d0,TEMP)+ARR2(2.700D-54,3137.d0,TEMP)*CFACTOR*1.0D6)
     RCT(nn,41) = (TYPE5(0.60d0,4.48d-23,-3.3d0,1.50d-12,0.0d0,TEMP,CFACTOR))
     RCT(nn,42) = (ARR2(4.20D-12,180.d0,TEMP))
     RCT(nn,43) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,44) = constant rate coefficient
   ! RCT(nn,45) = constant rate coefficient
   ! RCT(nn,46) = ! Photolysis J reaction 46: ROOH + hv =  HO2 + OH
   ! RCT(nn,47) = ! Photolysis J reaction 47: HCHO + hv {+ 2 O2} = 2 HO2 + CO
   ! RCT(nn,48) = ! Photolysis J reaction 48: HCHO +  hv = CO
     RCT(nn,49) = (ARR2(1.60D-11,-110.d0,TEMP))
     RCT(nn,50) = (ARR2(2.80D-12,-2518.d0,TEMP))
     RCT(nn,51) = (1.1d-13*(1.0d0-20.0d0/(20.0d0+ARR2(4.2D-18,180.d0,TEMP)* CNO/CFACTOR*1.d6*CFACTOR)))
     RCT(nn,52) = (ARR2(5.60D-12,311.d0,TEMP))
   ! RCT(nn,53) = ! Photolysis J reaction 53: ALD2 + hv {+ 2 O2} = HCHO + RO2 + RO2R  + CO + HO2
     RCT(nn,54) = (ARR2(1.40D-12,-1867.d0,TEMP))
     RCT(nn,55) = (ARR2(4.20D-12,180.0d0,TEMP))
     RCT(nn,56) = (TYPE5(0.19d0,6.29d-19,-4.1d0,4.92d-3,-3.6d0,TEMP,CFACTOR))
     RCT(nn,57) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,58) = constant rate coefficient
     RCT(nn,59) = (ARR2(2.00D+16,-13542.d0,TEMP))
   ! RCT(nn,60) = ! Photolysis J reaction 60:  MEK + hv = ALD2 + MCO3 + RO2R + RO2
     RCT(nn,61) = (ARR2(1.20D-11,-745.d0,TEMP))
   ! RCT(nn,62) = ! Photolysis J reaction 62: MGLY + hv = MCO3 + CO + HO2
   ! RCT(nn,63) = constant rate coefficient
     RCT(nn,64) = (ARR2(3.00D-13,-1427.d0,TEMP))
     RCT(nn,65) = (ARR2(2.40D-12,-1710.d0,TEMP))
     RCT(nn,66) = (ARR2(1.70D-11,-1232.d0,TEMP))
     RCT(nn,67) = (TEMP*TEMP*ARR2(1.27D-17,14.d0,TEMP))
     RCT(nn,68) = (ARR2(1.017D-11,-354.d0,TEMP)*dble(xx)+ARR2(2.312D-11,-289.d0,TEMP)*dble(xc))
     RCT(nn,69) = (ARR2(2.19D-11,-709.d0,TEMP))
     RCT(nn,70) = (ARR2(4.20D-12,180.d0,TEMP))
     RCT(nn,71) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,72) = constant rate coefficient
   ! RCT(nn,73) = constant rate coefficient
     RCT(nn,74) = (ARR2(4.20D-12,180.d0,TEMP))
     RCT(nn,75) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,76) = constant rate coefficient
   ! RCT(nn,77) = constant rate coefficient
     RCT(nn,78) = (ARR2(4.20D-12,180.d0,TEMP))
     RCT(nn,79) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,80) = constant rate coefficient
   ! RCT(nn,81) = constant rate coefficient
     RCT(nn,82) = (ARR2(2.15D-12,411.d0,TEMP))
     RCT(nn,83) = (ARR2(1.20D-14,-2634.d0,TEMP))
     RCT(nn,84) = (ARR2(1.04D-11,-792.d0,TEMP))
     RCT(nn,85) = (ARR2(3.70D-12,-2925.d0,TEMP))
     RCT(nn,86) = (ARR2(5.323D-12,504.d0,TEMP)*dble(y)+ARR2(1.074D-11,549.d0,TEMP)*dble(yc))
     RCT(nn,87) = (ARR2(1.323D-14,-2105.d0,TEMP)*dble(y)+ARR2(7.333D-15,-1137.d0,TEMP)*dble(yc))
     RCT(nn,88) = (ARR2(1.18D-11,-324.d0,TEMP)*dble(y)+ARR2(2.26D-11,10.d0,TEMP)*dble(yc))
     RCT(nn,89) = (ARR2(1.143D-11,-1935.d0,TEMP)*dble(y)+ARR2(3.23D-11,-975.d0,TEMP)*dble(yc))
   ! RCT(nn,90) = constant rate coefficient
   ! RCT(nn,91) = constant rate coefficient
   ! RCT(nn,92) = constant rate coefficient
   ! RCT(nn,93) = constant rate coefficient
   ! RCT(nn,94) = constant rate coefficient
   ! RCT(nn,95) = constant rate coefficient
   ! RCT(nn,96) = constant rate coefficient
   ! RCT(nn,97) = constant rate coefficient
     RCT(nn,98) = (ARR2(2.10D-12,322.d0,TEMP))
     RCT(nn,99) = (ARR2(1.407D-11,116.d0,TEMP)*dble(z)+4.77d-11*dble(zc))
   ! RCT(nn,100) = constant rate coefficient
   ! RCT(nn,101) = ! Photolysis J reaction 101 DIAL + hv = HO2 + CO + MCO3
   ! RCT(nn,102) = constant rate coefficient
   ! RCT(nn,103) = constant rate coefficient
   ! RCT(nn,104) = constant rate coefficient
     RCT(nn,105) = (ARR2(1.75D-13,1000.d0,TEMP))
   ! RCT(nn,106) = constant rate coefficient
     RCT(nn,107) = (ARR2(1.50D-11,500.d0,TEMP))
     RCT(nn,108) = (ARR2(7.00D-15,-1900.d0,TEMP))
   ! RCT(nn,109) = constant rate coefficient
   ! RCT(nn,110) = constant rate coefficient
     RCT(nn,111) = (ARR2(4.60D-11,230.d0,TEMP))
     RCT(nn,112) = (ARR2(4.00D-12,180.d0,TEMP))

  end do
  end do

  RCT(:,11) = 1d-21
  RCT(:,17) = 2.2d-10
  RCT(:,18) = 2.9d-11
  RCT(:,21) = 1d-24
  RCT(:,25) = 2.4d-13
  RCT(:,44) = 1d-15
  RCT(:,45) = 3d-12
  RCT(:,58) = 5.3d-12
  RCT(:,63) = 1.7d-11
  RCT(:,72) = 1d-15
  RCT(:,73) = 3d-12
  RCT(:,76) = 1d-15
  RCT(:,77) = 3d-12
  RCT(:,80) = 1d-13
  RCT(:,81) = 3d-12
  RCT(:,90) = 1d-13
  RCT(:,91) = 1d-13
  RCT(:,92) = 2.3d-17
  RCT(:,93) = 2.3d-17
  RCT(:,94) = 2.5d-14
  RCT(:,95) = 2.5d-14
  RCT(:,96) = 2.5d-14
  RCT(:,97) = 2.5d-14
  RCT(:,100) = 3d-11
  RCT(:,102) = 4d-11
  RCT(:,103) = 2.2d-11
  RCT(:,104) = 1.5d-11
  RCT(:,106) = 0.001d0
  RCT(:,109) = 1.8d-11
  RCT(:,110) = 3.2d-13

! Evaluate the bg product coefficients
  npt = gni * gnk
  call mach_adom2_bcs_coeffs(temp2d, p2d, bcs, npt, gni, gnk)
!
! Append the BCS product coefficients to the reaction rates, for use in the
! KPP solver (Fun and Jacobian subroutines)
  do jj = 1, nprcf
     nr = jj + nreact_kpp
     RCT(:, nr) = bcs(:, jj)
  end do

 CONTAINS

! Begin Rate Law Functions from KPP_HOME/util/UserRateLaws

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  User-defined Rate Law functions
!  Note: the default argument type for rate laws, as read from the equations file, is single precision
!        but all the internal calculations are performed in double precision
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~> Simplified Arrhenius, with two arguments
!~~~> Note: The argument B0 has a changed sign when compared to ARR
   real(kind=8) FUNCTION ARR2(A0, B0, TEMP)
      real(kind=8) A0,B0
      real(kind=8) TEMP
      ARR2 =  A0 * EXP(B0 / TEMP)
   END FUNCTION ARR2

   real(kind=8) FUNCTION TYPE5(A1, B1, C1, D1, E1, TEMP, CFACTOR)
      real(kind=8) :: A1, B1, C1, D1, E1
      real(kind=8) :: CFACTOR, TEMP
      real(kind=8) :: RK0, RKIF, RKM, RKMOK, EE
      RK0 = B1 * (TEMP**C1)             !B(T^C)
      RKIF = D1 * (TEMP**E1)            !D(T^E)
      RKM = RK0 * CFACTOR * 1.0d6       !B(T^C)M
      RKMOK = RKM / RKIF                ![B(T^C)M]/[D(T^E)]
      EE  = 1.0d0 / (1.0d0 + (LOG10(RKMOK))**2.0d0)    !EE
      TYPE5 = (RKM / (1.0d0 + RKMOK)) * (A1)**EE
   END FUNCTION TYPE5

end subroutine mach_kpp_adom2_rates
