! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Function of Chemical Model File
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
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!============================================================================!
!
! GEM-MACH Setup
! Fichier / File   : mach_kpp_saprc07cs_Function.ftn90
! Creation         : Diane Pendlebury (Feb. 2018) and Jack Chen (Feb. 2017)
! Modifications    : code adapted from KPP
!
!============================================================================

MODULE mach_kpp_saprc07cs_Function

  USE mach_pkg_gas_mod, ONLY: NVAR, NFIX, NREACT
  USE chm_utils_mod,    ONLY: DP
  IMPLICIT NONE

 CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

! Computation of equation rates
  A(1) = RCT(1)*V(51)
  A(2) = RCT(2)*V(39)*F(1)*F(3)
  A(3) = RCT(3)*V(39)*V(41)
  A(4) = RCT(4)*V(39)*V(50)
  A(5) = RCT(5)*V(39)*V(51)
  A(6) = RCT(6)*V(39)*V(51)
  A(7) = RCT(7)*V(41)*V(50)
  A(8) = RCT(8)*V(41)*V(51)
  A(9) = RCT(9)*V(48)*V(50)
  A(10) = RCT(10)*V(50)*V(50)*F(3)
  A(11) = RCT(11)*V(48)*V(51)
  A(12) = RCT(12)*V(19)
  A(13) = RCT(13)*V(19)*F(4)
  A(14) = RCT(14)*V(19)*F(4)*F(4)
  A(15) = RCT(15)*V(48)*V(51)
  A(16) = RCT(16)*V(48)
  A(17) = RCT(17)*V(48)
  A(18) = RCT(18)*V(41)
  A(19) = RCT(19)*V(41)
  A(20) = RCT(20)*V(21)*F(4)
  A(21) = RCT(21)*V(21)*F(1)
  A(22) = RCT(22)*V(47)*V(50)
  A(23) = RCT(23)*V(10)
  A(24) = RCT(24)*V(10)*V(47)
  A(25) = RCT(25)*V(47)*V(51)
  A(26) = RCT(26)*V(47)*V(48)
  A(27) = RCT(27)*V(25)*V(47)
  A(28) = RCT(28)*V(25)
  A(29) = RCT(29)*V(26)*V(47)
  A(30) = RCT(30)*V(41)*V(47)
  A(31) = RCT(31)*V(44)*V(50)
  A(32) = RCT(32)*V(44)*V(51)
  A(33) = RCT(33)*V(20)
  A(34) = RCT(34)*V(20)
  A(35) = RCT(35)*V(20)*V(47)
  A(36) = RCT(36)*V(41)*V(44)
  A(37) = RCT(37)*V(44)*V(44)
  A(38) = RCT(38)*V(44)*V(44)*F(4)
  A(39) = RCT(39)*V(44)*V(48)
  A(40) = RCT(40)*V(48)*V(48)
  A(41) = RCT(41)*V(9)
  A(42) = RCT(42)*V(9)*V(47)
  A(43) = RCT(43)*V(44)*V(47)
  A(44) = RCT(44)*V(6)*V(47)
  A(45) = RCT(45)*V(47)*F(5)
  A(46) = RCT(46)*V(43)*V(50)
  A(47) = RCT(47)*V(43)*V(44)
  A(48) = RCT(48)*V(43)*V(48)
  A(49) = RCT(49)*V(43)*V(43)
  A(50) = RCT(50)*V(45)*V(50)
  A(51) = RCT(51)*V(44)*V(45)
  A(52) = RCT(52)*V(45)*V(48)
  A(53) = RCT(53)*V(43)*V(45)
  A(54) = RCT(54)*V(45)*V(45)
  A(55) = RCT(55)*V(46)*V(51)
  A(56) = RCT(56)*V(17)
  A(57) = RCT(57)*V(17)
  A(58) = RCT(58)*V(46)*V(50)
  A(59) = RCT(59)*V(44)*V(46)
  A(60) = RCT(60)*V(46)*V(48)
  A(61) = RCT(61)*V(43)*V(46)
  A(62) = RCT(62)*V(45)*V(46)
  A(63) = RCT(63)*V(46)*V(46)
  A(64) = RCT(64)*V(51)*V(52)
  A(65) = RCT(65)*V(18)
  A(66) = RCT(66)*V(18)
  A(67) = RCT(67)*V(50)*V(52)
  A(68) = RCT(68)*V(44)*V(52)
  A(69) = RCT(69)*V(48)*V(52)
  A(70) = RCT(70)*V(43)*V(52)
  A(71) = RCT(71)*V(45)*V(52)
  A(72) = RCT(72)*V(46)*V(52)
  A(73) = RCT(73)*V(52)*V(52)
  A(74) = RCT(74)*V(22)*V(51)
  A(75) = RCT(75)*V(22)*V(44)
  A(76) = RCT(76)*V(22)
  A(77) = RCT(77)*V(30)*V(50)
  A(78) = RCT(78)*V(30)*V(44)
  A(79) = RCT(79)*V(30)*V(48)
  A(80) = RCT(80)*V(30)*V(43)
  A(81) = RCT(81)*V(30)*V(45)
  A(82) = RCT(82)*V(30)*V(30)
  A(83) = RCT(83)*V(30)*V(46)
  A(84) = RCT(84)*V(30)*V(52)
  A(85) = RCT(85)*V(33)
  A(86) = RCT(86)*V(33)
  A(87) = RCT(87)*V(33)*V(47)
  A(88) = RCT(88)*V(33)*V(48)
  A(89) = RCT(89)*V(40)*V(47)
  A(90) = RCT(90)*V(40)
  A(91) = RCT(91)*V(40)*V(48)
  A(92) = RCT(92)*V(35)*V(47)
  A(93) = RCT(93)*V(35)
  A(94) = RCT(94)*V(35)*V(48)
  A(95) = RCT(95)*V(12)*V(47)
  A(96) = RCT(96)*V(12)
  A(97) = RCT(97)*V(29)
  A(98) = RCT(98)*V(29)*V(47)
  A(99) = RCT(99)*V(29)*V(48)
  A(100) = RCT(100)*V(24)*V(47)
  A(101) = RCT(101)*V(24)*V(48)
  A(102) = RCT(102)*V(27)*V(47)
  A(103) = RCT(103)*V(27)*V(41)
  A(104) = RCT(104)*V(27)
  A(105) = RCT(105)*V(23)*V(47)
  A(106) = RCT(106)*V(23)*V(41)
  A(107) = RCT(107)*V(23)
  A(108) = RCT(108)*V(37)*V(47)
  A(109) = RCT(109)*V(37)*V(41)
  A(110) = RCT(110)*V(37)*V(48)
  A(111) = RCT(111)*V(37)
  A(112) = RCT(112)*V(47)*V(49)
  A(113) = RCT(113)*V(49)
  A(114) = RCT(114)*V(38)*V(47)
  A(115) = RCT(115)*V(38)
  A(116) = RCT(116)*V(47)*F(6)
  A(117) = RCT(117)*V(28)*V(47)
  A(118) = RCT(118)*V(28)*V(41)
  A(119) = RCT(119)*V(28)*V(48)
  A(120) = RCT(120)*V(28)*V(39)
  A(121) = RCT(121)*V(31)*V(47)
  A(122) = RCT(122)*V(31)*V(41)
  A(123) = RCT(123)*V(31)*V(48)
  A(124) = RCT(124)*V(31)*V(39)
  A(125) = RCT(125)*V(42)*V(47)
  A(126) = RCT(126)*V(16)*V(47)
  A(127) = RCT(127)*V(32)*V(47)
  A(128) = RCT(128)*V(32)*V(41)
  A(129) = RCT(129)*V(32)*V(48)
  A(130) = RCT(130)*V(32)*V(39)
  A(131) = RCT(131)*V(36)*V(47)
  A(132) = RCT(132)*V(36)*V(41)
  A(133) = RCT(133)*V(36)*V(48)
  A(134) = RCT(134)*V(36)*V(39)
  A(135) = RCT(135)*V(14)*V(47)
  A(136) = RCT(136)*V(15)*V(47)
  A(137) = RCT(137)*V(34)*V(47)
  A(138) = RCT(138)*V(34)*V(41)
  A(139) = RCT(139)*V(34)*V(48)
  A(140) = RCT(140)*V(34)*V(39)
  A(141) = RCT(141)*V(8)*V(50)
  A(142) = RCT(142)*V(8)*V(44)
  A(143) = RCT(143)*V(7)*V(47)
  A(144) = RCT(144)*V(7)*V(47)
  A(145) = RCT(145)*V(7)
  A(146) = RCT(146)*V(1)*V(47)
  A(147) = RCT(147)*V(50)
  A(148) = RCT(148)*V(11)
  A(149) = RCT(149)*V(19)
  A(150) = RCT(150)*V(19)
  A(151) = RCT(151)*F(3)
  A(152) = RCT(152)*F(6)
  A(153) = RCT(153)*V(13)*V(50)
  A(154) = RCT(154)*V(13)*F(3)
  A(155) = RCT(155)*V(13)*V(51)
  A(156) = RCT(156)*V(21)*V(41)
  A(157) = RCT(157)*V(11)*V(21)
  A(158) = RCT(158)*V(11)*V(21)

! Aggregate function
  Vdot(1) = A(143)-A(146)
  Vdot(2) = A(44)
  Vdot(3) = A(29)+0.4*A(57)+A(58)+A(60)+A(61)+A(62)+2*A(63)+0.4*A(66)+A(67)+A(69)+A(70)+A(71)+2*A(72)+2*A(73)+0.174&
              &*A(103)+0.174*A(106)+0.13*A(109)+0.12*A(118)+0.122*A(122)+0.125*A(128)+0.162*A(132)+0.045*A(138)
  Vdot(4) = A(74)+0.283*A(110)+0.174*A(114)+A(119)+0.813*A(123)+0.454*A(129)+0.269*A(133)+0.416*A(139)
  Vdot(5) = 0.856*A(59)+0.516*A(68)+2.5*A(74)-A(75)-A(76)+5.05*A(100)+A(101)+0.747*A(102)-0.796*A(103)+0.753*A(104)&
              &+0.747*A(105)-0.796*A(106)-A(107)+0.586*A(108)+0.842*A(109)-0.583*A(110)+0.123*A(111)+1.085*A(112)-0.091&
              &*A(113)+0.37*A(114)+0.786*A(115)+0.298*A(118)-A(119)+0.203*A(120)-0.717*A(121)-0.753*A(122)-0.064*A(123)-0.77&
              &*A(124)+1.216*A(125)+0.468*A(126)+0.803*A(127)+0.935*A(128)+0.632*A(129)+2.006*A(130)+0.566*A(131)+0.453&
              &*A(132)+0.372*A(133)+1.931*A(134)+0.453*A(135)+1.807*A(136)+5.154*A(137)+4.045*A(138)+4.922*A(139)+4.441&
              &*A(140)
  Vdot(6) = -A(44)
  Vdot(7) = A(142)-A(143)-A(144)-A(145)
  Vdot(8) = A(121)-A(141)-A(142)+0.387*A(144)
  Vdot(9) = A(37)+A(38)-A(41)-A(42)
  Vdot(10) = A(22)-A(23)-A(24)
  Vdot(11) = -A(148)+A(155)-A(157)-A(158)
  Vdot(12) = A(47)+A(51)-A(95)-A(96)
  Vdot(13) = A(147)-A(153)-A(154)-A(155)
  Vdot(14) = -A(135)
  Vdot(15) = -A(136)
  Vdot(16) = -A(126)
  Vdot(17) = A(55)-A(56)-A(57)
  Vdot(18) = A(64)-A(65)-A(66)
  Vdot(19) = A(11)-A(12)-A(13)-A(14)-A(149)-A(150)
  Vdot(20) = A(32)-A(33)-A(34)-A(35)
  Vdot(21) = A(18)-A(20)-A(21)+A(147)+A(148)-A(156)-A(157)-A(158)
  Vdot(22) = -A(74)-A(75)-A(76)+0.2*A(100)+A(101)
  Vdot(23) = -A(105)-A(106)-A(107)+0.454*A(135)+0.469*A(136)
  Vdot(24) = 0.5*A(74)+A(75)+A(76)-A(100)-A(101)+0.166*A(135)+0.108*A(136)
  Vdot(25) = 2*A(13)+2*A(14)+A(25)-A(27)-A(28)+0.2*A(39)+A(88)+A(91)+A(94)+A(99)+A(101)+0.158*A(110)
  Vdot(26) = -A(29)+A(85)+A(86)+A(87)+A(88)+A(90)+0.035*A(92)+A(93)+A(97)+A(98)+A(99)+0.334*A(102)+0.522*A(103)+0.695&
               &*A(104)+0.334*A(105)+0.522*A(106)+0.238*A(108)+0.483*A(109)+0.57*A(110)+0.867*A(111)+0.51*A(118)+0.788&
               &*A(120)+0.275*A(122)+0.001*A(126)+0.368*A(128)+0.297*A(132)+0.006*A(134)+0.001*A(137)+0.185*A(138)+0.01&
               &*A(139)+A(152)
  Vdot(27) = -0.998*A(102)-0.993*A(103)-0.997*A(104)+0.002*A(105)+0.007*A(106)+0.001*A(108)+0.167*A(135)+0.219*A(136)
  Vdot(28) = -A(117)-A(118)-A(119)-A(120)
  Vdot(29) = -A(97)-A(98)-A(99)+0.25*A(100)+0.329*A(102)+0.819*A(103)+0.418*A(104)+0.329*A(105)+0.819*A(106)+0.23*A(108)&
               &+0.863*A(109)+0.008*A(110)+0.007*A(120)+0.47*A(135)+0.709*A(136)+0.029*A(137)+0.083*A(138)+0.001*A(139)
  Vdot(30) = -A(77)-A(78)-A(79)-A(80)-A(81)-2*A(82)-A(83)-A(84)+0.202*A(102)+0.652*A(103)+0.202*A(105)+0.652*A(106)&
               &+0.238*A(108)+0.028*A(109)+0.084*A(111)+0.094*A(112)+0.677*A(113)+0.671*A(114)+0.167*A(115)+0.079*A(121)&
               &+0.192*A(122)+0.187*A(123)+0.24*A(124)+0.558*A(125)+0.943*A(126)+0.234*A(127)+0.023*A(128)+0.488*A(129)&
               &+0.052*A(131)+0.161*A(132)+0.762*A(133)+0.006*A(134)+0.388*A(137)+0.808*A(138)+1.347*A(139)
  Vdot(31) = -A(121)-A(122)-A(123)-A(124)
  Vdot(32) = -A(127)-A(128)-A(129)-A(130)
  Vdot(33) = 0.4*A(57)+A(58)+A(60)+A(61)+A(62)+2*A(63)+A(72)-A(85)-A(86)-A(87)-A(88)+A(90)+0.652*A(103)+0.173*A(104)&
               &+0.652*A(106)+0.149*A(108)+0.119*A(109)+0.222*A(110)+0.429*A(111)+0.213*A(112)+0.303*A(113)+0.011*A(114)&
               &+0.135*A(115)+A(116)+1.61*A(117)+A(118)+0.788*A(120)+0.624*A(121)+0.592*A(122)+0.49*A(124)+0.262*A(125)&
               &+0.039*A(126)+0.701*A(127)+0.604*A(128)+0.209*A(131)+0.55*A(132)+0.107*A(133)+0.264*A(137)+0.229*A(138)&
               &+0.017*A(139)
  Vdot(34) = -A(137)-A(138)-A(139)-A(140)
  Vdot(35) = -A(92)-A(93)-A(94)+0.407*A(102)+0.407*A(105)+0.243*A(108)+0.213*A(110)+0.545*A(112)+0.78*A(113)+0.037&
               &*A(114)+0.137*A(115)+A(119)+0.122*A(125)+0.226*A(126)+0.47*A(127)+0.384*A(128)+0.002*A(129)+0.45*A(130)&
               &+0.481*A(131)+0.333*A(132)+0.154*A(133)+0.074*A(134)+0.533*A(137)+0.22*A(138)+0.509*A(139)+0.147*A(140)
  Vdot(36) = -A(131)-A(132)-A(133)-A(134)
  Vdot(37) = -A(108)-A(109)-A(110)-A(111)+0.907*A(121)+0.7*A(122)+0.936*A(123)+0.04*A(127)+0.066*A(131)+0.049*A(132)&
               &+0.033*A(133)+0.006*A(134)+0.003*A(137)+0.002*A(138)+0.006*A(139)
  Vdot(38) = A(50)+0.559*A(110)-0.506*A(114)-A(115)+0.546*A(129)+0.322*A(133)+0.163*A(139)
  Vdot(39) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(17)+A(19)+A(21)-A(120)-A(124)-A(130)-A(134)-A(140)+A(150)+A(151)+A(153)&
               &+A(154)+A(155)+A(156)
  Vdot(40) = 0.086*A(59)+0.4*A(66)+A(67)+0.19*A(68)+A(69)+A(70)+A(71)+A(72)+2*A(73)-A(89)-A(90)-A(91)+0.035*A(92)+A(93)&
               &+0.051*A(108)+0.05*A(109)+0.184*A(111)+0.084*A(112)+0.163*A(113)+0.429*A(114)+0.444*A(115)+0.195*A(117)+0.1&
               &*A(120)+0.445*A(125)+0.31*A(126)+0.301*A(127)+0.177*A(128)+0.009*A(129)+0.788*A(131)+0.542*A(132)+0.546&
               &*A(133)+0.008*A(138)+0.001*A(139)
  Vdot(41) = A(2)-A(3)-A(7)-A(8)-A(18)-A(19)-A(30)-A(36)+0.3*A(59)+0.25*A(68)-A(103)-A(106)-A(109)-A(118)-A(122)-A(128)&
               &-A(132)-A(138)-A(156)
  Vdot(42) = 0.243*A(59)+0.526*A(68)+0.095*A(109)+0.018*A(118)+0.01*A(122)-A(125)+0.073*A(128)+0.133*A(132)+0.028*A(138)
  Vdot(43) = -A(46)-A(47)-A(48)-2*A(49)-A(53)+0.4*A(57)+A(58)+A(60)+A(62)+2*A(63)+0.4*A(66)+A(67)+A(69)+A(71)+2*A(72)+2&
               &*A(73)+A(76)+A(90)+0.035*A(92)+A(93)+0.8*A(100)+0.521*A(102)+0.173*A(104)+0.521*A(105)+0.496*A(108)+0.025&
               &*A(109)+0.792*A(110)+0.141*A(111)+0.379*A(112)+0.913*A(113)+0.305*A(114)+0.554*A(115)+A(116)+A(117)+A(119)&
               &+0.8*A(120)+0.907*A(121)+0.749*A(123)+0.25*A(124)+0.931*A(125)+0.766*A(126)+0.905*A(127)+0.144*A(128)+0.824&
               &*A(129)+0.914*A(131)+0.329*A(132)+0.456*A(133)+0.007*A(134)+0.482*A(135)+0.58*A(136)+0.759*A(137)+0.067&
               &*A(138)+0.162*A(139)
  Vdot(44) = A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)+A(46)&
               &-A(47)+A(48)+A(49)-A(51)+A(52)+A(53)+A(54)-A(59)+A(62)-A(68)+A(71)-A(75)+2*A(85)+A(87)+A(88)+A(90)+A(93)&
               &+A(96)+A(97)+0.522*A(103)+1.023*A(104)+0.522*A(106)+0.211*A(109)+0.655*A(111)+0.472*A(112)+0.189*A(114)&
               &+0.344*A(115)+0.16*A(118)+0.8*A(120)+0.066*A(122)+0.116*A(128)+0.093*A(132)+0.014*A(134)+0.166*A(135)+0.108&
               &*A(136)+0.052*A(138)+2*A(152)
  Vdot(45) = -A(50)-A(51)-A(52)-A(53)-2*A(54)-A(62)-A(71)+0.06*A(102)+0.06*A(105)+0.025*A(108)+0.05*A(110)+0.071*A(112)&
               &+0.087*A(113)+0.175*A(114)+0.102*A(115)+0.093*A(121)+0.008*A(122)+0.064*A(123)+0.01*A(124)+0.07*A(125)+0.227&
               &*A(126)+0.095*A(127)+0.004*A(128)+0.176*A(129)+0.086*A(131)+0.003*A(132)+0.136*A(133)+0.001*A(134)+0.068&
               &*A(135)+0.11*A(136)+0.2*A(137)+0.203*A(138)+0.397*A(139)
  Vdot(46) = -A(55)+A(56)+0.6*A(57)-A(58)-A(59)-A(60)-A(61)-A(62)-2*A(63)-A(72)+A(89)+A(91)+A(97)+A(98)+A(99)+0.201&
               &*A(102)+0.305*A(104)+0.201*A(105)+0.238*A(108)+0.354*A(111)+0.029*A(112)+0.4*A(113)+0.007*A(126)+0.147&
               &*A(132)+0.126*A(138)
  Vdot(47) = 2*A(20)-A(22)+A(23)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
               &*A(41)-A(42)-A(43)-A(44)-A(45)-A(87)-A(89)-A(92)+A(96)-A(98)-A(100)-A(102)+0.826*A(103)-A(105)+0.826*A(106)&
               &-A(108)+0.219*A(109)+0.084*A(111)-A(112)-A(114)-A(116)-A(117)+0.16*A(118)-A(121)+0.266*A(122)-A(125)-A(126)&
               &-A(127)+0.193*A(128)-A(131)+0.423*A(132)-0.716*A(135)-0.798*A(136)-A(137)+0.585*A(138)
  Vdot(48) = A(6)+A(8)-A(9)-A(11)+A(12)-A(15)-A(16)-A(17)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(52)+0.4*A(57)&
               &-A(60)+0.4*A(66)-A(69)-A(79)-A(88)-A(91)-A(94)-A(99)-A(101)-A(110)-A(119)-A(123)-A(129)-A(133)-A(139)+A(149)&
               &+A(150)
  Vdot(49) = A(51)+A(52)+A(53)+2*A(54)+A(62)+A(71)+0.048*A(102)+0.048*A(105)+A(107)+0.192*A(108)+0.033*A(109)+0.246&
               &*A(111)-0.622*A(112)-A(113)+0.106*A(114)+0.528*A(115)+0.1*A(122)+0.75*A(124)+0.141*A(125)+0.303*A(126)+0.119&
               &*A(127)+0.191*A(128)+0.002*A(129)+0.274*A(130)+0.06*A(131)+0.057*A(132)+0.009*A(133)+0.464*A(134)+0.077&
               &*A(135)+0.035*A(136)+0.259*A(137)+0.422*A(138)+0.012*A(139)+0.853*A(140)
  Vdot(50) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(15)+A(16)-A(22)+A(23)-A(31)-A(46)-A(50)-A(58)-A(67)-A(77)-A(147)+A(150)&
               &-A(153)+A(154)+2*A(157)
  Vdot(51) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(17)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)+A(33)&
               &+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(52)-A(55)+A(56)+0.6*A(57)+A(58)+A(60)-A(64)+A(65)+0.6&
               &*A(66)+A(67)+A(69)-A(74)+A(77)+A(79)+0.332*A(114)+A(115)+0.187*A(123)+0.409*A(133)+0.421*A(139)+A(149)&
               &-A(155)
  Vdot(52) = -A(64)+A(65)+0.6*A(66)-A(67)-A(68)-A(69)-A(70)-A(71)-A(72)-2*A(73)+0.965*A(92)+A(94)+0.217*A(102)+0.652&
               &*A(103)+0.5*A(104)+0.217*A(105)+0.652*A(106)+0.241*A(108)+0.053*A(109)+0.158*A(110)+0.343*A(111)+0.049&
               &*A(112)+0.6*A(113)+0.192*A(122)+0.24*A(124)+0.008*A(132)+0.007*A(134)+0.042*A(137)+0.149*A(138)+0.019*A(139)
      
END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE mach_kpp_saprc07cs_Function

