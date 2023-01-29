"""
Implements backward equations for IF97 region 2.

For the calculation of properties as function of p,h or of p,s without any 
iteration, the two backward equations require extremely good numerical 
consistency with the basic equation. The backward equations for this region
are expressed in dimensionless form, theta = T / Ts, and read as follows, 
where pi = P / Ps, theta = T / Ts, eta = h / hs, and sigma = s / ss.

    a) T / Ts = theta(pi, eta  ) = SUM( ni pi^Ii (eta - 2.1)^Ji )
       T / Ts = theta(pi, sigma) = SUM( ni pi^Ii (sigma - 2)^Ji )

    b) T / Ts = theta(pi, eta  ) = SUM( ni (pi - 2)^Ii (eta - 2.6)^Ji )
       T / Ts = theta(pi, sigma) = SUM( ni  pi^Ii (10 - sigma)^Ji )

    c) T / Ts = theta(pi, eta  ) = SUM( ni (pi + 25)^Ii (eta - 1.8)^Ji )
       T / Ts = theta(pi, sigma) = SUM( ni  pi^Ii (2 - sigma)^Ji )

This module is restricted to providing the dimensionless temperature, 
theta(pi, [eta, sigma]) in region 2. The principle reference for this module
is the IAPWS Industrial Formulation 1997:
    
    International Association for the Properties of Water and Steam, IAPWS R7-97(2012),
    Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic 
    Properties of Water and Steam (2012), available from: http://www.iapws.org.

The dimensional forward and backward functions and their derivatives, 
f(P, T) and f(P, [h, s]) respectively, are defined elsewhere in this package.
"""

# type annotations
from __future__ import annotations

###########################################################
#####       Range of Validity (Boundary Constants)    #####
###########################################################
Pbnd0 = 6.54670 # [MPa      ]
Pbnd1 = 100.0   # [MPa      ]
Tbnd0 = 554.485 # [K        ]
Tbnd1 = 1019.32 # [K        ]
Pbnd2a2b = 4.0  # [MPa      ]
sbnd2b2c = 5.85 # [kJ / kg K]

###########################################################
#####       Subregion Identification Functions        #####
###########################################################

# Constants and domain scale factors
_n = [ 0.90584275814723e+03, -0.67955786399241e+00,  0.12809002730136e-03,
       0.26526571908428e+04,  0.45257578905948e+01]
_Ps = 1.0 # [MPa     ]
_hs = 1.0 # [kJ / kg ]

# Auxiliary functions
def bnd2b2cP(h: float) -> float:
    """Auxiliary equation for subregion 2b and 2c boundary pressure [MPa].
    Reference: Equation (20) from R7-97(2012)"""
    eta = h / _hs
    return _Ps * (_n[0] + _n[1] * eta + _n[2] * eta**2)

def bnd2b2ch(P: float) -> float:
    """Auxiliary equation for subregion 2b and 2c boundary enthalpy [kJ / kg].
    Reference: Equation (21) from R7-97(2012)"""
    pi = P / _Ps
    return _hs * (_n[3] + ((pi - _n[4]) / _n[2])**0.5)

# Subregion identification functions
def region_h(P: float, h: float) -> int:
    """Identification of subregion of region 2 from IF97 specification,
    using pressure and enthalpy as primary variables.
    Reference: Figure (2) from R7-97(2012)"""
    if P <= Pbnd2a2b:
        return 1
    elif (P <= Pbnd0) or (h >= bnd2b2ch(P)):
        return 2
    else:
        return 3

def region_s(P: float, s: float) -> int:
    """Identification of subregion of region 2 from IF97 specification,
    using pressure and entropy as primary variables.
    Reference: Figure (2) from R7-97(2012)"""
    if P <= Pbnd2a2b:
        return 1
    elif (P <= Pbnd0) or (s >= sbnd2b2c):
        return 2
    else:
        return 3

###########################################################
#####       Constants and Auxiliary Functions         #####
###########################################################

# Region 2, backwards equations for f(P, h)
Ps_h = 1.0  # [Mpa    ]
Ts_h = 1.0  # [K      ]
hs_h = 2000 # [kJ / kg]

def theta2a_h(pi: float, eta: float) -> float:
    """Dimensionless temperature as a function of pressure and enthalpy (2a).
    Reference: Equation (22) from R7-97(2012)"""
    I = [ 0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,
          2,  2,  2,  2,  2,  2,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  7]
    J = [ 0,  1,  2,  3,  7, 20,  0,  1,  2,  3,  7,  9, 11, 18, 44,  0,  2,
          7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
    n = [ 0.10898952318288e+04,  0.84951654495535e+03, -0.10781748091826e+03,  0.33153654801263e+02,
         -0.74232016790248e+01,  0.11765048724356e+02,  0.18445749355790e+01, -0.41792700549624e+01,
          0.62478196935812e+01, -0.17344563108114e+02, -0.20058176862096e+03,  0.27196065473796e+03,
         -0.45511318285818e+03,  0.30919688604755e+04,  0.25226640357872e+06, -0.61707422868339e-02,
         -0.31078046629583e+00,  0.11670873077107e+02,  0.12812798404046e+09, -0.98554909623276e+09,
          0.28224546973002e+10, -0.35948971410703e+10,  0.17227349913197e+10, -0.13551334240775e+05,
          0.12848734664650e+08,  0.13865724283226e+01,  0.23598832556514e+06, -0.13105236545054e+08,
          0.73999835474766e+04, -0.55196697030060e+06,  0.37154085996233e+07,  0.19127729239660e+05,
         -0.41535164835634e+06, -0.62459855192507e+02]
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (eta - 2.1)**Ji
    return sum

def theta2b_h(pi: float, eta: float) -> float:
    """Dimensionless temperature as a function of pressure and enthalpy (2b).
    Reference: Equation (23) from R7-97(2012)"""
    I = [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,
          2,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  5,  5,  5,  6,  7,  7,  9,  9]
    J = [ 0,  1,  2, 12, 18, 24, 28, 40,  0,  2,  6, 12, 18, 24, 28, 40,  2,  8, 18,
         40,  1,  2, 12, 24,  2, 12, 18, 24, 28, 40, 18, 24, 40, 28,  2, 28,  1, 40]
    n = [ 0.14895041079516e+04,  0.74307798314034e+03, -0.97708318797837e+02,  0.24742464705674e+01,
         -0.63281320016026e+00,  0.11385952129658e+01, -0.47811863648625e+00,  0.85208123431544e-02,
          0.93747147377932e+00,  0.33593118604916e+01,  0.33809355601454e+01,  0.16844539671904e+00,
          0.73875745236695e+00, -0.47128737436186e+00,  0.15020273139707e+00, -0.21764114219750e-02,
         -0.21810755324761e-01, -0.10829784403677e+00, -0.46333324635812e-01,  0.71280351959551e-04,
          0.11032831789999e-03,  0.18955248387902e-03,  0.30891541160537e-02,  0.13555504554949e-02,
          0.28640237477456e-06, -0.10779857357512e-04, -0.76462712454814e-04,  0.14052392818316e-04,
         -0.31083814331434e-04, -0.10302738212103e-05,  0.28217281635040e-06,  0.12704902271945e-05,
          0.73803353468292e-07, -0.11030139238909e-07, -0.81456365207833e-13, -0.25180545682962e-10,
         -0.17565233969407e-17,  0.86934156344163e-14]
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (pi - 2)**Ii * (eta - 2.6)**Ji
    return sum

def theta2c_h(pi: float, eta: float) -> float:
    """Dimensionless temperature as a function of pressure and enthalpy (2c).
    Reference: Equation (24) from R7-97(2012)"""
   Ic_h = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1,  0,  0,  1,  1,  2,  6,  6,  6,  6,  6,  6,  6,  6]
   Jc_h = [ 0,  4,  0,  2,  0,  2,  0,  1,  0,  2,  0,  1,  4,  8,  4,  0,  1,  4, 10, 12, 16, 20, 22]
   nc_h = [-0.32368398555242e+13,  0.73263350902181e+13,  0.35825089945447e+12, -0.58340131851590e+12,
           -0.10783068217470e+11,  0.20825544563171e+11,  0.61074783564516e+06,  0.85977722535580e+06,
           -0.25745723604170e+05,  0.31081088422714e+05,  0.12082315865936e+04,  0.48219755109255e+03,
            0.37966001272486e+01, -0.10842984880077e+02, -0.45364172676660e-01,  0.14559115658698e-12,
            0.11261597407230e-11, -0.17804982240686e-10,  0.12324579690832e-06, -0.11606921130984e-05,
            0.27846367088554e-04, -0.59270038474176e-03,  0.12918582991878e-02]
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (pi + 25)**Ii * (eta - 1.8)**Ji
    return sum

# Region 2, backwards equations for f(P, s)
Ps_s  = 1.0    # [Mpa      ]
Ts_s  = 1.0    # [K        ]
ss_sa = 2.0    # [kJ / kg K]
ss_sb = 0.7853 # [kJ / kg K]
ss_sc = 2.9251 # [kJ / kg K]


def theta2a_s(pi: float, sigma: float) -> float:
    """Dimensionless temperature as a function of pressure and entropy (2a).
    Reference: Equation (25) from R7-97(2012)"""
    I = [ -1.5,  -1.5,  -1.5,  -1.5,  -1.5,  -1.5, -1.25, -1.25, -1.25,  -1.0,  -1.0,  -1.0,  -1.0,
          -1.0,  -1.0, -0.75, -0.75,  -0.5,  -0.5,  -0.5,  -0.5, -0.25, -0.25, -0.25, -0.25,  0.25,  
          0.25,  0.25,  0.25,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,  0.75,  0.75,  0.75,
          0.75,   1.0,   1.0,  1.25,  1.25,   1.5,   1.5]
    J = [  -24,   -23,   -19,   -13,   -11,   -10,   -19,   -15,    -6,   -26,   -21,   -17,   -16,
            -9,    -8,   -15,   -14,   -26,   -13,    -9,    -7,   -27,   -25,   -11,    -6,     1,
             4,     8,    11,     0,     1,     5,     6,    10,    14,    16,     0,     4,     9,
            17,     7,    18,     3,    15,     5,    18]
    n = [-0.39235983861984e+06,  0.51526573827270e+06,  0.40482443161048e+05, -0.32193790923902e+03,
          0.96961424218694e+02, -0.22867846371773e+02, -0.44942914124357e+06, -0.50118336020166e+04,
          0.35684463560015e+00,  0.44235335848190e+05, -0.13673388811708e+05,  0.42163260207864e+06,
          0.22516925837475e+05,  0.47442144865646e+03, -0.14931130797647e+03, -0.19781126320452e+06, 
         -0.23554399470760e+05, -0.19070616302076e+05,  0.55375669883164e+05,  0.38293691437363e+04,
         -0.60391860580567e+03,  0.19363102620331e+04,  0.42660643698610e+04, -0.59780638872718e+04,
         -0.70401463926862e+03,  0.33836784107553e+03,  0.20862786635187e+02,  0.33834172656196e-01,
         -0.43124428414893e-04,  0.16653791356412e+03, -0.13986292055898e+03, -0.78849547999872e+00,
          0.72132411753872e-01, -0.59754839398283e-02, -0.12141358953904e-04,  0.23227096733871e-06,
         -0.10538463566194e+02,  0.20718925496502e+01, -0.72193155260427e-01,  0.20749887081120e-06,
         -0.18340657911379e-01,  0.29036272348696e-06,  0.21037527893619e+00,  0.25681239729999e-03,
         -0.12799002933781e-01, -0.82198102652018e-05]
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (sigma - 2)**Ji
    return sum

def theta2b_s(pi: float, sigma: float) -> float:
    """Dimensionless temperature as a function of pressure and entropy (2b).
    Reference: Equation (26) from R7-97(2012)"""
    I = [-6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0,
          0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  5,  5,  5]
    J = [ 0, 11,  0, 11,  0,  1, 11,  0,  1, 11, 12,  0,  1,  6, 10,  0,  1,  5,  8,  9,  0,  1,
          2,  4,  5,  6,  9,  0,  1,  2,  3,  7,  8,  0,  1,  5,  0,  1,  3,  0,  1,  0,  1,  2]
    n = [ 0.31687665083497e+06,  0.20864175881858e+02, -0.39859399803599e+06, -0.21816058518877e+02,
          0.22369785194242e+06, -0.27841703445817e+04,  0.99207436071480e+01, -0.75197512299157e+05,
          0.29708605951158e+04, -0.34406878548526e+01,  0.38815564249115e+00,  0.17511295085750e+05,
         -0.14237112854449e+04,  0.10943803364167e+01,  0.89971619308495e+00, -0.33759740098958e+04,
          0.47162885818355e+03, -0.19188241993679e+01,  0.41078580492196e+00, -0.33465378172097e+00,
          0.13870034777505e+04, -0.40663326195838e+03,  0.41727347159610e+02,  0.21932549434532e+01,
         -0.10320050009077e+01,  0.35882943516703e+00,  0.52511453726066e-02,  0.12838916450705e+02,
         -0.28642437219381e+01,  0.56912683664855e+00, -0.99962954584931e-01, -0.32632037778459e-02,
          0.23320922576723e-03, -0.15334809857450e+00,  0.29072288239902e-01,  0.37534702741167e-03,
          0.17296691702411e-02, -0.38556050844504e-03, -0.35017712292608e-04, -0.14566393631492e-04,
          0.56420857267269e-05,  0.41286150074605e-07, -0.20684671118824e-07,  0.16409393674725e-08]
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (10 - sigma)**Ji
    return sum

def theta2c_s(pi: float, sigma: float) -> float:
    """Dimensionless temperature as a function of pressure and entropy (2c).
    Reference: Equation (27) from R7-97(2012)"""
    I = [-2, -2, -1,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,
          5,  6,  6,  7,  7,  7,  7,  7]
    J = [ 0,  1,  0,  0,  1,  2,  3,  0,  1,  3,  4,  0,  1,  2,  0,  1,  5,  0,  1,  4,  0,  1,
          2,  0,  1,  0,  1,  3,  4,  5]
    n = [ 0.90968501005365e+03,  0.24045667088420e+04, -0.59162326387130e+03,  0.54145404128074e+03,
         -0.27098308411192e+03,  0.97976525097926e+03, -0.46966772959435e+03,  0.14399274604723e+02,
         -0.19104204230429e+02,  0.53299167111971e+01, -0.21252975375934e+02, -0.31147334413760e+00,
          0.60334840894623e+00, -0.42764839702509e-01,  0.58185597255259e-02, -0.14597008284753e-01,
          0.56631175631027e-02, -0.76155864584577e-04,  0.22440342919332e-03, -0.12561095013413e-04,
          0.63323132660934e-06, -0.20541989675375e-05,  0.36405370390082e-07, -0.29759897789215e-08,
          0.10136618529763e-07,  0.59925719692351e-11, -0.20677870105164e-10, -0.20874278181886e-10,
          0.10162166825089e-09, -0.16429828281347e-09]

    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (2 - sigma)**Ji
    return sum
