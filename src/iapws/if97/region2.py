from math import log

###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# constants and non-dimenionalization;
# Region 2, forwards equations for (P, T)
J0 = [0,  1, -5, -4, -3, -2, -1,  2,  3]
n0 = [-0.96927686500217e+01,  0.10086655968018e+02, -0.56087911283020e-02,  0.71452738081455e-01, -0.40710498223928e+00,  0.14240819171444e+01,
      -0.43839511319450e+01, -0.28408632460772e+00,  0.21268463753307e-01]
Ir = [   1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  5,  6,  6,  6,
         7,  7,  7,  8,  8,  9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
Jr = [   0,  1,  2,  3,  6,  1,  2,  4,  7, 36,  0,  1,  3,  6, 35,  1,  2,  3,  7,  3, 16, 35,
         0, 11, 25,  8, 36, 13,  4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
nr = [-0.17731742473213e-02, -0.17834862292358e-01, -0.45996013696365e-01, -0.57581259083432e-01, -0.50325278727930e-01, -0.33032641670203e-04,
      -0.18948987516315e-03, -0.39392777243355e-02, -0.43797295650573e-01, -0.26674547914087e-04,  0.20481737692309e-07,  0.43870667284435e-06,
      -0.32277677238570e-04, -0.15033924542148e-02, -0.40668253562649e-01, -0.78847309559367e-09,  0.12790717852285e-07,  0.48225372718507e-06,
       0.22922076337661e-05, -0.16714766451061e-10, -0.21171472321355e-02, -0.23895741934104e+02, -0.59059564324270e-17, -0.12621808899101e-05,
      -0.38946842435739e-01,  0.11256211360459e-10, -0.82311340897998e+01,  0.19809712802088e-07,  0.10406965210174e-18, -0.10234747095929e-12,
      -0.10018179379511e-08, -0.80882908646985e-10,  0.10693031879409e+00, -0.33662250574171e+00,  0.89185845355421e-24,  0.30629316876232e-12,
      -0.42002467698208e-05, -0.59056029685639e-25,  0.37826947613457e-05, -0.12768608934681e-14,  0.73087610595061e-28,  0.55414715350778e-16,
      -0.94369707241210e-06]
Ps = 1.0      # [Mpa      ]
Ts = 540.0    # [K        ]
R  = 0.461526 # [kJ / kg K]

# constants for subregion boundaries;
n  = [ 0.90584275814723e+03, -0.67955786399241e+00,  0.12809002730136e-03,  0.26526571908428e+04,  0.45257578905948e+01]
Ps_br = 1.0 # [Mpa    ]
hs_br = 1.0 # [kJ / kg]

# constants and non-dimenionalization;
# Region 2, backwards equations for (P, h)
Ia_h = [0, 0, 0, 0, 0,  0, 1, 1, 1, 1, 1, 1,  1,  1,  1, 2, 2, 2,  2,  2,  2,  2,  2,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  7]
Ja_h = [0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
na_h = [ 0.10898952318288e+04,  0.84951654495535e+03, -0.10781748091826e+03,  0.33153654801263e+02, -0.74232016790248e+01,  0.11765048724356e+02,
         0.18445749355790e+01, -0.41792700549624e+01,  0.62478196935812e+01, -0.17344563108114e+02, -0.20058176862096e+03,  0.27196065473796e+03,
        -0.45511318285818e+03,  0.30919688604755e+04,  0.25226640357872e+06, -0.61707422868339e-02, -0.31078046629583e+00,  0.11670873077107e+02,
         0.12812798404046e+09, -0.98554909623276e+09,  0.28224546973002e+10, -0.35948971410703e+10,  0.17227349913197e+10, -0.13551334240775e+05,
         0.12848734664650e+08,  0.13865724283226e+01,  0.23598832556514e+06, -0.13105236545054e+08,  0.73999835474766e+04, -0.55196697030060e+06,
         0.37154085996233e+07,  0.19127729239660e+05, -0.41535164835634e+06, -0.62459855192507e+02]
Ib_h = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9]
Jb_h = [0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40]
nb_h = [ 0.14895041079516e+04,  0.74307798314034e+03, -0.97708318797837e+02,  0.24742464705674e+01, -0.63281320016026e+00,  0.11385952129658e+01,
        -0.47811863648625e+00,  0.85208123431544e-02,  0.93747147377932e+00,  0.33593118604916e+01,  0.33809355601454e+01,  0.16844539671904e+00,
         0.73875745236695e+00, -0.47128737436186e+00,  0.15020273139707e+00, -0.21764114219750e-02, -0.21810755324761e-01, -0.10829784403677e+00,
        -0.46333324635812e-01,  0.71280351959551e-04,  0.11032831789999e-03,  0.18955248387902e-03,  0.30891541160537e-02,  0.13555504554949e-02,
         0.28640237477456e-06, -0.10779857357512e-04, -0.76462712454814e-04,  0.14052392818316e-04, -0.31083814331434e-04, -0.10302738212103e-05,
         0.28217281635040e-06,  0.12704902271945e-05,  0.73803353468292e-07, -0.11030139238909e-07, -0.81456365207833e-13, -0.25180545682962e-10,
        -0.17565233969407e-17,  0.86934156344163e-14]
Ic_h = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6,  6,  6,  6,  6,  6]
Jc_h = [ 0,  4,  0,  2,  0,  2,  0,  1,  0,  2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22]
nc_h = [-0.32368398555242e+13,  0.73263350902181e+13,  0.35825089945447e+12, -0.58340131851590e+12, -0.10783068217470e+11,  0.20825544563171e+11,
         0.61074783564516e+06,  0.85977722535580e+06, -0.25745723604170e+05,  0.31081088422714e+05,  0.12082315865936e+04,  0.48219755109255e+03,
         0.37966001272486e+01, -0.10842984880077e+02, -0.45364172676660e-01,  0.14559115658698e-12,  0.11261597407230e-11, -0.17804982240686e-10,
         0.12324579690832e-06, -0.11606921130984e-05,  0.27846367088554e-04, -0.59270038474176e-03,  0.12918582991878e-02]
Ps_bh = 1.0    # [Mpa]
Ts_bh = 1.0    # [K   ]
hs_bh = 2000.0 # [kJ / kg]

# constants and non-dimenionalization;
# Region 2, backwards equations for (P, s)
Ia_s = [ -1.5,  -1.5,  -1.5,  -1.5,  -1.5,  -1.5, -1.25, -1.25, -1.25,  -1.0,  -1.0,  -1.0,  -1.0,  -1.0,  -1.0, -0.75, -0.75,  -0.5,
         -0.5,  -0.5,  -0.5, -0.25, -0.25, -0.25, -0.25,  0.25,  0.25,  0.25,  0.25,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5, 
         0.75,  0.75,  0.75,  0.75,   1.0,   1.0,  1.25,  1.25,   1.5,   1.5]
Ja_s = [  -24,   -23,   -19,   -13,   -11,   -10,   -19,   -15,    -6,   -26,   -21,   -17,   -16,    -9,    -8,   -15,   -14,   -26, 
          -13,    -9,    -7,   -27,   -25,   -11,    -6,     1,     4,     8,    11,      0,    1,     5,     6,    10,    14,    16, 
            0,     4,     9,    17,     7,    18,     3,    15,     5,    18]
na_s = [-0.39235983861984e+06,  0.51526573827270e+06,  0.40482443161048e+05, -0.32193790923902e+03,  0.96961424218694e+02, -0.22867846371773e+02,
        -0.44942914124357e+06, -0.50118336020166e+04,  0.35684463560015e+00,  0.44235335848190e+05, -0.13673388811708e+05,  0.42163260207864e+06,
         0.22516925837475e+05,  0.47442144865646e+03, -0.14931130797647e+03, -0.19781126320452e+06, -0.23554399470760e+05, -0.19070616302076e+05,
         0.55375669883164e+05,  0.38293691437363e+04, -0.60391860580567e+03,  0.19363102620331e+04,  0.42660643698610e+04, -0.59780638872718e+04,
        -0.70401463926862e+03,  0.33836784107553e+03,  0.20862786635187e+02,  0.33834172656196e-01, -0.43124428414893e-04,  0.16653791356412e+03,
        -0.13986292055898e+03, -0.78849547999872e+00,  0.72132411753872e-01, -0.59754839398283e-02, -0.12141358953904e-04,  0.23227096733871e-06,
        -0.10538463566194e+02,  0.20718925496502e+01, -0.72193155260427e-01,  0.20749887081120e-06, -0.18340657911379e-01,  0.29036272348696e-06,
         0.21037527893619e+00,  0.25681239729999e-03, -0.12799002933781e-01, -0.82198102652018e-05]
Ib_s = [   -6,    -6,    -5,    -5,    -4,    -4,    -4,    -3,    -3,    -3,    -3,    -2,    -2,    -2,    -2,    -1,    -1,    -1,
           -1,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     1,     1,     1,     1,     2,     2,     2, 
            3,     3,     3,     4,     4,     5,     5,     5]
Jb_s = [    0,    11,     0,    11,     0,     1,    11,     0,     1,    11,    12,     0,     1,     6,    10,     0,     1,     5, 
            8,     9,     0,     1,     2,     4,     5,     6,     9,     0,     1,     2,     3,     7,     8,     0,     1,     5, 
            0,     1,     3,     0,     1,     0,     1,     2]
nb_s = [ 0.31687665083497e+06,  0.20864175881858e+02, -0.39859399803599e+06, -0.21816058518877e+02,  0.22369785194242e+06, -0.27841703445817e+04,
         0.99207436071480e+01, -0.75197512299157e+05,  0.29708605951158e+04, -0.34406878548526e+01,  0.38815564249115e+00,  0.17511295085750e+05,
        -0.14237112854449e+04,  0.10943803364167e+01,  0.89971619308495e+00, -0.33759740098958e+04,  0.47162885818355e+03, -0.19188241993679e+01,
         0.41078580492196e+00, -0.33465378172097e+00,  0.13870034777505e+04, -0.40663326195838e+03,  0.41727347159610e+02,  0.21932549434532e+01,
        -0.10320050009077e+01,  0.35882943516703e+00,  0.52511453726066e-02,  0.12838916450705e+02, -0.28642437219381e+01,  0.56912683664855e+00,
        -0.99962954584931e-01, -0.32632037778459e-02,  0.23320922576723e-03, -0.15334809857450e+00,  0.29072288239902e-01,  0.37534702741167e-03,
         0.17296691702411e-02, -0.38556050844504e-03, -0.35017712292608e-04, -0.14566393631492e-04,  0.56420857267269e-05,  0.41286150074605e-07,
        -0.20684671118824e-07,  0.16409393674725e-08]
Ic_s = [-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7]
Jc_s = [ 0,  1,  0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5]
nc_s = [ 0.90968501005365e+03,  0.24045667088420e+04, -0.59162326387130e+03,  0.54145404128074e+03, -0.27098308411192e+03,  0.97976525097926e+03,
        -0.46966772959435e+03,  0.14399274604723e+02, -0.19104204230429e+02,  0.53299167111971e+01, -0.21252975375934e+02, -0.31147334413760e+00,
         0.60334840894623e+00, -0.42764839702509e-01,  0.58185597255259e-02, -0.14597008284753e-01,  0.56631175631027e-02, -0.76155864584577e-04,
         0.22440342919332e-03, -0.12561095013413e-04,  0.63323132660934e-06, -0.20541989675375e-05,  0.36405370390082e-07, -0.29759897789215e-08,
         0.10136618529763e-07,  0.59925719692351e-11, -0.20677870105164e-10, -0.20874278181886e-10,  0.10162166825089e-09, -0.16429828281347e-09]
Ps_bs = 1.0     # [Mpa      ]
Ts_bs = 1.0     # [K        ]
ss_bsa = 2.0    # [kJ / kg K]
ss_bsb = 0.7853 # [kJ / kg K]
ss_bsc = 2.9251 # [kJ / kg K]

# Boundaries defining Region 2, and subregions 2a, 2b, 2c
def bnd2b2c(P):
    """ Boundary between region 2b and 2c"""
    pi = P / Ps_br
    eta = n[3] + ((pi - n[4]) / n[2])**0.5 
    return eta * hs_br

def idRegion_h(P, h):
    if P <= 4.0:
        return 1
    elif P <= 6.54670:
        return 2
    elif h >= bnd2b2c(P):
        return 2
    else:
        return 3

def idRegion_s(P, s):
    if P <= 4.0:
        return 1
    elif P <= 6.54670:
        return 2
    elif s >= 5.85:
        return 2
    else:
        return 3
Tbnd25 = 1073.15

#### dimensionless functions ####
def gamma(pi, tau):
    """ Dimensionless form for the specific Gibbs free energy"""
    sum = log(pi)
    for Ji, ni in zip(J0, n0):
        sum += ni * tau**Ji
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * (tau - 0.5)**Ji
    return sum

def gamma_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""
    sum = 1 / pi
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gamma_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""   
    sum = -1 / pi**2
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gamma_tau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * tau**(Ji - 1)
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gamma_tautau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * (Ji - 1) * tau**(Ji - 2)
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * (Ji - 1) * (tau - 0.5)**(Ji - 2)
    return sum

def gamma_pitau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and temperature (tau)"""
    sum = 0
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gammaR_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""
    sum = 0
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gammaR_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""   
    sum = 0
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gammaR_pitau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""
    sum = 0
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def theta2a_h(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2a)"""
    sum = 0
    for Ii, Ji, ni in zip(Ia_h, Ja_h, na_h):
        sum += ni * pi**Ii * (eta - 2.1)**Ji
    return sum

def theta2b_h(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2b)"""
    sum = 0
    for Ii, Ji, ni in zip(Ib_h, Jb_h, nb_h):
        sum += ni * (pi - 2)**Ii * (eta - 2.6)**Ji
    return sum

def theta2c_h(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2c)"""
    sum = 0
    for Ii, Ji, ni in zip(Ic_h, Jc_h, nc_h):
        sum += ni * (pi + 25)**Ii * (eta - 1.8)**Ji
    return sum

def theta_h(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy"""
    P = pi * Ps_bh
    h = eta * hs_bh
    region = idRegion_h(P, h)
    if region is 1:
        return theta2a_h(pi, eta)
    elif region is 2:
        return theta2b_h(pi, eta)
    else:
        return theta2c_h(pi, eta)

def theta2a_s(pi, sigma):
    """ Dimensionless form for the temperature 
        as a function of pressure and entropy (2a)"""
    sum = 0
    for Ii, Ji, ni in zip(Ia_s, Ja_s, na_s):
        sum += ni * pi**Ii * (sigma - 2)**Ji
    return sum

def theta2b_s(pi, sigma):
    """ Dimensionless form for the temperature 
        as a function of pressure and entropy (2b)"""
    sum = 0
    for Ii, Ji, ni in zip(Ib_s, Jb_s, nb_s):
        sum += ni * pi**Ii * (10 - sigma)**Ji
    return sum

def theta2c_s(pi, sigma):
    """ Dimensionless form for the temperature 
        as a function of pressure and entropy (2c)"""
    sum = 0
    for Ii, Ji, ni in zip(Ic_s, Jc_s, nc_s):
        sum += ni * pi**Ii * (2 - sigma)**Ji
    return sum

def theta_s(pi, sigma):
    """ Dimensionless form for the temperature 
        as a function of pressure and entropy"""
    P = pi * Ps_bs
    s = sigma * ss_bsa
    region = idRegion_s(P, s)
    if region is 1:
        return theta2a_s(pi, sigma)
    elif region is 2:
        return theta2b_s(pi, s / ss_bsb)
    else:
        return theta2c_s(pi, s / ss_bsc)

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

#### region 2 properties ####
def g(P, T):
    """ Specific Gibbs free energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return gamma(pi, tau) * R * T

def v(P, T):
    """ Specific volume [m^3 / kg]"""
    pi = P / Ps
    tau = Ts / T
    return pi * gamma_pi(pi, tau) * R * T / (P * 1e+03)

def u(P, T):
    """ Specific internal energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - pi * gamma_pi(pi, tau)) * R * T    

def s(P, T):
    """ Specific entropy [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - gamma(pi, tau)) * R

def h(P, T):
    """ Specific enthalpy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return tau * gamma_tau(pi, tau) * R * T

def cp(P, T):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return -tau**2 * gamma_tautau(pi, tau) * R

def cv(P, T):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return (-tau**2 * gamma_tautau(pi, tau) - (1 + pi * gammaR_pi(pi, tau) - tau * gammaR_pitau(pi, tau))**2 / (1 - pi**2 * gammaR_pipi(pi, tau))) * R

def w(P, T):
    """ Speed of sound [m / s]"""
    pi = P / Ps
    tau = Ts / T
    return (R * T * 1e+03 * (1 + 2 * pi * gammaR_pi(pi, tau) + pi**2 * gammaR_pi(pi, tau)**2) / ((1 - pi**2 * gammaR_pipi(pi, tau)) + (1 + pi * gammaR_pi(pi, tau) - tau * pi * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau))))**0.5

def av(P, T):
    """Isobaric cubic expansion coefficient [1 / K]"""
    pi = P / Ps
    tau = Ts / T
    return ((1 + pi * gammaR_pi(pi, tau) - tau * pi* gammaR_pitau(pi, tau)) / (1 + pi * gammaR_pi(pi, tau))) / T

def kT(P, T):
    """Isothermal compressibility [m^3 / kJ]"""
    pi = P / Ps
    tau = Ts / T
    return ((1 - pi**2 * gammaR_pipi(pi, tau)) / (1 + pi * gammaR_pi(pi, tau))) / (P * 1e+03)

#### region 2 property derivatives ####
def dgdP(P, T):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T)

def dvdP(P, T):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return -v(P, T) * kT(P, T)

def dudP(P, T):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T) * (P * 1e+03 * kT(P, T) - T * av(P, T))

def dsdP(P, T):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""
    return -v(P, T) * av(P, T)

def dhdP(P, T):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T) * (1 - T * av(P, T))

def dgdT(P, T):
    """ Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return -s(P, T)

def dvdT(P, T):
    """ Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""
    return v(P, T) * av(P, T)

def dudT(P, T):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return cp(P, T) - P * 1e+03 * v(P, T) * av(P, T)

def dsdT(P, T):
    """ Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant pressure"""
    return cp(P, T) / T

def dhdT(P, T):
    """ Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return cp(P, T)

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

#### region 2 properties ####
def g_h(P, h):
    """ Specific gibbs free energy [kJ / kg]"""
    return g(P, T_h(P, h))

def v_h(P, h):
    """ Specific volume [m^3 / kg]"""
    return v(P, T_h(P, h))

def u_h(P, h):
    """ Specific internal energy [kJ / kg]"""
    return u(P, T_h(P, h))   

def s_h(P, h):
    """ Specific entropy [kJ / kg K]"""
    return s(P, T_h(P, h))

def T_h(P, h):
    """ Temperature [K]"""
    pi = P / Ps_bh
    eta = h / hs_bh
    return theta_h(pi, eta) * Ts_bh

def cp_h(P, h):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    return cp(P, T_h(P, h))

def cv_h(P, h):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    return cv(P, T_h(P, h))

def w_h(P, h):
    """ Speed of sound [m / s]"""
    return w(P, T_h(P, h))

def av_h(P, h):
    """Isobaric cubic expansion coefficient [1 / K]"""
    return av(P, T_h(P, h))

def kT_h(P, h):
    """Isothermal compressibility [m^3 / kJ]"""
    return kT(P, T_h(P, h))

#### region 2 property derivatives ####
def dgdP_h(P, h):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dgdP(P, T) * dhdT(P, T) - dgdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dvdP_h(P, h):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dvdP(P, T) * dhdT(P, T) - dvdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dudP_h(P, h):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dudP(P, T) * dhdT(P, T) - dudT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dsdP_h(P, h):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dsdP(P, T) * dhdT(P, T) - dsdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dTdP_h(P, h):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return -dhdP(P, T) / dhdT(P, T)

def dgdh_h(P, h):
    """ Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dgdT(P, T) / dhdT(P, T)

def dvdh_h(P, h):
    """ Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dvdT(P, T) / dhdT(P, T)

def dudh_h(P, h):
    """ Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dudT(P, T) / dhdT(P, T)

def dsdh_h(P, h):
    """ Derivative of specific entropy [kJ kg / kg K kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dsdT(P, T) / dhdT(P, T)

def dTdh_h(P, h):
    """ Derivative of Temperature [K kg / kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return 1 / dhdT(P, T)

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################

#### region 2 properties ####
def g_s(P, s):
    """ Specific gibbs free energy [kJ / kg]"""
    return g(P, T_s(P, s))

def v_s(P, s):
    """ Specific volume [m^3 / kg]"""
    return v(P, T_s(P, s))

def u_s(P, s):
    """ Specific internal energy [kJ / kg]"""
    return u(P, T_s(P, s))   

def T_s(P, s):
    """ Temperature [K]"""
    pi = P / Ps_bs
    sigma = s / ss_bsa
    return theta_s(pi, sigma) * Ts_bs

def h_s(P, s):
    """ Specific enthalpy [kJ / kg]"""
    return h(P, T_s(P, s))

def cp_s(P, s):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    return cp(P, T_s(P, s))

def cv_s(P, s):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    return cv(P, T_s(P, s))

def w_s(P, s):
    """ Speed of sound [m / s]"""
    return w(P, T_s(P, s))

def av_s(P, s):
    """Isobaric cubic expansion coefficient [1 / K]"""
    return av(P, T_s(P, s))

def kT_s(P, s):
    """Isothermal compressibility [m^3 / kJ]"""
    return kT(P, T_s(P, s))

#### region 2 property derivatives ####
def dgdP_s(P, s):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dgdP(P, T) * dsdT(P, T) - dgdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dvdP_s(P, s):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dvdP(P, T) * dsdT(P, T) - dvdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dudP_s(P, s):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dudP(P, T) * dsdT(P, T) - dudT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dTdP_s(P, s):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return -dsdP(P, T) / dsdT(P, T)

def dhdP_s(P, s):
    """ Derivative of specific entropy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dhdP(P, T) * dsdT(P, T) - dhdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dgds_s(P, s):
    """ Derivative of specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dgdT(P, T) / dsdT(P, T)

def dvds_s(P, s):
    """ Derivative of specific volume [m^3 kg K/ kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dvdT(P, T) / dsdT(P, T)

def duds_s(P, s):
    """ Derivative of specific internal energy [kJ kg K/ kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dudT(P, T) / dsdT(P, T)

def dTds_s(P, s):
    """ Derivative of Temperature [K kg K / kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return 1 / dsdT(P, T)

def dhds_s(P, s):
    """ Derivative of specific enthalpy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dhdT(P, T) / dsdT(P, T)
