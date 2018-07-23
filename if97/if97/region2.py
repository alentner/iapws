from math import log

# Constants for region 2
J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
n0 = [-0.96927686500217e1,  0.10086655968018e2, -0.56087911283020e-2, 0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1,
      -0.43839511319450e1, -0.28408632460772,    0.21268463753307e-1]
Ir = [1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  5,  6,  6, 6, 
      7,  7,  7,  8,  8,  9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
Jr = [0,  1,  2,  3,  6,  1,  2,  4,  7, 36,  0,  1,  3,  6, 35,  1,  2,  3,  7,  3, 16, 35,
      0, 11, 25,  8, 36, 13,  4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
nr = [-0.0017731742473212999,  -0.017834862292357999,   -0.045996013696365003,   -0.057581259083432,     -0.050325278727930002,   -3.3032641670203e-05,
      -0.00018948987516315,    -0.0039392777243355001,  -0.043797295650572998,   -2.6674547914087001e-05, 2.0481737692308999e-08,  4.3870667284435001e-07, 
      -3.2277677238570002e-05, -0.0015033924542148,     -0.040668253562648998,   -7.8847309559367001e-10, 1.2790717852285001e-08,  4.8225372718507002e-07,
       2.2922076337661001e-06, -1.6714766451061001e-11, -0.0021171472321354998, -23.895741934103999,     -5.9059564324270004e-18, -1.2621808899101e-06, 
      -0.038946842435739003,    1.1256211360459e-11,    -8.2311340897998004,      1.9809712802088e-08,    1.0406965210174e-19,    -1.0234747095929e-13, 
      -1.0018179379511e-09,    -8.0882908646984998e-11,  0.10693031879409,       -0.33662250574170999,    8.9185845355420999e-25,  3.0629316876231997e-13, 
      -4.2002467698208001e-06, -5.9056029685639003e-26,  3.7826947613457002e-06, -1.2768608934681e-15,    7.3087610595061e-29,     5.5414715350778001e-17,
      -9.4369707241209998e-07]
n  = [0.90584275814723e3, -0.67955786399241, 0.12809002730136e-3, 0.26526571908428e4, 0.45257578905948e1]
Ia = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7]
Ja = [0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
na = [ 0.10898952318288e4,  0.84951654495535e3, -0.10781748091826e3,   0.33153654801263e2,  -0.74232016790248e1,   0.11765048724356e2,
       0.18445749355790e1, -0.41792700549624e1,  0.62478196935812e1,  -0.17344563108114e2,  -0.20058176862096e3,   0.27196065473796e3,
      -0.45511318285818e3,  0.30919688604755e4,  0.25226640357872e6,  -0.61707422868339e-2, -0.31078046629583,     0.11670873077107e2,
       0.12812798404046e9, -0.98554909623276e9,  0.28224546973002e10, -0.35948971410703e10,  0.17227349913197e10, -0.13551334240775e5,
       0.12848734664650e8,  0.13865724283226e1,  0.23598832556514e6,  -0.13105236545054e8,   0.73999835474766e4,  -0.55196697030060e6,
       0.37154085996233e7,  0.19127729239660e5, -0.41535164835634e6,  -0.62459855192507e2]
Ib = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9]
Jb = [0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40]
nb = [ 0.14895041079516e4,   0.74307798314034e3,  -0.97708318797837e2,   0.24742464705674e1,  -0.63281320016026,      0.11385952129658e1,
      -0.47811863648625,     0.85208123431544e-2,  0.93747147377932,     0.33593118604916e1,   0.33809355601454e1,    0.16844539671904,
       0.73875745236695,    -0.47128737436186,     0.15020273139707,    -0.21764114219750e-2, -0.21810755324761e-1,  -0.10829784403677,
      -0.46333324635812e-1,  0.71280351959551e-4,  0.11032831789999e-3,  0.18955248387902e-3,  0.30891541160537e-2,   0.13555504554949e-2,
       0.28640237477456e-6, -0.10779857357512e-4, -0.76462712454814e-4,  0.14052392818316e-4, -0.31083814331434e-4,  -0.10302738212103e-5,
       0.28217281635040e-6,  0.12704902271945e-5,  0.73803353468292e-7, -0.11030139238909e-7, -0.81456365207833e-13, -0.25180545682962e-10,
      -0.17565233969407e-17, 0.86934156344163e-14]
Ic = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6]
Jc = [0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22]
nc = [-0.32368398555242e13,  0.73263350902181e13, 0.35825089945447e12, -0.58340131851590e12, -0.10783068217470e11,   0.20825544563171e11,
       0.61074783564516e6,   0.85977722535580e6, -0.25745723604170e5,   0.31081088422714e5,   0.12082315865936e4,    0.48219755109255e3,
       0.37966001272486e1,  -0.10842984880077e2, -0.45364172676660e-1,  0.14559115658698e-12, 0.11261597407230e-11, -0.17804982240686e-10,
       0.12324579690832e-6, -0.11606921130984e-5, 0.27846367088554e-4, -0.59270038474176e-3,  0.12918582991878e-2]

# Non-dimenionalization for region 2
Ps = 1.0        #[Mpa]
Ts = 540.0      #[K]
hs = 1.0        #[kJ / kg]
R  = 0.461526   #[kJ / kg K]
Pb = 1.0        #[Mpa]
Tb = 1.0        #[K]
hb = 2000       #[kJ / kg]

# Boundaries defining subregions 2a, 2b, & 2c
def idRegion(P, h):
    if P <= 4.0:
        return 1
    elif P <= 6.54670:
        return 2
    elif h >= bnd2b2c(P):
        return 2
    else:
        return 3
def bnd2b2c(P):
    """ Boundary between region 2b and 2c"""
    pi = P / Ps
    eta = n[3] + ((pi - n[4]) / n[2])**0.5 
    return eta * hs
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
def theta2a(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2a)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ia, Ja, na):
        sum += ni * pib**Ii * (eta - 2.1)**Ji
    return sum
def theta2b(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2b)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += ni * (pib - 2)**Ii * (eta - 2.6)**Ji
    return sum
def theta2c(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy (2c)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ic, Jc, nc):
        sum += ni * (pib + 25)**Ii * (eta - 1.8)**Ji
    return sum
def theta(pi, eta):
    """ Dimensionless form for the temperature 
        as a function of pressure and enthalpy"""
    P = pi * Ps
    h = eta * hb
    region = idRegion(P, h)

    if region is 1:
        return theta2a(pi, eta)
    elif region is 2:
        return theta2b(pi, eta)
    else:
        return theta2c(pi, eta)

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

    return pi * gamma_pi(pi, tau) * R * T / (P * 10**6 / 1000)
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

    return (R * T * 1000 * (1 + 2 * pi * gammaR_pi(pi, tau) + pi**2 * gammaR_pi(pi, tau)**2) / ((1 - pi**2 * gammaR_pipi(pi, tau)) + (1 + pi * gammaR_pi(pi, tau) - tau * pi * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau))))**0.5
def a(P, T):
    """Isobaric cubic expansion coefficient [1 / K]"""
    pi = P / Ps
    tau = Ts / T

    return ((1 + pi * gammaR_pi(pi, tau) - tau * pi* gammaR_pitau(pi, tau)) / (1 + pi * gammaR_pi(pi, tau))) / T
def k(P, T):
    """Isothermal compressibility [kg / kJ]"""
    pi = P / Ps
    tau = Ts / T

    return ((1 - pi**2 * gammaR_pipi(pi, tau)) / (1 + pi * gammaR_pi(pi, tau))) / (P * 10**6 / 1000)

#### region 2 property derivatives ####
def dgdP(P, T):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""

    return v(P, T)
def dvdP(P, T):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""

    return -v(P, T) * k(P, T)
def dudP(P, T):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""

    return v(P, T) * ((P * 10**6 / 1000) * k(P, T) - T * a(P, T))
def dsdP(P, T):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""

    return -v(P, T) * a(P, T)
def dhdP(P, T):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""

    return v(P, T) * (1 - T * a(P, T))

def dgdT(P, T):
    """ Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""

    return -s(P, T)
def dvdT(P, T):
    """ Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""

    return v(P, T) * a(P, T)
def dudT(P, T):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""

    return cp(P, T) - (P * 10**6 / 1000) * v(P, T) * a(P, T)
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
    pi = P / Ps
    eta = h / hb

    return theta(pi, eta) * Tb
def cp_h(P, h):
    """ Specific isobaric heat capacity [kJ / kg K]"""

    return cp(P, T_h(P, h))
def cv_h(P, h):
    """ Specific isochoric heat capacity [kJ / kg K]"""

    return cv(P, T_h(P, h))
def w_h(P, h):
    """ Speed of sound [m / s]"""

    return w(P, T_h(P, h))
def a_h(P, h):
    """Isobaric cubic expansion coefficient [1 / K]"""

    return a(P, T_h(P, h))
def k_h(P, h):
    """Isothermal compressibility [kg / kJ]"""

    return k(P, T_h(P, h))

#### region 2 property derivatives ####
def dgdP_h(P, h):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant enthalpy"""
    T = T_h(P, h)

    return (dgdP(P, T) * dhdT(P, T) - dgdT(P, T) * dhdP(P, T)) / dhdT(P, T)
def dvdP_h(P, h):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant enthalpy"""
    T = T_h(P, h)

    return (dvdP(P, T) * dhdT(P, T) - dvdT(P, T) * dhdP(P, T)) / dhdT(P, T)
def dudP_h(P, h):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant enthalpy"""
    T = T_h(P, h)

    return (dudP(P, T) * dhdT(P, T) - dudT(P, T) * dhdP(P, T)) / dhdT(P, T)
def dsdP_h(P, h):
    """ Derivative of specific entropy [kJ m^3 / kg kJ]
    w.r.t pressure at constant enthalpy"""
    T = T_h(P, h)

    return (dsdP(P, T) * dhdT(P, T) - dsdT(P, T) * dhdP(P, T)) / dhdT(P, T)
def dTdP_h(P, h):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant enthalpy"""
    T = T_h(P, h)

    return -dhdP(P, T) / dhdT(P, T)

def dgdh_h(P, h):
    """ Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t enthalpy at constant pressure"""
    T = T_h(P, h)

    return dgdT(P, T) / dhdT(P, T)
def dvdh_h(P, h):
    """ Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t enthalpy at constant pressure"""
    T = T_h(P, h)

    return dvdT(P, T) / dhdT(P, T)
def dudh_h(P, h):
    """ Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t enthalpy at constant pressure"""
    T = T_h(P, h)

    return dvdT(P, T) / dhdT(P, T)
def dsdh_h(P, h):
    """ Derivative of specific entropy [kJ kg / kg kJ]
    w.r.t enthalpy at constant pressure"""
    T = T_h(P, h)

    return dvdT(P, T) / dhdT(P, T)
def dTdh_h(P, h):
    """ Derivative of Temperature [K kg / kJ]
    w.r.t enthalpy at constant pressure"""
    T = T_h(P, h)

    return 1 / dhdT(P, T)