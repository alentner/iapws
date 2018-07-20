#####################################################
###    Extension for p, h as primary variables    ###
#####################################################

from if97 import region2
from math import log

# Constants for region 2
J0 = region2.J0
n0 = region2.n0
Ir = region2.Ir
Jr = region2.Jr
nr = region2.nr
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
Ps = region2.Ps #[Mpa]
Ts = region2.Ts #[K]
hs = 1.0        #[kJ / kg]
R  = region2.R  #[kJ / kg K]
Pb = 1.0        #[Mpa]
Tb = 1.0        #[K]
hb = 2000       #[kJ / kg]

def bnd2b2c(P):
    """ Boundary between region 2b and 2c"""
    pi = P / Ps
    eta = n[3] + ((pi - n[4]) / n[2])**0.5 
    return eta * hs

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

    if P <= 4.0:
        return theta2a(pi, eta)
    elif h > bnd2b2c(P):
        return theta2b(pi, eta)
    else:
        return theta2c(pi, eta)

def theta2a_pi(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2a)
        as a function of pressure and enthalpy w.r.t. pressure (pi)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ia, Ja, na):
        sum += (Ps / Pb) * ni * Ii * pib**(Ii - 1) * (eta - 2.1)**Ji
    return sum

def theta2b_pi(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2b)
        as a function of pressure and enthalpy w.r.t. pressure (pi)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += (Ps / Pb) * ni * Ii * (pib - 2)**(Ii - 1) * (eta - 2.6)**Ji
    return sum

def theta2c_pi(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2c)
        as a function of pressure and enthalpy w.r.t. pressure (pi)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ic, Jc, nc):
        sum += (Ps / Pb) * ni * Ii * (pib + 25)**(Ii - 1) * (eta - 1.8)**Ji
    return sum

def theta_pi(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature
        as a function of pressure and enthalpy w.r.t. pressure (pi)"""
    P = pi * Ps
    h = eta * hb

    if P <= 4.0:
        return theta2a_pi(pi, eta)
    elif h > bnd2b2c(P):
        return theta2b_pi(pi, eta)
    else:
        return theta2c_pi(pi, eta)

def theta2a_eta(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2a)
        as a function of pressure and enthalpy w.r.t. enthalpy (eta)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ia, Ja, na):
        sum += ni * pib**Ii * Ji * (eta - 2.1)**(Ji - 1)
    return sum

def theta2b_eta(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2b)
        as a function of pressure and enthalpy w.r.t. enthalpy (eta)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += ni * (pib - 2)**Ii * Ji * (eta - 2.6)**(Ji - 1)
    return sum

def theta2c_eta(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature (2c)
        as a function of pressure and enthalpy w.r.t. enthalpy (eta)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ic, Jc, nc):
        sum += ni * (pib + 25)**Ii * Ji * (eta - 1.8)**(Ji - 1)
    return sum

def theta_eta(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature
        as a function of pressure and enthalpy w.r.t. enthalpy (eta)"""
    P = pi * Ps
    h = eta * hb

    if P <= 4.0:
        return theta2a_eta(pi, eta)
    elif h > bnd2b2c(P):
        return theta2b_eta(pi, eta)
    else:
        return theta2c_eta(pi, eta)

def gamma(pi, eta):
    """ --- Not Tested ---
        Dimensionless form for the specific Gibbs free energy;
        using backward relation for temperature"""
    sum = log(pi)
    for Ji, ni in zip(J0, n0):
        sum += ni * (Ts / (Tb * theta(pi, eta)))**Ji

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * ((Ts / (Tb * theta(pi, eta))) - 0.5)**Ji
    return sum

def gamma_pi(pi, eta):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 1 / pi

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * ((Ts / (Tb * theta(pi, eta))) - 0.5)**Ji
    return sum

def gamma_pi_pi(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = -1 / pi**2

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += (ni * Ii * (Ii - 1) * pi**(Ii - 2) * ((Ts / (Tb * theta(pi, eta))) - 0.5)**Ji) + \
              (-ni * Ii * pi**(Ii - 1) * Ji * ((Ts * theta_pi(pi, eta)) / (Tb * theta(pi, eta)**2)) * ((Ts / (Tb * theta(pi, eta))) - 0.5)**(Ji - 1))
    return sum

def gamma_pi_eta(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and enthalpy (eta); using backward 
        relation for temperature"""
    sum = 1 / pi

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += -ni * Ii * pi**(Ii - 1) * Ji * ((Ts * theta_eta(pi, eta)) / (Tb * theta(pi, eta)**2)) * ((Ts / (Tb * theta(pi, eta))) - 0.5)**(Ji - 1)
    return sum

def gamma_tau(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau); using backward relation for 
        temperature"""
    if alpha == 1:
        theta = theta2a(pi, eta)
    elif alpha == 2:
        theta = theta2b(pi, eta)
    else:
        theta = theta2c(pi, eta)

    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * (Ts / (Tb * theta(pi, eta)))**(Ji - 1)

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * ((Ts / (Tb * theta(pi, eta))) - 0.5)**(Ji - 1)
    return sum

def gamma_tau_pi(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau) and pressure (pi); using backward 
        relation for temperature"""
    sum = 0 
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += (ni * Ii * pi**(Ii - 1) * Ji * ((Ts / (Tb * theta(pi, eta))) - 0.5)**(Ji - 1)) + \
              (-ni * pi**Ii * Ji * (Ji - 1) * ((Ts * theta_pi(pi, eta)) / (Tb * theta(pi, eta)**2)) * ((Ts / (Tb * theta(pi, eta))) - 0.5)**(Ji - 2))
    return sum

#### Insert the remaining backward relavent relations ####

def T(P, h):
    """ Temperature [K]"""
    pi = P / Ps
    eta = h / hb

    return theta(pi, eta) * Tb

def v(P, h):
    """ Specific volume [m^3 / kg]"""
    pi = P / Ps
    eta = h / hb

    return theta(pi, eta) * gamma_pi(pi, eta) * R * Tb / (Ps * 10**6 / 1000)

def dvdp(P, h):
    """ --- Not Tested ---
        Derivative of specific volume [m^3 m^3 / kg kJ]
        w.r.t. pressure"""
    pi = P / Ps
    eta = h / hb

    return (theta_pi(pi, eta) * gamma_pi(pi, eta) + theta(pi, eta) * gamma_pi_pi(pi, eta)) * R * Tb / (Ps * 10**6 / 1000)**2

def dvdh(P, h):
    """ --- Not Tested ---
        Derivative of specific volume [m^3 kg / kg kJ]
        w.r.t. enthalpy"""
    pi = P / Ps
    eta = h / hb

    return (theta_eta(pi, eta) * gamma_pi(pi, eta) + theta(pi, eta) * gamma_pi_eta(pi, eta)) * R * Tb / (hb * Ps * 10**6 / 1000)