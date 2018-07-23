# Constants for region 1
I = [ 0,  0, 0, 0, 0, 0, 0, 0,  1,  1,  1, 1, 1, 1,  2, 2, 2, 2,  2,  3, 3, 3,  4,  4,  4,  5,   8,  8,  21,  23,  29,  30,  31,  32]
J = [ -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
n = [ 0.14632971213167,     -0.84548187169114,    -0.37563603672040e1,    0.33855169168385e1,  -0.95791963387872,     0.15772038513228,    
     -0.16616417199501e-1,   0.81214629983568e-3,  0.28319080123804e-3,  -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1, 
     -0.21841717175414e-1,  -0.52838357969930e-4, -0.47184321073267e-3,  -0.30001780793026e-3,  0.47661393906987e-4, -0.44141845330846e-5, 
     -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5,  -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6, 
     -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,  -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19, 
      0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23, -0.93537087292458e-25]
Ib = [0, 0, 0, 0,  0,  0, 1, 1, 1, 1, 1,  1,  1,  2,  2,  3,  3,  4,  5,  6]
Jb = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
nb = [-0.23872489924521e3,    0.40421188637945e3,   0.11349746881718e3, -0.58457616048039e1, -0.15285482413140e-3,  -0.10866707695377e-5,
     -0.13391744872602e2,    0.43211039183559e2,  -0.54010067170506e2,  0.30535892203916e2, -0.65964749423638e1,    0.93965400878363e-2,
      0.11573647505340e-6,  -0.25858641282073e-4, -0.40644363084799e-8, 0.66456186191635e-7, 0.80670734103027e-10, -0.93477771213947e-12,
      0.58265442020601e-14, -0.15020185953503e-16]

# Non-dimenionalization for region 1
Ps = 16.53      #[Mpa]
Ts = 1386.0     #[K]
R  = 0.461526   #[kJ / kg K]
Pb = 1.0        #[Mpa]
Tb = 1.0        #[K]
hb = 2500       #[kJ / kg]

# Boundaries defining region 1
Tbnd01 = 273.15     #[K]
Tbnd13 = 623.15     #[K]
Pbnd0  = 0.0        #[MPa]
Pbnd1  = 100        #[MPa]

#### dimensionless functions ####
def gamma(pi, tau):
    """ Dimensionless form for the specific Gibbs free energy"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * (tau - 1.222)**Ji
    return sum
def gamma_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (tau - 1.222)**Ji
    return sum
def gamma_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""   
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (tau - 1.222)**Ji
    return sum
def gamma_tau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (tau - 1.222)**(Ji - 1)
    return sum
def gamma_tautau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * (tau - 1.222)**(Ji - 2)
    return sum
def gamma_pitau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and temperature (tau)"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (tau - 1.222)**(Ji - 1)
    return sum
def theta(pi, eta):
    """ Dimensionless form for the temperature as a function of pressure and enthalpy"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += ni * pib**Ii * (eta + 1.0)**Ji
    return sum

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

#### region 1 properties ####
def g(P, T):
    """ Specific gibbs free energy [kJ / kg]"""
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

    return (-tau**2 * gamma_tautau(pi, tau) + (gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / gamma_pipi(pi, tau)) * R
def w(P, T):
    """ Speed of sound [m / s]"""
    pi = P / Ps
    tau = Ts / T

    return (gamma_pi(pi, tau)**2 * R * T * 1000 / ((gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau)) - gamma_pipi(pi, tau)))**0.5
def a(P, T):
    """Isobaric cubic expansion coefficient [1 / K]"""
    pi = P / Ps
    tau = Ts / T

    return (1 - tau * gamma_pitau(pi, tau) / gamma_pi(pi, tau)) / T
def k(P, T):
    """Isothermal compressibility [kg / kJ]"""
    pi = P / Ps
    tau = Ts / T

    return -(pi * gamma_pipi(pi, tau) / gamma_pi(pi, tau)) / (P * 10**6 / 1000)

#### region 1 property derivatives ####
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

#### region 1 properties ####
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

#### region 1 property derivatives ####
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