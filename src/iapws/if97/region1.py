###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# constants and non-dimenionalization;
# Region 1, forwards equations for (P, T)
I = [  0,  0, 0, 0, 0, 0, 0, 0,  1,  1,  1, 1, 1, 1,  2, 2, 2, 2,  2,  3, 3, 3,  4,  4,  4,  5,   8,  8,  21,  23,  29,  30,  31,  32]
J = [ -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
n = [ 0.14632971213167e+00, -0.84548187169114e+00, -0.37563603672040e+01,  0.33855169168385e+01, -0.95791963387872e+00,  0.15772038513228e+00,    
     -0.16616417199501e-01,  0.81214629983568e-03,  0.28319080123804e-03, -0.60706301565874e-03, -0.18990068218419e-01, -0.32529748770505e-01, 
     -0.21841717175414e-01, -0.52838357969930e-04, -0.47184321073267e-03, -0.30001780793026e-03,  0.47661393906987e-04, -0.44141845330846e-05, 
     -0.72694996297594e-15, -0.31679644845054e-04, -0.28270797985312e-05, -0.85205128120103e-09, -0.22425281908000e-05, -0.65171222895601e-06, 
     -0.14341729937924e-12, -0.40516996860117e-06, -0.12734301741641e-08, -0.17424871230634e-09, -0.68762131295531e-18,  0.14478307828521e-19, 
      0.26335781662795e-22, -0.11947622640071e-22,  0.18228094581404e-23, -0.93537087292458e-25]
Ps = 16.53    # [Mpa      ]
Ts = 1386.0   # [K        ]
R  = 0.461526 # [kJ / kg K]

# constants and non-dimenionalization;
# Region 1, backwards equations for (P, h)
I_bh = [0, 0, 0, 0,  0,  0, 1, 1, 1, 1, 1,  1,  1,  2,  2,  3,  3,  4,  5,  6]
J_bh = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
n_bh = [-0.23872489924521e+03,  0.40421188637945e+03,  0.11349746881718e+03, -0.58457616048039e+01, -0.15285482413140e-03, -0.10866707695377e-05,
        -0.13391744872602e+02,  0.43211039183559e+02, -0.54010067170506e+02,  0.30535892203916e+02, -0.65964749423638e+01,  0.93965400878363e-02,
         0.11573647505340e-06, -0.25858641282073e-04, -0.40644363084799e-08,  0.66456186191635e-07,  0.80670734103027e-10, -0.93477771213947e-12,
         0.58265442020601e-14, -0.15020185953503e-16]
Ps_bh = 1.0  # [Mpa    ]
Ts_bh = 1.0  # [K      ]
hs_bh = 2500 # [kJ / kg]

# constants and non-dimenionalization;
# Region 1, backwards equations for (P, s)
I_bs = [0, 0, 0, 0,  0,  0, 1, 1, 1, 1,  1,  1, 2, 2, 2, 2,  2,  3,  3,  4]
J_bs = [0, 1, 2, 3, 11, 31, 0, 1, 2, 3, 12, 31, 0, 1, 2, 9, 31, 10, 32, 32]
n_bs = [ 0.17478268058307e+03,  0.34806930892873e+02,  0.65292584978455e+01,  0.33039981775489e+00, -0.19281382923196e-06, -0.24909197244573e-22,
        -0.26107636489332e+00,  0.22592965981586e+00, -0.64256463395226e-01,  0.78876289270526e-02,  0.35672110607366e-09,  0.17332496994895e-23,
         0.56608900654837e-03, -0.32635483139717e-03,  0.44778286690632e-04, -0.51322156908507e-09, -0.42522657042207e-25,  0.26400441360689e-12,
         0.78124600459723e-28, -0.30732199903668e-30]
Ps_bs = 1.0 # [Mpa      ]
Ts_bs = 1.0 # [K        ]
ss_bs = 1.0 # [kJ / kg K]

# boundaries defining Region 1
Tbnd01 = 273.16 # [K  ]
Tbnd13 = 623.15 # [K  ]
Pbnd0  = 1.0e-6 # [MPa]
Pbnd1  = 100    # [MPa]

#### dimensionless functions ####
def gamma(pi, tau):
    """ Dimensionless form for the specific Gibbs free energy"""
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * (tau - 1.222)**Ji
    return sum

def gamma_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (tau - 1.222)**Ji
    return sum

def gamma_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""   
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (tau - 1.222)**Ji
    return sum

def gamma_tau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (tau - 1.222)**(Ji - 1)
    return sum

def gamma_tautau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * (tau - 1.222)**(Ji - 2)
    return sum

def gamma_pitau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and temperature (tau)"""
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (tau - 1.222)**(Ji - 1)
    return sum

def theta_T(pi, eta):
    """ Dimensionless form for the temperature as a function of pressure and enthalpy"""
    sum = 0.0
    for Ii, Ji, ni in zip(I_bh, J_bh, n_bh):
        sum += ni * pi**Ii * (eta + 1.0)**Ji
    return sum

def theta_s(pi, sigma):
    """ Dimensionless form for the temperature as a function of pressure and entropy"""
    sum = 0.0
    for Ii, Ji, ni in zip(I_bs, J_bs, n_bs):
        sum += ni * pi**Ii * (sigma + 2.0)**Ji
    return sum

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

#### region 1 properties ####
def g(P, T):
    """Specific gibbs free energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return gamma(pi, tau) * R * T

def v(P, T):
    """Specific volume [m^3 / kg]"""
    pi = P / Ps
    tau = Ts / T
    return pi * gamma_pi(pi, tau) * R * T / (P * 10**6 / 1000)

def u(P, T):
    """Specific internal energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - pi * gamma_pi(pi, tau)) * R * T    

def s(P, T):
    """Specific entropy [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - gamma(pi, tau)) * R

def h(P, T):
    """Specific enthalpy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T
    return tau * gamma_tau(pi, tau) * R * T

def cp(P, T):
    """Specific isobaric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return -tau**2 * gamma_tautau(pi, tau) * R

def cv(P, T):
    """Specific isochoric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T
    return (-tau**2 * gamma_tautau(pi, tau) + (gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / gamma_pipi(pi, tau)) * R

def w(P, T):
    """Speed of sound [m / s]"""
    pi = P / Ps
    tau = Ts / T
    return (gamma_pi(pi, tau)**2 * R * T * 1e+03 / ((gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau)) - gamma_pipi(pi, tau)))**0.5

def av(P, T):
    """Isobaric cubic expansion coefficient [1 / K]"""
    pi = P / Ps
    tau = Ts / T
    return (1 - tau * gamma_pitau(pi, tau) / gamma_pi(pi, tau)) / T

def kT(P, T):
    """Isothermal compressibility [m^3 / kJ]"""
    pi = P / Ps
    tau = Ts / T
    return -(pi * gamma_pipi(pi, tau) / gamma_pi(pi, tau)) / (P * 1e+03)

#### region 1 property derivatives ####
def dgdP(P, T):
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T)

def dvdP(P, T):
    """Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return -v(P, T) * kT(P, T)

def dudP(P, T):
    """Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T) * (P * 1e+03 * kT(P, T) - T * av(P, T))

def dsdP(P, T):
    """Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""
    return -v(P, T) * av(P, T)

def dhdP(P, T):
    """Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    return v(P, T) * (1 - T * av(P, T))

def dgdT(P, T):
    """Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return -s(P, T)

def dvdT(P, T):
    """Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""
    return v(P, T) * av(P, T)

def dudT(P, T):
    """Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return cp(P, T) - P * 1e+03 * v(P, T) * av(P, T)

def dsdT(P, T):
    """Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant pressure"""
    return cp(P, T) / T

def dhdT(P, T):
    """Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    return cp(P, T)

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

#### region 1 properties ####
def g_h(P, h):
    """Specific gibbs free energy [kJ / kg]"""
    return g(P, T_h(P, h))

def v_h(P, h):
    """Specific volume [m^3 / kg]"""
    return v(P, T_h(P, h))

def u_h(P, h):
    """Specific internal energy [kJ / kg]"""
    return u(P, T_h(P, h))   

def s_h(P, h):
    """Specific entropy [kJ / kg K]"""
    return s(P, T_h(P, h))

def T_h(P, h):
    """Temperature [K]"""
    pi = P / Ps_bh
    eta = h / hs_bh
    return theta_T(pi, eta) * Ts_bh

def cp_h(P, h):
    """Specific isobaric heat capacity [kJ / kg K]"""
    return cp(P, T_h(P, h))

def cv_h(P, h):
    """Specific isochoric heat capacity [kJ / kg K]"""
    return cv(P, T_h(P, h))

def w_h(P, h):
    """Speed of sound [m / s]"""
    return w(P, T_h(P, h))

def av_h(P, h):
    """Isobaric cubic expansion coefficient [1 / K]"""
    return av(P, T_h(P, h))

def kT_h(P, h):
    """Isothermal compressibility [m^3 / kJ]"""
    return kT(P, T_h(P, h))

#### region 1 property derivatives ####
def dgdP_h(P, h):
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dgdP(P, T) * dhdT(P, T) - dgdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dvdP_h(P, h):
    """Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dvdP(P, T) * dhdT(P, T) - dvdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dudP_h(P, h):
    """Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dudP(P, T) * dhdT(P, T) - dudT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dsdP_h(P, h):
    """Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return (dsdP(P, T) * dhdT(P, T) - dsdT(P, T) * dhdP(P, T)) / dhdT(P, T)

def dTdP_h(P, h):
    """Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific enthalpy"""
    T = T_h(P, h)
    return -dhdP(P, T) / dhdT(P, T)

def dgdh_h(P, h):
    """Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dgdT(P, T) / dhdT(P, T)

def dvdh_h(P, h):
    """Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dvdT(P, T) / dhdT(P, T)

def dudh_h(P, h):
    """Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dudT(P, T) / dhdT(P, T)

def dsdh_h(P, h):
    """Derivative of specific entropy [kJ kg / kg K kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return dsdT(P, T) / dhdT(P, T)

def dTdh_h(P, h):
    """Derivative of Temperature [K kg / kJ]
    w.r.t specific enthalpy at constant pressure"""
    T = T_h(P, h)
    return 1 / dhdT(P, T)

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################

#### region 1 properties ####
def g_s(P, s):
    """Specific gibbs free energy [kJ / kg]"""
    return g(P, T_s(P, s))

def v_s(P, s):
    """Specific volume [m^3 / kg]"""
    return v(P, T_s(P, s))

def u_s(P, s):
    """Specific internal energy [kJ / kg]"""
    return u(P, T_s(P, s))   

def T_s(P, s):
    """Temperature [K]"""
    pi = P / Ps_bs
    sigma = s / ss_bs
    return theta_s(pi, sigma) * Ts_bs

def h_s(P, s):
    """Specific enthalpy [kJ / kg]"""
    return h(P, T_s(P, s))

def cp_s(P, s):
    """Specific isobaric heat capacity [kJ / kg K]"""
    return cp(P, T_s(P, s))

def cv_s(P, s):
    """Specific isochoric heat capacity [kJ / kg K]"""
    return cv(P, T_s(P, s))

def w_s(P, s):
    """Speed of sound [m / s]"""
    return w(P, T_s(P, s))

def av_s(P, s):
    """Isobaric cubic expansion coefficient [1 / K]"""
    return av(P, T_s(P, s))

def kT_s(P, s):
    """Isothermal compressibility [m^3 / kJ]"""
    return kT(P, T_s(P, s))

#### region 1 property derivatives ####
def dgdP_s(P, s):
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dgdP(P, T) * dsdT(P, T) - dgdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dvdP_s(P, s):
    """Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dvdP(P, T) * dsdT(P, T) - dvdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dudP_s(P, s):
    """Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dudP(P, T) * dsdT(P, T) - dudT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dTdP_s(P, s):
    """Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return -dsdP(P, T) / dsdT(P, T)

def dhdP_s(P, s):
    """Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    T = T_s(P, s)
    return (dhdP(P, T) * dsdT(P, T) - dhdT(P, T) * dsdP(P, T)) / dsdT(P, T)

def dgds_s(P, s):
    """Derivative of specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dgdT(P, T) / dsdT(P, T)

def dvds_s(P, s):
    """Derivative of specific volume [m^3 kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dvdT(P, T) / dsdT(P, T)

def duds_s(P, s):
    """Derivative of specific internal energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dudT(P, T) / dsdT(P, T)

def dTds_s(P, s):
    """Derivative of Temperature [K kg K/ kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return 1 / dsdT(P, T)

def dhds_s(P, s):
    """Derivative of specific enthalpy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    T = T_s(P, s)
    return dhdT(P, T) / dsdT(P, T)
