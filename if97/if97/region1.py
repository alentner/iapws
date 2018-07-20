# Constants for region 1
I = [ 0,  0, 0, 0, 0, 0, 0, 0,  1,  1,  1, 1, 1, 1,  2, 2, 2, 2,  2,  3, 3, 3,  4,  4,  4,  5,   8,  8,  21,  23,  29,  30,  31,  32]
J = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
n = [ 0.14632971213167,     -0.84548187169114,    -0.37563603672040e1,    0.33855169168385e1,  -0.95791963387872,     0.15772038513228,    
     -0.16616417199501e-1,   0.81214629983568e-3,  0.28319080123804e-3,  -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1, 
     -0.21841717175414e-1,  -0.52838357969930e-4, -0.47184321073267e-3,  -0.30001780793026e-3,  0.47661393906987e-4, -0.44141845330846e-5, 
     -0.72694996297594e-15, -0.31679644845054e-4, -0.28270797985312e-5,  -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6, 
     -0.14341729937924e-12, -0.40516996860117e-6, -0.12734301741641e-8,  -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19, 
      0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23, -0.93537087292458e-25]

# Non-dimenionalization for region 1
Ps = 16.53      #[Mpa]
Ts = 1386.0     #[K]
R  = 0.461526   #[kJ / kg K]

# Boundaries defining region 1
Tbnd13 = 623.15     #[K]
Pbnd1 = 100         #[MPa]

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

    return (-tau**2 * gamma_tautau(pi, tau) + (gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / gamma_pipi(pi, tau)) * R

def w(P, T):
    """ Speed of sound [m / s]"""
    pi = P / Ps
    tau = Ts / T

    return (gamma_pi(pi, tau)**2 * R * T * 1000 / ((gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau)) - gamma_pipi(pi, tau)))**0.5