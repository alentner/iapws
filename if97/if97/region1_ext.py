#####################################################
###    Extension for p, h as primary variables    ###
#####################################################

from if97 import region1

# Constants for region 1
I = region1.I
J = region1.J
n = region1.n
Ib = [0, 0, 0, 0,  0,  0, 1, 1, 1, 1, 1,  1,  1,  2,  2,  3,  3,  4,  5,  6]
Jb = [0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
nb = [-0.23872489924521e3,    0.40421188637945e3,   0.11349746881718e3, -0.58457616048039e1, -0.15285482413140e-3,  -0.10866707695377e-5,
     -0.13391744872602e2,    0.43211039183559e2,  -0.54010067170506e2,  0.30535892203916e2, -0.65964749423638e1,    0.93965400878363e-2,
      0.11573647505340e-6,  -0.25858641282073e-4, -0.40644363084799e-8, 0.66456186191635e-7, 0.80670734103027e-10, -0.93477771213947e-12,
      0.58265442020601e-14, -0.15020185953503e-16]

# Non-dimenionalization for region 1
Ps = region1.Ps #[Mpa]
Ts = region1.Ts #[K]
R  = region1.R  #[kJ / kg K]
Pb = 1.0        #[Mpa]
Tb = 1.0        #[K]
hb = 2500       #[kJ / kg]

def theta(pi, eta):
    """ Dimensionless form for the temperature as a function of pressure and enthalpy"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += ni * pib**Ii * (eta + 1.0)**Ji
    return sum

def theta_pi(pi, eta):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature 
        as a function of pressure and enthalpy w.r.t. pressure (pi)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += (Ps / Pb) * ni * Ii * pib**(Ii - 1) * (eta + 1.0)**Ji
    return sum

def theta_eta(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature 
        as a function of pressure and enthalpy w.r.t. enthalpy (eta)"""
    pib = pi * Ps / Pb
    sum = 0
    for Ii, Ji, ni in zip(Ib, Jb, nb):
        sum += ni * pib**Ii * Ji * (eta + 1.0)**(Ji - 1)
    return sum

def gammaB(pi, eta):
    """ --- Not Tested ---
        Dimensionless form for the specific Gibbs free energy;
        using backward relation for temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * (Ts / (Tb * theta(pi, eta)) - 1.222)**Ji
    return sum

def gamma_pi(pi, eta):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (Ts / (Tb * theta(pi, eta)) - 1.222)**Ji
    return sum

def gamma_pi_pi(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += (ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (Ts / (Tb * theta(pi, eta)) - 1.222)**Ji) + \
               (ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * ((Ts * theta_pi(pi, eta)) / (Tb * theta(pi, eta)**2)) * (Ts / (Tb * theta(pi, eta)) - 1.222)**(Ji - 1))
    return sum

def gamma_pi_eta(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and enthalpy (eta); using backward 
        relation for temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * ((Ts * theta_eta(pi, eta)) / (Tb * theta(pi, eta)**2)) * (Ts / (Tb * theta(pi, eta)) - 1.222)**(Ji - 1)
    return sum

def gamma_tau(pi, eta):
    """ --- Not Tested ---
        Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ts / (Tb * theta(pi, eta)) - 1.222)**(Ji - 1)
    return sum

def gamma_tau_pi(pi, eta):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau) and pressure (pi); using backward 
        relation for temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += (-ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (Ts / (Tb * theta(pi, eta)) - 1.222)**(Ji - 1)) + \
               (-ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * ((Ts * theta_pi(pi, eta)) / (Tb * theta(pi, eta)**2)) * (Ts / (Tb * theta(pi, eta)) - 1.222)**(Ji - 2))
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