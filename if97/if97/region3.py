from math import log

# Constants for region 3
I = [ 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,  2,  2,  2, 2, 3, 3, 3, 3, 3, 4,
      4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
J = [ 0, 1, 2, 7, 10, 12, 23,  2, 6, 15, 17, 0,  2, 6,  7, 22, 26,  0, 2, 4, 16,
     26, 0, 2, 4, 26,  1,  3, 26, 0,  2, 26, 2, 26, 2, 26,  0,  1, 26]
n = [-0.15732845290239e2,   0.20944396974307e2,  -0.76867707878716e1,   0.26185947787954e1, -0.28080781148620e1,   0.12053369696517e1,
     -0.84566812812502e-2, -0.12654315477714e1,  -0.11524407806681e1,   0.88521043984318,   -0.64207765181607,     0.38493460186671,
     -0.85214708824206,     0.48972281541877e1,  -0.30502617256965e1,   0.39420536879154e-1, 0.12558408424308,    -0.27999329698710,
      0.13899799569460e1,  -0.20189915023570e1,  -0.82147637173963e-2, -0.47596035734923,    0.43984074473500e-1, -0.44476435428739,
      0.90572070719733,     0.70522450087967,     0.10770512626332,    -0.32913623258954,   -0.50871062041158,    -0.22175400873096e-1, 
      0.94260751665092e-1,  0.16436278447961,    -0.13503372241348e-1, -0.14834345352472e-1, 0.57922953628084e-3,  0.32308904703711e-2,
      0.80964802996215e-4, -0.16557679795037e-3, -0.44923899061815e-4]
n0 = 0.10658070028513e1
n23 = [ 0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2]

# Non-dimenionalization for region 3
Ps = 22.064     #[Mpa]
rhos = 322      #[kg / m^3]
Ts = 647.096    #[K]
R  = 0.461526   #[kJ / kg K]
R  = 0.461526   #[kJ / kg K]
Pb = 1.00       #[Mpa]
Tb = 1.00       #[K]

# Boundaries defining region 3
def bnd23P(T):
    theta = T / Tb

    return Pb * (n23[0] + n23[1] * theta + n23[2] * theta**2)
def bnd23T(P):
    pi = P / Pb

    return Tb * (n23[3] + ((pi - n23[4]) / n23[2])**0.5)

#### dimensionless functions ####
def phi(delta, tau):
    """ Dimensionless form for the Specific Helmholtz free energy"""
    sum = n0 * log(delta)

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * delta**Ii * tau**Ji
    return sum
def phi_delta(delta, tau):
    """ Derivative of dimensionless form for the 
        Specific Helmholtz free energy w.r.t. dimensionless 
        density (delta)"""
    sum = n0 / delta

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * delta**(Ii - 1) * tau**Ji
    return sum
def phi_deltadelta(delta, tau):
    """ Derivative (second) of dimensionless form for the 
        Specific Helmholtz free energy w.r.t. dimensionless 
        density (delta)"""
    sum = -n0 / delta**2

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * (Ii - 1) * delta**(Ii - 2) * tau**Ji
    return sum
def phi_tau(delta, tau):
    """ Derivative of dimensionless form for the 
        Specific Helmholtz free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * delta**Ii * Ji * tau**(Ji - 1)
    return sum
def phi_tautau(delta, tau):
    """ Derivative (second) of dimensionless form for the 
        Specific Helmholtz free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * delta**Ii * Ji * (Ji - 1) * tau**(Ji - 2)
    return sum
def phi_deltatau(delta, tau):
    """ Derivative (second) of dimensionless form for the 
        Specific Helmholtz free energy w.r.t. dimensionless 
        density (delta) and temperature (tau)"""
    sum = 0

    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * delta**(Ii - 1) * Ji * tau**(Ji - 1)
    return sum

###########################################################
#####     Specific Volume-Temperature Formulation     #####
###########################################################

#### region 3 properties ####
def f(nu, T):
    """ Specific Helmholtz free energy [kJ / kg]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return phi(delta, tau) * R * T
def P(nu, T):
    """ Pressure [Mpa]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return delta * phi_delta(delta, tau) * R * T * rho / (10**6 / 1000)
def u(nu, T):
    """ Specific internal energy [kJ / kg]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return tau * phi_tau(delta, tau) * R * T    
def s(nu, T):
    """ Specific entropy [kJ / kg K]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (tau * phi_tau(delta, tau) - phi(delta, tau)) * R
def h(nu, T):
    """ Specific enthalpy [kJ / kg]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (tau * phi_tau(delta, tau) + delta * phi_delta(delta, tau)) * R * T
def cp(nu, T):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (-tau**2 * phi_tautau(delta, tau) + (delta * phi_delta(delta, tau) - delta * tau * phi_deltatau(delta, tau))**2 / (2 * delta * phi_delta(delta, tau) + delta**2 * phi_deltadelta(delta, tau))) * R
def cv(nu, T):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (-tau**2 * phi_tautau(delta, tau)) * R
def w(nu, T):
    """ Speed of sound [m / s]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return ((2 * delta * phi_delta(delta, tau) + delta**2 * phi_deltadelta(delta, tau) - (delta * phi_delta(delta, tau) - delta * tau * phi_deltatau(delta, tau))**2 / (tau**2 * phi_tautau(delta, tau))) * R * T * 1000)**0.5
def a(nu, T):
    """Relative pressure coefficient [1 / K]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (1 - tau * phi_deltatau(delta, tau) / phi_delta(delta, tau)) / T
def b(nu, T):
    """Isothermal stress coefficient [kg / m^3]"""
    rho = 1 / nu
    delta = rho / rhos
    tau = Ts / T

    return (2 + delta * phi_deltadelta(delta, tau) / phi_delta(delta, tau)) * rho

#### region 3 property derivatives ####
def dfdv(nu, T):
    """ Derivative of specific helmholtz free energy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""

    return -P(nu, T) * (10**6 / 1000)
def dPdv(nu, T):
    """ Derivative of pressure [kJ kg / m^3 m^3]
    w.r.t specific volume at constant temperature"""

    return -P(nu, T) * b(nu, T)
def dudv(nu, T):
    """ Derivative of specific internal energy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""

    return P(nu, T) * (10**6 / 1000) * (T * a(nu, T) - 1)
def dsdv(nu, T):
    """ Derivative of specific entropy [kJ kg / kg K m^3]
    w.r.t specific volume at constant temperature"""

    return P(nu, T) * (10**6 / 1000) * a(nu, T)
def dhdv(nu, T):
    """ Derivative of specific enthalpy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""

    return P(nu, T) * (10**6 / 1000) * (T * a(nu, T) - nu * b(nu, T))

def dfdT(nu, T):
    """ Derivative of specific helmholtz free energy [kJ / kg K]
    w.r.t temperature at constant specific volume"""

    return -s(nu, T)
def dPdT(nu, T):
    """ Derivative of pressure [kJ / m^3 K]
    w.r.t temperature at constant specific volume"""

    return P(nu, T) * a(nu, T)
def dudT(nu, T):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant specific volume"""

    return cv(nu, T)
def dsdT(nu, T):
    """ Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant specific volume"""

    return cv(nu, T) / T
def dhdT(nu, T):
    """ Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant specific volume"""

    return cv(nu, T) + P(nu, T) * (10**6 / 1000) * nu * a(nu, T)

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

#### region 3 properties ####


#### region 3 property derivatives ####


###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

#### region 3 properties ####


#### region 3 property derivatives ####
