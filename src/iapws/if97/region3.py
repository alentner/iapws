"""Provides the definitions for thermodymic properties and auxillaries for region 3"""

# type annotations
from __future__ import annotations
from typing import cast, TYPE_CHECKING

# system libraries
from math import log

# internal libraries
from .

# static analysis
if TYPE_CHECKING:
    from typing import Any, Callable, TypeVar
    F = TypeVar('F', bound = Callable[..., Any])
    D = Callable[[F], F]

# deal w/ runtime cast
else:
    F = None

###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# constants and non-dimenionalization;
# Region 3, forwards equations for (v, T)
I   = [ 0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3, 
        3,  3,  4,  4,  4,  4,  5,  5,  5,  6,  6,  6,  7,  8,  9,  9, 10, 10, 11]
J   = [ 0,  1,  2,  7, 10, 12, 23,  2,  6, 15, 17,  0,  2,  6,  7, 22, 26,  0,  2,  4,
       16, 26,  0,  2,  4, 26,  1,  3, 26,  0,  2, 26,  2, 26,  2, 26,  0,  1, 26]
n   = [-0.15732845290239e+02,  0.20944396974307e+02, -0.76867707878716e+01,  0.26185947787954e+01, -0.28080781148620e+01,
        0.12053369696517e+01, -0.84566812812502e-02, -0.12654315477714e+01, -0.11524407806681e+01,  0.88521043984318e+00,
       -0.64207765181607e+00,  0.38493460186671e+00, -0.85214708824206e+00,  0.48972281541877e+01, -0.30502617256965e+01,
        0.39420536879154e-01,  0.12558408424308e+00, -0.27999329698710e+00,  0.13899799569460e+01, -0.20189915023570e+01,
       -0.82147637173963e-02, -0.47596035734923e+00,  0.43984074473500e-01, -0.44476435428739e+00,  0.90572070719733e+00,
        0.70522450087967e+00,  0.10770512626332e+00, -0.32913623258954e+00, -0.50871062041158e+00, -0.22175400873096e-01,
        0.94260751665092e-01,  0.16436278447961e+00, -0.13503372241348e-01, -0.14834345352472e-01,  0.57922953628084e-03,
        0.32308904703711e-02,  0.80964802996215e-04, -0.16557679795037e-03, -0.44923899061815e-04]
n0  =   0.10658070028513e+01
Ps  = 22.064   # [Mpa      ]
Ts  = 647.096  # [K        ]
R   = 0.461526 # [kJ / kg K]
rhos = 322     # [kg / m^3 ]

# constants for region boundary;
# Regions 2 and 3, forwards equations for P or T
n23 = [ 0.34805185628969e+03, -0.11671859879975e+01,  0.10192970039326e-02,  0.57254459862746e+03,  0.13918839778870e+02]
Pb  = 1.00    # [Mpa      ]
Tb  = 1.00    # [K        ]

# constants for subregion boundaries;
# Subregions 3ab-3rx, backwards equations for v(P, T)
Iab = [ 0,  1,  2, -1, -2]
Icd = [ 0,  1,  2,  3]
Igh = [ 0,  1,  2,  3,  4]
Iij = [ 0,  1,  2,  3,  4]
Ijk = [ 0,  1,  2,  3,  4]
Imn = [ 0,  1,  2,  3]
Iop = [ 0,  1,  2, -1, -2]
Iqu = [ 0,  1,  2,  3]
Irx = [ 0,  1,  2,  3]

nab = [ 0.154793642129415e+00, -0.187661219490113e+03,  0.213144632222113e+02, -0.191887498864292e+04,  0.918419702359447e+03]
ncd = [ 0.585276966696349e+03,  0.278233532206915e+01, -0.127283549295878e-01,  0.159090746562729e-03] 
ngh = [-0.249284240900418e+05,  0.428143584791546e+04, -0.269029173140130e+03,  0.751608051114157e+01, -0.787105249910383e-01]
nij = [ 0.584814781649163e+03, -0.616179320924617e+00,  0.260763050899562e+00, -0.587071076864459e-02,  0.515308185433082e-04]
njk = [ 0.617229772068439e+03, -0.770600270141675e+00,  0.697072596851896e+00, -0.157391839848015e+01,  0.137897492684194e-03]
nmn = [ 0.535339483742384e+03,  0.761978122720128e+01, -0.158365725441648e+00,  0.137897492648194e-03]
nop = [ 0.969461372400213e+03, -0.332500170441278e+03,  0.642859598466067e+02,  0.773845935768222e+03, -0.152313732937084e+04]
nqu = [ 0.565603648239126e+03,  0.529062258221222e+01, -0.102020639611016e+00,  0.122240301070145e-02]
nrx = [ 0.584561202520006e+03, -0.102961025163669e+01,  0.243293362700452e+00, -0.294905044740799e-02]
nuv = [ 0.528199646263062e+03,  0.890579602135307e+01, -0.222814134903755e+00,  0.286791682263697e-02]

# Boundaries defining region 3, and subregions 3a-3t
def bnd23P(T: float) -> float:
    """Auxiliary Equation for the Boundary between Regions 2 and 3;
    determine the boundary pressure line using temperature.
    Reference: Equation (5) from R7-97(2012)"""
    theta = T / Tb
    return Pb * (n23[0] + n23[1] * theta + n23[2] * theta**2)

def bnd23T(P: float) -> float:
    """Auxiliary Equation for the Boundary between Regions 2 and 3;
    determine the boundary temperature line using pressure.
    Reference: Equation (6) from R7-97(2012)"""
    pi = P / Pb
    return Tb * (n23[3] + ((pi - n23[4]) / n23[2])**0.5)

def _bnd3xy(Ixy: list[int], nxy: list[float]) -> D:
    """Usefull decorator factory to implement subregion x and y boundaries;
    determine the boundary temperature line using pressure.
    Reference: Equation (1) from SR5-05(2016)"""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(P: float) -> float:
            pi = P / Pb
            sum = 0.0
            for Ii, ni in zip(Ixy, nxy):
                sum += ni * pi**Ii
            return Tb * sum
        return cast(F, wrapper)
    return decorator

def bnd3abT(P: float) -> float:
    """Subregion boundary equation between subregions a and b;
    determine the boundary temperature line using pressure.
    Reference: Equation (2) from SR5-05(2016)"""
    pi = P / Pb
    sum = 0.0
    for Ii, ni in zip(Iab, nab):
        sum += ni * log(pi)**Ii
    return Tb * sum

@_bnd3xy(Icd, ncd)
def bnd3cdT(P: float) -> float:
    """Subregion boundary equation between subregions c and d"""
    pass

def bnd3efT(P: float) -> float:
    """Subregion boundary equation between subregions e and f;
    determine the boundary temperature line using pressure.
    Reference: Equation (3) from SR5-05(2016)"""
    pi = P / Pb
    return Tb * (3.727888004 * (pi - 22.064) + 647.096)

@_bnd3xy(Igh, ngh)
def bnd3ghT(P: float) -> float:
    """Subregion boundary equation between subregions g and h"""
    pass

@_bnd3xy(Iij, nij)
def bnd3ijT(P: float) -> float:
    """Subregion boundary equation between subregions i and j"""
    pass

@_bnd3xy(Ijk, njk)
def bnd3jkT(P: float) -> float:
    """Subregion boundary equation between subregions j and k"""
    pass

@_bnd3xy(Imn, nmn)
def bnd3mnT(P: float) -> float:
    """Subregion boundary equation between subregions m and n"""
    pass

def bnd3opT(P: float) -> float:
    """Subregion boundary equation between subregions o and p;
    determine the boundary temperature line using pressure.
    Reference: Equation (2) from SR5-05(2016)"""
    pi = P / Pb
    sum = 0.0
    for Ii, ni in zip(Iop, nop):
        sum += ni * log(pi)**Ii
    return Tb * sum

@_bnd3xy(Iqu, nqu)
def bnd3quT(P: float) -> float:
    """Subregion boundary equation between subregions q and u"""
    pass

@_bnd3xy(Irx, nrx)
def bnd3rxT(P: float) -> float:
    """Subregion boundary equation between subregions r and x"""
    pass

def idRegion_T(P: float, T: float) -> float:
    """Identification of subregion from IF97 specification
    using pressure and temperature as primary varibles.
    Reference: Table (2) from SR5-05(2016)"""

    # Constant boundaries
    Pbnd3cd = 1.900881189173929e+01

    region = 0

    if 40.0 < P <= 100.0:
        region = 1 if T <= bnd3abT(P) else 2

    elif 25.0 < P <= 40.0:
        Tbnd3ab = bnd3abT(P)
        Tbnd3cd = bnd3cdT(P)
        Tbnd3ef = bnd3efT(P)
        if T <= Tbnd3cd: region = 3
        elif Tbnd3cd < T <= Tbnd3ab: region = 4
        elif Tbnd3ab < T <= Tbnd3ef: region = 5
        elif T > Tbnd3ef: region = 6
        else: pass

    elif 23.5 < P <= 25.0


    elif 23.0 < P <= 23.5


    elif 22.5 < P <= 23.0


    elif 

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
def f(v, T):
    """Specific Helmholtz free energy [kJ / kg]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return phi(delta, tau) * R * T

def P(v, T):
    """Pressure [Mpa]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return delta * phi_delta(delta, tau) * R * T * rho / 1e+03

def u(v, T):
    """Specific internal energy [kJ / kg]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return tau * phi_tau(delta, tau) * R * T    

def s(v, T):
    """Specific entropy [kJ / kg K]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return (tau * phi_tau(delta, tau) - phi(delta, tau)) * R

def h(v, T):
    """Specific enthalpy [kJ / kg]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return (tau * phi_tau(delta, tau) + delta * phi_delta(delta, tau)) * R * T

def cp(v, T):
    """Specific isobaric heat capacity [kJ / kg K]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return (-tau**2 * phi_tautau(delta, tau) + (delta * phi_delta(delta, tau) - delta * tau * phi_deltatau(delta, tau))**2 / (2 * delta * phi_delta(delta, tau) + delta**2 * phi_deltadelta(delta, tau))) * R

def cv(v, T):
    """Specific isochoric heat capacity [kJ / kg K]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return (-tau**2 * phi_tautau(delta, tau)) * R

def w(v, T):
    """Speed of sound [m / s]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return ((2 * delta * phi_delta(delta, tau) + delta**2 * phi_deltadelta(delta, tau) - (delta * phi_delta(delta, tau) - delta * tau * phi_deltatau(delta, tau))**2 / (tau**2 * phi_tautau(delta, tau))) * R * T * 1000)**0.5

def av(v, T):
    """Isobaric cubic expansion coefficient [1 / K]"""
    return ap(v, T) / (v * bp(v, T))

def kT(v, T):
    """Isothermal compressibility [m^3 / kJ]"""
    return 1 / (v * P(v * T) * 1e+03 * bp(v, T))

def ap(v, T):
    """Relative pressure coefficient [1 / K]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T
    return (1 - tau * phi_deltatau(delta, tau) / phi_delta(delta, tau)) / T

def bp(v, T):
    """Isothermal stress coefficient [kg / m^3]"""
    rho = 1 / v
    delta = rho / rhos
    tau = Ts / T

    return (2 + delta * phi_deltadelta(delta, tau) / phi_delta(delta, tau)) * rho

#### region 3 property derivatives ####
def dfdv(v, T):
    """ Derivative of specific helmholtz free energy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""
    return -P(v, T) * 1e+03

def dPdv(v, T):
    """ Derivative of pressure [kJ kg / m^3 m^3]
    w.r.t specific volume at constant temperature"""
    return -P(v, T) * bp(v, T)

def dudv(v, T):
    """ Derivative of specific internal energy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""
    return P(v, T) * 1e+03 * (T * ap(v, T) - 1)

def dsdv(v, T):
    """ Derivative of specific entropy [kJ kg / kg K m^3]
    w.r.t specific volume at constant temperature"""
    return P(v, T) * 1e+03 * ap(v, T)

def dhdv(v, T):
    """ Derivative of specific enthalpy [kJ kg / kg m^3]
    w.r.t specific volume at constant temperature"""
    return P(v, T) * 1e+03 * (T * ap(v, T) - v * bp(v, T))

def dfdT(v, T):
    """ Derivative of specific helmholtz free energy [kJ / kg K]
    w.r.t temperature at constant specific volume"""
    return -s(v, T)

def dPdT(v, T):
    """ Derivative of pressure [kJ / m^3 K]
    w.r.t temperature at constant specific volume"""
    return P(v, T) * ap(v, T)

def dudT(v, T):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant specific volume"""
    return cv(v, T)

def dsdT(v, T):
    """ Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant specific volume"""
    return cv(v, T) / T

def dhdT(v, T):
    """ Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant specific volume"""
    return cv(v, T) + P(v, T) * 1e+03 * v * ap(v, T)

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
