"""Provides basic equations for dimensionless specific Gibbs free energy, gamma(pi, theta), and derivatives in region 1."""

# type annotations
from __future__ import annotations

###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# Region 1, forwards equations for f(P, T)
_I = [  0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   2,   2,   2,
        2,   2,   3,   3,   3,   4,   4,   4,   5,   8,   8,  21,  23,  29,  30,  31,  32]
_J = [ -2,  -1,   0,   1,   2,   3,   4,   5,  -9,  -7,  -1,   0,   1,   3,  -3,   0,   1,
        3,  17,  -4,   0,   6,  -5,  -2,  10,  -8, -11,  -6, -29, -31, -38, -39, -40, -41]
_n = [ 0.14632971213167e+00, -0.84548187169114e+00, -0.37563603672040e+01,  0.33855169168385e+01,
      -0.95791963387872e+00,  0.15772038513228e+00, -0.16616417199501e-01,  0.81214629983568e-03,
       0.28319080123804e-03, -0.60706301565874e-03, -0.18990068218419e-01, -0.32529748770505e-01,
      -0.21841717175414e-01, -0.52838357969930e-04, -0.47184321073267e-03, -0.30001780793026e-03,
       0.47661393906987e-04, -0.44141845330846e-05, -0.72694996297594e-15, -0.31679644845054e-04,
      -0.28270797985312e-05, -0.85205128120103e-09, -0.22425281908000e-05, -0.65171222895601e-06,
      -0.14341729937924e-12, -0.40516996860117e-06, -0.12734301741641e-08, -0.17424871230634e-09,
      -0.68762131295531e-18,  0.14478307828521e-19,  0.26335781662795e-22, -0.11947622640071e-22,
       0.18228094581404e-23, -0.93537087292458e-25]
Ps = 16.53    # [Mpa      ]
Ts = 1386.0   # [K        ]
R  = 0.461526 # [kJ / kg K]

def gamma(pi: float, tau: float) -> float:
    """Dimensionless form for the specific Gibbs free energy.
    Reference: Equation (7) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += ni * (7.1 - pi)**Ii * (tau - 1.222)**Ji
    return sum

def gamma_pi(pi: float, tau: float) -> float:
    """Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi).
    Reference: Table (4) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (tau - 1.222)**Ji
    return sum

def gamma_pipi(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi).
    Reference: Table (4) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (tau - 1.222)**Ji
    return sum

def gamma_tau(pi: float, tau: float) -> float:
    """Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless temperature (tau).
    Reference: Table (4) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += ni * (7.1 - pi)**Ii * Ji * (tau - 1.222)**(Ji - 1)
    return sum

def gamma_tautau(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless temperature (tau).
    Reference: Table(4) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * (tau - 1.222)**(Ji - 2)
    return sum

def gamma_pitau(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi) and temperature (tau).
    Reference: Table(4) from R7-97(2012)"""
    sum = 0.0
    for Ii, Ji, ni in zip(_I, _J, _n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (tau - 1.222)**(Ji - 1)
    return sum
