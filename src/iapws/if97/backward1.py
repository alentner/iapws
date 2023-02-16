"""
Implements backward equations for IF97 region 1.

For the calculation of properties as function of p,h or of p,s without any 
iteration, the two backward equations require extremely good numerical 
consistency with the basic equation. The backward equations for this region
are expressed in dimensionless form, theta = T / Ts, and read as follows, 
where pi = P / Ps, theta = T / Ts, eta = h / hs, and sigma = s / ss.

    T / Ts = theta(pi, eta  ) = SUM( ni pi^Ii (eta   + 1)^Ji )
    T / Ts = theta(pi, sigma) = SUM( ni pi^Ii (sigma + 2)^Ji )

This module is restricted to providing the backward temperature relation,
theta(pi, [eta, sigma]), and the backward pressure and temperarure relations,
pi(eta, sigma) and theta(eta, sigme), in region 1. 

The principle reference for this module is the IAPWS Industrial Formulation 1997:
    
    International Association for the Properties of Water and Steam, IAPWS R7-97(2012),
    Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic 
    Properties of Water and Steam (2012), available from: http://www.iapws.org.

The dimensional forward functions and their derivatives, g(P, T), are defined 
elsewhere in this package.
"""

# Type annotations
from __future__ import annotations

# Define public interface
__all__ = ["T_Ph", "T_Ps", "P_hs", "T_hs"]

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

# Region 1 dimensionless backward equations

def _theta_Ph(pi: float, eta: float) -> float:
    """Dimensionless temperature,
    w.r.t. dimensionless enthalpy (eta).
    Reference: Equation (11) from R7-97(2012)"""
    I = [  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  3,  3,  4,  5,  6]
    J = [  0,  1,  2,  6, 22, 32,  0,  1,  2,  3,  4, 10, 32, 10, 32, 10, 32, 32, 32, 32]
    n = [-0.23872489924521e+03,  0.40421188637945e+03,  0.11349746881718e+03, -0.58457616048039e+01,
         -0.15285482413140e-03, -0.10866707695377e-05, -0.13391744872602e+02,  0.43211039183559e+02,
         -0.54010067170506e+02,  0.30535892203916e+02, -0.65964749423638e+01,  0.93965400878363e-02,
          0.11573647505340e-06, -0.25858641282073e-04, -0.40644363084799e-08,  0.66456186191635e-07,
          0.80670734103027e-10, -0.93477771213947e-12,  0.58265442020601e-14, -0.15020185953503e-16]
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (eta + 1.0)**Ji
    return sum

# Region 1 backward properties

def T_Ph(P: float, h: float) -> float:
    """Temperature [K].
    Reference: Equation (11) from R7-97(2012)"""
    Ps = 1.0  # [Mpa    ]
    Ts = 1.0  # [K      ]
    hs = 2500 # [kJ / kg]
    pi = P / Ps
    eta = h / hs
    return _theta_Ph(pi, eta) * Ts

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################

# Region 1 dimensionless backward equations

def _theta_Ps(pi: float, sigma: float) -> float:
    """Dimensionless temperature,
    w.r.t. dimensionless entropy (sigma).
    Reference: Equation (13) from R7-97(2012)"""
    I = [  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  4]
    J = [  0,  1,  2,  3, 11, 31,  0,  1,  2,  3, 12, 31,  0,  1,  2,  9, 31, 10, 32, 32]
    n = [ 0.17478268058307e+03,  0.34806930892873e+02,  0.65292584978455e+01,  0.33039981775489e+00,
         -0.19281382923196e-06, -0.24909197244573e-22, -0.26107636489332e+00,  0.22592965981586e+00,
         -0.64256463395226e-01,  0.78876289270526e-02,  0.35672110607366e-09,  0.17332496994895e-23,
          0.56608900654837e-03, -0.32635483139717e-03,  0.44778286690632e-04, -0.51322156908507e-09,
         -0.42522657042207e-25,  0.26400441360689e-12,  0.78124600459723e-28, -0.30732199903668e-30]
    sum = 0.0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * pi**Ii * (sigma + 2.0)**Ji
    return sum

# Region 1 backward properties

def T_Ps(P: float, s: float) -> float:
    """Temperature [K].
    Reference: Equation (13) from R7-97(2012)"""
    Ps = 1.0 # [Mpa      ]
    Ts = 1.0 # [K        ]
    ss = 1.0 # [kJ / kg K]
    pi = P / Ps
    sigma = s / ss
    return _theta_Ps(pi, sigma) * Ts

###########################################################
#####           Enthalpy-Entropy Formulation          #####
###########################################################

# Region 1 backward properties

def P_hs(h: float, s: float) -> float:
    """Placeholder for future functionality."""
    assert False, "Backward Equation P(h, s) is not implemented!"
    return 0.0

def T_hs(h: float, s: float) -> float:
    """Placeholder for future functionality."""
    assert False, "Backward Equation T(h, s) is not implemented!"
    return 0.0
