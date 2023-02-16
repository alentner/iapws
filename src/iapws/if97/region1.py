"""
Implements water properties for IF97 region 1.

The basic equation for this region is a fundamental equation for the specific
Gibbs free energy g. This equation is expressed in dimensionless form,
gamma = g / RT, and reads as follows, where pi = P / Ps and tau = Ts / T.

    g / RT = gamma(pi, tau) = SUM( ni (7.1 - pi)^Ii (tau - 1.222)^Ji )

This module is restricted to providing the following water properties 
given a pressure and temperature defined thermodynamic state in region 1. 

The principle references for this module are:
    
    International Association for the Properties of Water and Steam, 
    IAPWS R7-97(2012), Revised Release on the IAPWS Industrial Formulation 1997
    for the Thermodynamic Properties of Water and Steam (2012),
    available from: http://www.iapws.org.

The dimensionless forward and backward basic equations and their derivatives,
f(pi, tau) and fx(pi, [eta, sigma]) respectively, are defined elsewhere.
"""

# Type annotations
from __future__ import annotations

# Internal libraries
from .basic1 import (Ps, Ts, R, gamma, gamma_pi, gamma_tau, 
                     gamma_pipi, gamma_tautau, gamma_pitau)
from .unit import _kilo_mega

# Define public interface
__all__ = ["Pbnd0", "Pbnd1", "Tbnd0", "Tbnd1",
           "g", "v", "u", "s", "h", "cp", "cv", "w", "av", "kT", "f"]

###########################################################
#####       Range of Validity (Boundary Constants)    #####
###########################################################
Pbnd0 = 611.213e-03 # [MPa ]
Pbnd1 = 100.0       # [MPa ]
Tbnd0 = 273.15      # [K   ]
Tbnd1 = 623.15      # [K   ]

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

# Region 1 forward properties

def g(P: float, T: float) -> float:
    """Specific gibbs free energy [kJ / kg].
    Reference: Equation (7) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return gamma(pi, tau) * R * T

def v(P: float, T: float) -> float:
    """Specific volume [m^3 / kg].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return pi * gamma_pi(pi, tau) * R * T / (P * _kilo_mega)

def u(P: float, T: float) -> float:
    """Specific internal energy [kJ / kg].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - pi * gamma_pi(pi, tau)) * R * T    

def s(P: float, T: float) -> float:
    """Specific entropy [kJ / kg K].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return (tau * gamma_tau(pi, tau) - gamma(pi, tau)) * R

def h(P: float, T: float) -> float:
    """Specific enthalpy [kJ / kg].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return tau * gamma_tau(pi, tau) * R * T

def cp(P: float, T: float) -> float:
    """Specific isobaric heat capacity [kJ / kg K].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return -tau**2 * gamma_tautau(pi, tau) * R

def cv(P: float, T: float) -> float:
    """Specific isochoric heat capacity [kJ / kg K].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return (-tau**2 * gamma_tautau(pi, tau) + (gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / gamma_pipi(pi, tau)) * R

def w(P: float, T: float) -> float:
    """Speed of sound [m / s].
    Reference: Table (3) from R7-97(2012)"""
    pi = P / Ps
    tau = Ts / T
    return (gamma_pi(pi, tau)**2 * R * T * _kilo_mega / ((gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau)) - gamma_pipi(pi, tau)))**0.5

def av(P: float, T: float) -> float:
    """Isobaric cubic expansion coefficient [1 / K].
    Reference: Equation (6) from AN3-07(2018)"""
    pi = P / Ps
    tau = Ts / T
    return (1 - tau * gamma_pitau(pi, tau) / gamma_pi(pi, tau)) / T

def kT(P: float, T: float) -> float:
    """Isothermal compressibility [m^3 / kJ].
    Reference: Equation (6) from AN3-07(2018)"""
    pi = P / Ps
    tau = Ts / T
    return -(pi * gamma_pipi(pi, tau) / gamma_pi(pi, tau)) / (P * _kilo_mega)

def f(P: float, T: float) -> float:
    """Specific helmholtz free energy [kJ / kg].
    Reference: Basic relation -> f = u - T * s"""
    return u(P, T) - T * s(P, T)
