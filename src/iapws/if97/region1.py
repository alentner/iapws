"""Provides the thermodynamic properties and derivatives for region 1"""

# type annotations
from __future__ import annotations

# internal libraries
from .basic1 import gamma, gamma_pi, gamma_pipi, gamma_tau, gamma_tautau, gamma_pitau
from .basic1 import Ps, Ts, R
from .backward1 import theta_h, theta_s
from .backward1 import Ps_h, Ts_h, hs_h, Ps_s, Ts_s, ss_s
from .unit import _kilo_mega

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

#### region 1 properties ####
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

#### region 1 property derivatives ####
def dPdP(P: float, T: float) -> float:
    """Derivative of pressure [kJ m^3 / m^3 kJ],
    w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return 1.0

def dTdP(P: float, T: float) -> float:
    """Derivative of temperature [K m^3 / kJ],
    w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return 0.0

def dvdP(P: float, T: float) -> float:
    """Derivative of specific volume [m^3 m^3 / kg kJ],
    w.r.t pressure at constant temperature
    Reference: Table (2) from AN3-07(2018)"""
    return -v(P, T) * kT(P, T)

def dudP(P: float, T: float) -> float:
    """Derivative of specific internal energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return v(P, T) * (P * _kilo_mega * kT(P, T) - T * av(P, T))

def dhdP(P: float, T: float) -> float:
    """Derivative of specific enthalpy [kJ m^3 / kg kJ],
    w.r.t pressure at constant temperature
    Reference: Table (2) from AN3-07(2018)"""
    return v(P, T) * (1 - T * av(P, T))

def dsdP(P: float, T: float) -> float:
    """Derivative of specific entropy [kJ m^3 / kg K kJ],
    w.r.t pressure at constant temperature
    Reference: Table (2) from AN3-07(2018)"""
    return -v(P, T) * av(P, T)

def dgdP(P: float, T: float) -> float:
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return v(P, T)

def dfdP(P: float, T: float) -> float:
    """Derivative of specific helmholtz free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return (P * _kilo_mega) * v(P, T) * kT(P, T)

def dPdT(P: float, T: float) -> float:
    """Derivative of pressure [kJ / m^3 K],
    w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return 0.0

def dTdT(P: float, T: float) -> float:
    """Derivative of temperature [K / K],
    w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return 1.0

def dvdT(P: float, T: float) -> float:
    """Derivative of specific volume [m^3 / kg K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return v(P, T) * av(P, T)

def dudT(P: float, T: float) -> float:
    """Derivative of specific internal energy [kJ / kg K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return cp(P, T) - (P * _kilo_mega) * v(P, T) * av(P, T)

def dhdT(P: float, T: float) -> float:
    """Derivative of specific enthalpy [kJ / kg K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return cp(P, T)

def dsdT(P: float, T: float) -> float:
    """Derivative of specific entropy [kJ / kg K K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return cp(P, T) / T

def dgdT(P: float, T: float) -> float:
    """Derivative of specific gibbs free energy [kJ / kg K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return -s(P, T)

def dfdT(P: float, T: float) -> float:
    """Derivative of specific helmholtz free energy [kJ / kg K],
    w.r.t temperature at constant pressure
    Reference: Table (2) from AN3-07(2018)"""
    return -(P * _kilo_mega) * v(P, T) * av(P, T) - s(P, T)

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

#### region 1 properties ####
def g_h(P: float, h: float) -> float:
    """Specific gibbs free energy [kJ / kg]."""
    return g(P, T_h(P, h))

def v_h(P: float, h: float) -> float:
    """Specific volume [m^3 / kg]."""
    return v(P, T_h(P, h))

def u_h(P: float, h: float) -> float:
    """Specific internal energy [kJ / kg]."""
    return u(P, T_h(P, h))

def s_h(P: float, h: float) -> float:
    """Specific entropy [kJ / kg K]."""
    return s(P, T_h(P, h))

def T_h(P: float, h: float) -> float:
    """Temperature [K].
    Reference: Equation (11) from R7-97(2012)"""
    pi = P / Ps_h
    eta = h / hs_h
    return theta_h(pi, eta) * Ts_h

def cp_h(P: float, h: float) -> float:
    """Specific isobaric heat capacity [kJ / kg K]."""
    return cp(P, T_h(P, h))

def cv_h(P: float, h: float) -> float:
    """Specific isochoric heat capacity [kJ / kg K]."""
    return cv(P, T_h(P, h))

def w_h(P: float, h: float) -> float:
    """Speed of sound [m / s]."""
    return w(P, T_h(P, h))

def av_h(P: float, h: float) -> float:
    """Isobaric cubic expansion coefficient [1 / K]."""
    return av(P, T_h(P, h))

def kT_h(P: float, h: float) -> float:
    """Isothermal compressibility [m^3 / kJ]."""
    return kT(P, T_h(P, h))

def f_h(P: float, h: float) -> float:
    """Specific Helmholtz free energy [kJ / kg]."""
    return f(P, T_h(P, h))

#### region 1 property derivatives ####
def dPdP_h(P: float, h: float) -> float:
    """Derivative of pressure [kJ m^3 / m^3 kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    return 1.0

def dTdP_h(P: float, h: float) -> float:
    """Derivative of temperature [K m^3 / kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return -dhdP(P, T) / dhdT(P, T)

def dvdP_h(P: float, h: float) -> float:
    """Derivative of specific volume [m^3 m^3 / kg kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dvdP(P, T) - dvdT(P, T) * dhdP(P, T) / dhdT(P, T)

def dudP_h(P: float, h: float) -> float:
    """Derivative of specific internal energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dudP(P, T) - dudT(P, T) * dhdP(P, T) / dhdT(P, T)

def dhdP_h(P: float, h: float) -> float:
    """Derivative of specific enthalpy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    return 0.0

def dsdP_h(P: float, h: float) -> float:
    """Derivative of specific entropy [kJ m^3 / kg K kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dsdP(P, T) - dsdT(P, T) * dhdP(P, T) / dhdT(P, T)

def dgdP_h(P: float, h: float) -> float:
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dgdP(P, T) - dgdT(P, T) * dhdP(P, T) / dhdT(P, T)

def dfdP_h(P: float, h: float) -> float:
    """Derivative of specific helmholtz free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific enthalpy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dfdP(P, T) - dfdT(P, T) * dhdP(P, T) / dhdT(P, T)

def dPdh_h(P: float, h: float) -> float:
    """Derivative of pressure [kJ kg /m^3 kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    return 0.0

def dTdh_h(P: float, h: float) -> float:
    """Derivative of temperature [K kg / kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return 1 / dhdT(P, T)

def dvdh_h(P: float, h: float) -> float:
    """Derivative of specific volume [m^3 kg / kg kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dvdT(P, T) / dhdT(P, T)

def dudh_h(P: float, h: float) -> float:
    """Derivative of specific internal energy [kJ kg / kg kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dudT(P, T) / dhdT(P, T)

def dhdh_h(P: float, h: float) -> float:
    """Derivative of specific enthalpy [kJ kg / kg kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    return 1.0

def dsdh_h(P: float, h: float) -> float:
    """Derivative of specific entropy [kJ kg / kg K kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dsdT(P, T) / dhdT(P, T)

def dgdh_h(P: float, h: float) -> float:
    """Derivative of specific gibbs free energy [kJ kg / kg kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dgdT(P, T) / dhdT(P, T)

def dfdh_h(P: float, h: float) -> float:
    """Derivative of specific helmholtz free energy [kJ kg / kg kJ],
    w.r.t specific enthalpy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_h(P, h)
    return dfdT(P, T) / dhdT(P, T)

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################

#### region 1 properties ####
def g_s(P, s):
    (P: float, s: float) -> float:
    """Specific gibbs free energy [kJ / kg]."""
    return g(P, T_s(P, s))

def v_s(P: float, s: float) -> float:
    """Specific volume [m^3 / kg]."""
    return v(P, T_s(P, s))

def u_s(P: float, s: float) -> float:
    """Specific internal energy [kJ / kg]."""
    return u(P, T_s(P, s))   

def T_s(P: float, s: float) -> float:
    """Temperature [K].
    Reference: Equation (13) from R7-97(2012)"""
    pi = P / Ps_s
    sigma = s / ss_s
    return theta_s(pi, sigma) * Ts_bs

def h_s(P: float, s: float) -> float:
    """Specific enthalpy [kJ / kg]."""
    return h(P, T_s(P, s))

def cp_s(P: float, s: float) -> float:
    """Specific isobaric heat capacity [kJ / kg K]."""
    return cp(P, T_s(P, s))

def cv_s(P: float, s: float) -> float:
    """Specific isochoric heat capacity [kJ / kg K]."""
    return cv(P, T_s(P, s))

def w_s(P: float, s: float) -> float:
    """Speed of sound [m / s]."""
    return w(P, T_s(P, s))

def av_s(P: float, s: float) -> float:
    """Isobaric cubic expansion coefficient [1 / K]."""
    return av(P, T_s(P, s))

def kT_s(P: float, s: float) -> float:
    """Isothermal compressibility [m^3 / kJ]."""
    return kT(P, T_s(P, s))

def f_s(P: float, s: float) -> float:
    """Specific Helmholtz free energy [kJ / kg]."""
    return f(P, T_s(P, s))

#### region 1 property derivatives ####
def dPdP_s(P: float, s: float) -> float:
    """Derivative of pressure [kJ m^3 /m^3 kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    return 1.0

def dTdP_s(P: float, s: float) -> float:
    """Derivative of temperature [K m^3 / kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return -dsdP(P, T) / dsdT(P, T)

def dvdP_s(P: float, s: float) -> float:
    """Derivative of specific volume [m^3 m^3 / kg kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dvdP(P, T) - dvdT(P, T) * dsdP(P, T) / dsdT(P, T)

def dudP_s(P: float, s: float) -> float:
    """Derivative of specific internal energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dudP(P, T) - dudT(P, T) * dsdP(P, T) / dsdT(P, T)

def dhdP_s(P: float, s: float) -> float:
    """Derivative of specific enthalpy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dhdP(P, T) - dhdT(P, T) * dsdP(P, T) / dsdT(P, T)

def dsdP_s(P: float, s: float) -> float:
    """Derivative of specific entropy [kJ m^3 / kg K kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    return 0.0

def dgdP_s(P: float, s: float) -> float:
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dgdP(P, T) - dgdT(P, T) * dsdP(P, T) / dsdT(P, T)

def dfdP_s(P: float, s: float) -> float:
    """Derivative of specific helmholtz free energy [kJ m^3 / kg kJ],
    w.r.t pressure at constant specific entropy.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dfdP(P, T) - dfdT(P, T) * dsdP(P, T) / dsdT(P, T)

def dPds_s(P: float, s: float) -> float:
    """Derivative of pressure [kJ kg K/m^3 kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    return 0.0

def dTds_s(P: float, s: float) -> float:
    """Derivative of temperature [K kg K/ kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return 1 / dsdT(P, T)

def dvds_s(P: float, s: float) -> float:
    """Derivative of specific volume [m^3 kg K / kg kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dvdT(P, T) / dsdT(P, T)

def duds_s(P: float, s: float) -> float:
    """Derivative of specific internal energy [kJ kg K / kg kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dudT(P, T) / dsdT(P, T)

def dhds_s(P: float, s: float) -> float:
    """Derivative of specific enthalpy [kJ kg K / kg kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dhdT(P, T) / dsdT(P, T)

def dsds_s(P: float, s: float) -> float:
    """Derivative of specific entropy [kJ kg K / kg K kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    return 1.0

def dgds_s(P: float, s: float) -> float:
    """Derivative of specific gibbs free energy [kJ kg K / kg kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dgdT(P, T) / dsdT(P, T)

def dfds_s(P: float, s: float) -> float:
    """Derivative of specific helmholtz free energy [kJ kg K / kg kJ],
    w.r.t specific entropy at constant pressure.
    Reference: Equation (5) from AN3-07(2018)"""
    T = T_s(P, s)
    return dfdT(P, T) / dsdT(P, T)
