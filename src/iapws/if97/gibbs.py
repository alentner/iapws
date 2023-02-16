
# type annotations
from __future__ import annotations

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

def dPdP(region: M, P: float, T: float) -> float:
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
