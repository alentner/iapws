"""
Implements partial derivatives for IF97 gibbs free energy regions.

The basic equation for these regions is a fundamental equation for the specific 
Gibbs free energy g. This partial derivatives are expressed in general form, 
where z is the property, x is the partial derivative, and y is held constant.

                      (dz/dP)T * (dy/dT)P - (dz/dT)P * (dy/dP)T
    (dz/dx)y (P, T) = -----------------------------------------
                      (dx/dP)T * (dy/dT)P - (dx/dT)P * (dy/dP)T

This module is restricted to providing the partial derivatives for the water
properties provided elsewhere based on the Gibbs free energy, regions 1, 2, and 5. 

The principle references for this module are:
    
    International Association for the Properties of Water and Steam, IAPWS R7-97(2012),
    Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic 
    Properties of Water and Steam (2012), available from: http://www.iapws.org.

    International Association for the Properties of Water and Steam,
    AN3-07(2018), Revised Advisory Note No. 3: Thermodynamic Derivatives from
    IAPWS Formulations (2018), available from: http://www.iapws.org.

The dimensional forward and backward functions and their derivatives, 
g(P, T) and [T, P]([P], [h, s]) respectively, are defined elsewhere in this package.
"""
# type annotations
from __future__ import annotations

# system libraries
import sys

# static analysis
from typing import TypeVar
from types import ModuleType
M = TypeVar('M', bound=ModuleType)

# module access and module level @property(s)
this = sys.modules[__name__]

def dzdx_y(P: float, T: float, *, equation: M, z: str, x: str, y: str) -> float:
    """General Derivative Procedure for Gibbs Energy Equations."""
    dxdT_P = getattr(this, f'd{x}dT_P')(P, T, equation=equation)
    dxdP_T = getattr(this, f'd{x}dP_T')(P, T, equation=equation)
    dydT_P = getattr(this, f'd{y}dT_P')(P, T, equation=equation)
    dydP_T = getattr(this, f'd{y}dP_T')(P, T, equation=equation)
    dzdT_P = getattr(this, f'd{z}dT_P')(P, T, equation=equation)
    dzdP_T = getattr(this, f'd{z}dP_T')(P, T, equation=equation)
    return (dzdP_T * dydT_P - dzdT_P * dydP_T) / (dxdP_T * dydT_P - dxdT_P * dydP_T)

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

def dPdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of pressure [kJ / m^3 K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return 0.0

def dTdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of temperature [K / K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return 1.0

def dvdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific volume [m^3 / kg K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.v(P, T) * equation.av(P, T)

def dudT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific internal energy [kJ / kg K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.cp(P, T) - (P * _kilo_mega) * v(P, T) * av(P, T)

def dhdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific enthalpy [kJ / kg K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.cp(P, T)

def dsdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific entropy [kJ / kg K K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.cp(P, T) / T

def dgdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific gibbs free energy [kJ / kg K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return -equation.s(P, T)

def dfdT_P(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific helmholtz free energy [kJ / kg K], w.r.t temperature at constant pressure.
    Reference: Table (2) from AN3-07(2018)"""
    return -(P * _kilo_mega) * equation.v(P, T) * equation.av(P, T) - equation.s(P, T)

def dPdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of pressure [kJ m^3 / m^3 kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return 1.0

def dTdP_T(equation: float, T: float, *, equation: M) -> float:
    """Derivative of temperature [K m^3 / kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return 0.0

def dvdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific volume [m^3 m^3 / kg kJ], w.r.t pressure at constant temperature
    Reference: Table (2) from AN3-07(2018)"""
    return -equation.v(P, T) * equation.kT(P, T)

def dudP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific internal energy [kJ m^3 / kg kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.v(P, T) * (P * _kilo_mega * equation.kT(P, T) - T * equation.av(P, T))

def dhdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific enthalpy [kJ m^3 / kg kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.v(P, T) * (1 - T * equation.av(P, T))

def dsdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific entropy [kJ m^3 / kg K kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return -equation.v(P, T) * equation.av(P, T)

def dgdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return equation.v(P, T)

def dfdP_T(P: float, T: float, *, equation: M) -> float:
    """Derivative of specific helmholtz free energy [kJ m^3 / kg kJ], w.r.t pressure at constant temperature.
    Reference: Table (2) from AN3-07(2018)"""
    return (P * _kilo_mega) * equation.v(P, T) * equation.kT(P, T)
