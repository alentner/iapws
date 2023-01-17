""" Provides unit definitions for physical quantities and conversions """

# type annotations
from __future__ import annotations

# internal libraries
from .support import _conversion

# basic conversion constants
_ftlbm_btu = 7.37562E+2 # [ft lbf] -> [Btu]  
_sqft_sqin = 6.94444E-3 # [ft2   ] -> [in2]
_kilo_mega = 1.00000E+3 # [k_    ] -> [M_ ] 

# conversion of basic physical quantities
_mass   = 2.20462E+0 # [lbm] -> [kg ]
_length = 3.28084E+0 # [ft ] -> [m  ]
_volume = 3.53147E+1 # [ft3] -> [m3 ]
_energy = 9.47817E-1 # [Btu] -> [kJ ]
_temperature = 9 / 5 # [F  ] -> [K  ]

# conversion of compound physical quantities
_pressure = _ftlbm_btu / _volume * _sqft_sqin * _kilo_mega # [psi] -> [Mpa]

#### unit conversion functions ####
def T(T: float, /, *, english: bool = True) -> float:
    """Temperature (quantity) [F] -> [K]"""
    if english:
        return 5 * (T - 32) / 9 + 273.15
    else:
        return 9 * (T - 273.15) / 5 + 32

@_conversion(_pressure)
def P(P: float, /, *, english: bool = True) -> float:
    """Pressure [psi] -> [MPa]"""
    return P

@_conversion(_energy / (_mass * _temperature))
def g(g: float, /, *, english: bool = True) -> float:
    """Specific gibbs free energy [Btu / lbm F] -> [kJ / kg K]"""
    return g

@_conversion(_volume / _mass)
def v(v: float, /, *, english: bool = True) -> float:
    """Specific volume [ft3 / lbm] -> [m^3 / kg]"""
    return v

@_conversion(_energy / _mass)
def u(u: float, /, *, english: bool = True) -> float:
    """Specific internal energy [Btu / lbm] -> [kJ / kg]"""
    return u

@_conversion(_energy / (_mass * _temperature))
def s(s: float, /, *, english: bool = True) -> float:
    """Specific entropy [Btu / lbm F] -> [kJ / kg K]"""
    return s

@_conversion(_energy / _mass)
def h(h: float, /, *, english: bool = True) -> float:
    """Specific enthalpy [Btu / lbm] -> [kJ / kg]"""
    return h

@_conversion(_energy / (_mass * _temperature))
def cp(cp: float, /, *, english: bool = True) -> float:
    """Specific isobaric heat capacity [Btu / lbm F] -> [kJ / kg K]"""
    return cp

@_conversion(_energy / (_mass * _temperature))
def cv(cv: float, /, *, english: bool = True) -> float:
    """Specific isochoric heat capacity [Btu / lbm F] -> [kJ / kg K]"""
    return cv

@_conversion(_length)
def w(w: float, /, *, english: bool = True) -> float:
    """Speed of sound [ft / s] -> [m / s]"""
    return w

@_conversion( 1 / _temperature)
def av(av: float, /, *, english: bool = True) -> float:
    """Isobaric cubic expansion coefficient [1 / F] -> [1 / K]"""
    return av

@_conversion(_kilo_mega / _pressure)
def kT(kT: float, /, *, english: bool = True) -> float:
    """Isothermal compressibility [1 / psi] -> [m^3 / kJ]"""
    return kT

#### unit conversion (derivatives) functions ####
@_conversion(_temperature)
def dT(T: float, /, *, english: bool = True) -> float:
    """Temperature (scale) [F] -> [K]"""
    return T

@_conversion(_pressure)
def dP(P: float, /, *, english: bool = True) -> float:
    """Pressure (scale) [psi] -> [MPa]"""
    return P

@_conversion(_energy / (_mass * _temperature))
def dg(g: float, /, *, english: bool = True) -> float:
    """Specific gibbs free energy (scale) [Btu / lbm F] -> [kJ / kg K]"""
    return g

@_conversion(_volume / _mass)
def dv(v: float, /, *, english: bool = True) -> float:
    """Specific volume (scale) [ft3 / lbm] -> [m^3 / kg]"""
    return v

@_conversion(_energy / _mass)
def du(u: float, /, *, english: bool = True) -> float:
    """Specific internal energy (scale) [Btu / lbm] -> [kJ / kg]"""
    return u

@_conversion(_energy / (_mass * _temperature))
def ds(s: float, /, *, english: bool = True) -> float:
    """Specific entropy (scale) [Btu / lbm F] -> [kJ / kg K]"""
    return s

@_conversion(_energy / _mass)
def dh(h: float, /, *, english: bool = True) -> float:
    """Specific enthalpy (scale) [Btu / lbm] -> [kJ / kg]"""
    return h

#### unit conversion (reciprical derivatives) functions ####
@_conversion(1 / _temperature)
def d_dT(T: float, /, *, english: bool = True) -> float:
    """Temperature (inverse) [1 / F] -> [1 / K]"""
    return 1 / T

@_conversion(1 / _pressure)
def d_dP(P: float, /, *, english: bool = True) -> float:
    """Pressure (inverse) [1 / psi] -> [1 / MPa]"""
    return 1 / P

@_conversion(_mass / _energy)
def d_dh(h: float, /, *, english: bool = True) -> float:
    """Specific enthalpy (inverse) [lbm / Btu] -> [kg / kJ]"""
    return 1 / h

@_conversion(_mass * _temperature / _energy)
def d_ds(s: float, /, *, english: bool = True) -> float:
    """Specific entropy (inverse) [lbm F / Btu] -> [kg K / kJ]"""
    return 1 / s
