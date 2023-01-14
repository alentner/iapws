# type annotations
from __future__ import annotations
from typing import cast, TYPE_CHECKING

# internal libraries
from .support import conversion

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

@conversion(_pressure)
def P(P: float) -> float:
    """Pressure [psi] -> [MPa]"""
    return P

def T(T: float, english: bool = True) -> float:
    """Temperature [F] -> [K]"""
    if english:
        return 5 * (T - 32) / 9 + 273.15
    else:
        return 9 * (T - 273.15) / 5 + 32

@conversion(_energy / (_mass * _temperature))
def g(g: float) -> float:
    """Specific gibbs free energy [Btu / lbm F] -> [kJ / kg K]"""
    return g

@conversion(_volume / _mass)
def v(v: float) -> float:
    """Specific volume [ft3 / lbm] -> [m^3 / kg]"""
    return v

@conversion(_energy / _mass)
def u(u: float) -> float:
    """Specific internal energy [Btu / lbm] -> [kJ / kg]"""
    return u

@conversion(_energy / (_mass * _temperature))
def s(s: float) -> float:
    """Specific entropy [Btu / lbm F] -> [kJ / kg K]"""
    return s

@conversion(_energy / _mass)
def h(h: float) -> float:
    """Specific enthalpy [Btu / lbm] -> [kJ / kg]"""
    return h

@conversion(_energy / (_mass * _temperature))
def cp(cp: float) -> float:
    """Specific isobaric heat capacity [Btu / lbm F] -> [kJ / kg K]"""
    return cp

@conversion(_energy / (_mass * _temperature))
def cv(cv: float) -> float:
    """Specific isochoric heat capacity [Btu / lbm F] -> [kJ / kg K]"""
    return cv

@conversion(_length)
def w(w: float) -> float:
    """Speed of sound [ft / s] -> [m / s]"""
    return w

@conversion( 1 / _temperature)
def a(a: float) -> float:
    """Isobaric cubic expansion coefficient [1 / F] -> [1 / K]"""
    return a

@conversion(1 / _pressure)
def k(k: float) -> float:
    """Isothermal compressibility [1 / psi] -> [1 / MPa]"""
    return k

