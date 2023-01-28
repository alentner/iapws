"""Provides the Auxiliary equation definitions for IF97 water properties."""

# type annotations
from __future__ import annotations

###########################################################
#####       Constants and Auxiliary Functions         #####
###########################################################

# Auxiliary Equation for the boundary between region2 and 3
_n = [ 0.34805185628969e+03, -0.11671859879975e+01,  0.10192970039326e-02,  0.57254459862746e+03,  0.13918839778870e+02]
Ps = 1.0 # [MPa]
Ts = 1.0 # [K  ]

def bnd23P(T: float) -> float:
    """Auxiliary equation for region 2 and 3 boundary pressure [MPa].
    Reference: Equation (5) from R7-97(2012)"""
    theta = T / Ts
    return Ps * (_n[0] + _n[1] * theta + _n[2] * theta**2)

def bnd23T(P: float) -> float:
    """Auxiliary equation for region 2 and 3 boundary temperature [K].
    Reference: Equation (5) from R7-97(2012)"""
    pi = P / Ps
    return Ts * (_n[3] + ((pi - _n[4]) / _n[2])**0.5)
