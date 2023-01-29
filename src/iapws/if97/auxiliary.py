"""
Implements auxiliary basic equations for IF97.


The boundary between regions 2 and 3 is defined by the following simple
quadratic pressure-temperature relation, the B23-equation and its equivalent.
    
    pi = n1 + n2 theta + n3 theta^2
    theta = n4 + ( (pi - n5) / n3 )^0.5


This module is restricted to providing the dimensional auxiliary equation
for the boundary between regions 2 and 3. The principle reference for this module
is the IAPWS Industrial Formulation 1997:
    
    International Association for the Properties of Water and Steam, IAPWS R7-97(2012),
    Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic 
    Properties of Water and Steam (2012), available from: http://www.iapws.org.
"""

# type annotations
from __future__ import annotations

###########################################################
#####       Range of Validity (Boundary Constants)    #####
###########################################################
Pbnd0 = 16.5292 # [MPa ]
Pbnd1 = 100.0   # [MPa ]
Tbnd0 = 623.15  # [K   ]
Tbnd1 = 863.15  # [K   ]

###########################################################
#####       Constants and Auxiliary Functions         #####
###########################################################

# Region 2 and 3 boundary constants and domain scale factors
_n = [ 0.34805185628969e+03, -0.11671859879975e+01,  0.10192970039326e-02,
       0.57254459862746e+03,  0.13918839778870e+02]
_Ps = 1.0 # [MPa]
_Ts = 1.0 # [K  ]

# Auxiliary functions
def bnd23P(T: float) -> float:
    """Auxiliary equation for region 2 and 3 boundary pressure [MPa].
    Reference: Equation (5) from R7-97(2012)"""
    theta = T / _Ts
    return _Ps * (_n[0] + _n[1] * theta + _n[2] * theta**2)

def bnd23T(P: float) -> float:
    """Auxiliary equation for region 2 and 3 boundary temperature [K].
    Reference: Equation (5) from R7-97(2012)"""
    pi = P / _Ps
    return _Ts * (_n[3] + ((pi - _n[4]) / _n[2])**0.5)
