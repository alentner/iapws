"""
Implements basic equations for IF97 region 2.

The basic equation for this region is a fundamental equation for the specific 
Gibbs free energy g. This equation is expressed in dimensionless form, 
gamma = g / RT, and reads as follows, where pi = P / Ps and tau = Ts / T.

    g / RT = gamma(pi, tau) = ln(pi) + SUM( n0i tau^J0i ) + SUM( ni pi^Ii (tau - 0.5)^Ji )

This module is restricted to providing the dimensionless specific Gibbs free
energy, gamma(pi, tau), and the associated first and second derivatives with
respect to pi and tau in region 2. The principle reference for this module
is the IAPWS Industrial Formulation 1997:
    
    International Association for the Properties of Water and Steam, IAPWS R7-97(2012),
    Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic 
    Properties of Water and Steam (2012), available from: http://www.iapws.org.

The dimensional forward and backward functions and their derivatives, 
f(P, T) and f(P, [h, s]) respectively, are defined elsewhere in this package.
"""

# type annotations
from __future__ import annotations

# system libraries
from math import log

###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# Region 2 fitting constants and domain scale factors
_J0 = [  0,  1, -5, -4, -3, -2, -1,  2,  3]
_n0 = [-0.96927686500217e+01,  0.10086655968018e+02, -0.56087911283020e-02,  0.71452738081455e-01,
      -0.40710498223928e+00,  0.14240819171444e+01, -0.43839511319450e+01, -0.28408632460772e+00,
       0.21268463753307e-01]
_Ir = [  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  5,  6,  6,  6,
        7,  7,  7,  8,  8,  9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
_Jr = [  0,  1,  2,  3,  6,  1,  2,  4,  7, 36,  0,  1,  3,  6, 35,  1,  2,  3,  7,  3, 16, 35,
        0, 11, 25,  8, 36, 13,  4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
_nr = [-0.17731742473213e-02, -0.17834862292358e-01, -0.45996013696365e-01, -0.57581259083432e-01,
      -0.50325278727930e-01, -0.33032641670203e-04, -0.18948987516315e-03, -0.39392777243355e-02,
      -0.43797295650573e-01, -0.26674547914087e-04,  0.20481737692309e-07,  0.43870667284435e-06, 
      -0.32277677238570e-04, -0.15033924542148e-02, -0.40668253562649e-01, -0.78847309559367e-09,
       0.12790717852285e-07,  0.48225372718507e-06,  0.22922076337661e-05, -0.16714766451061e-10,
      -0.21171472321355e-02, -0.23895741934104e+02, -0.59059564324270e-17, -0.12621808899101e-05,
      -0.38946842435739e-01,  0.11256211360459e-10, -0.82311340897998e+01,  0.19809712802088e-07,
       0.10406965210174e-18, -0.10234747095929e-12, -0.10018179379511e-08, -0.80882908646985e-10,
       0.10693031879409e+00, -0.33662250574171e+00,  0.89185845355421e-24,  0.30629316876232e-12,
      -0.42002467698208e-05, -0.59056029685639e-25,  0.37826947613457e-05, -0.12768608934681e-14,
       0.73087610595061e-28,  0.55414715350778e-16, -0.94369707241210e-06]
Ps = 1.0      # [Mpa      ]
Ts = 540.0    # [K        ]
R  = 0.461526 # [kJ / kg K]

# Region 2 dimensionless functions
def gamma(pi: float, tau: float) -> float:
    """Dimensionless form for the specific Gibbs free energy.
    Reference: Equation (16) and (17) from R7-97(2012)"""
    sum = log(pi)
    for Ji, ni in zip(_J0, _n0):
        sum += ni * tau**Ji
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * pi**Ii * (tau - 0.5)**Ji
    return sum

def gamma_pi(pi: float, tau: float) -> float:
    """Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi).
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 1 / pi
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gamma_pipi(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi).
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = -1 / pi**2
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gamma_tau(pi: float, tau: float) -> float:
    """Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless temperature (tau).
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0 
    for Ji, ni in zip(_J0, _n0):
        sum += ni * Ji * tau**(Ji - 1)
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * pi**Ii * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gamma_tautau(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless temperature (tau)
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0 
    for Ji, ni in zip(_J0, _n0):
        sum += ni * Ji * (Ji - 1) * tau**(Ji - 2)
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * pi**Ii * Ji * (Ji - 1) * (tau - 0.5)**(Ji - 2)
    return sum

def gamma_pitau(pi: float, tau: float) -> float:
    """Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi) and temperature (tau).
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gammaR_pi(pi: float, tau: float) -> float:
    """ Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi); residual part.
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gammaR_pipi(pi: float, tau: float) -> float:
    """ Derivative (second) of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi); residual part.
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gammaR_pitau(pi: float, tau: float) -> float:
    """ Derivative of dimensionless specific Gibbs free energy,
    w.r.t. dimensionless pressure (pi); residual part.
    Reference: Table (13) and (14) from R7-97(2012)"""
    sum = 0
    for Ii, Ji, ni in zip(_Ir, _Jr, _nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum
