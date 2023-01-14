from ..if97 import h2o
from scipy.optimize import brentq as root
from matplotlib import pyplot

def _dGdP(P, h0, s0):
    h = h2o.h_s(P, s0)
    if (h >= h0):
        return 100.0
    dhdp = h2o.dhdP_s(P, s0)
    v = h2o.v_s(P, s0)
    dvdp = h2o.dvdP_s(P, s0)

    return (dvdp / v**2) * (2 * (h0 - h))**0.5 + dhdp / (v * (2 * (h0 - h))**0.5)
def Pcritical(P0, h0):
    """Critical Pressure; [MPa]
    assuming homogenious equilibrium model"""
    s0 = h2o.s_h(P0, h0)

    return root(_dGdP, 611.213e-6, P0 - 1e-6, (h0, s0))
def Gcritical(P0, h0):
    """Critical mass flux; [kg / m**2 s]
    assuming homogenious equilibrium model"""   
    Pc = Pcritical(P0, h0)
    s0 = h2o.s_h(P0, h0)
    h = h2o.h_s(Pc, s0)
    v = h2o.v_s(Pc, s0)
    
    return 1000**0.5 * (2 * (h0 - h))**0.5 / v
def dGdPcritical(P0, h0):
    """Derivative of critical mass flux [kg / m**2 s]
    w.r.t. stagnation pressure at constant stagnation enthalpy;
    assuming homogenious equilibrium model"""   
    Pc = Pcritical(P0, h0)
    s0 = h2o.s_h(P0, h0)

    h = h2o.h_s(Pc, s0)
    v = h2o.v_s(Pc, s0)

    dhdP0 = h2o.dhds_s(Pc, s0) * h2o.dsdP_h(P0, h0)
    dvdP0 = h2o.dvds_s(Pc, s0) * h2o.dsdP_h(P0, h0)
    
    return -1000**0.5 * (dvdP0 / v**2) * (2 * (h0 - h))**0.5 + dhdP0 / (v * (2 * (h0 - h))**0.5)
