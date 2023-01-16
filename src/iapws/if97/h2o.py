""" Provides the public interface to the water properties """

# internal libraries
from . import region1, region2, region3, region4
from . import unit
from .support import _english, _first, _one, _output, _region, _zero

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################
@_english((unit.P, unit.T), (_output, ))
def idRegion(P: float, T: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and temperature as primary varibles"""

    # Constant boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1
    Tbnd01 = region1.Tbnd01
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13

    # Con-constant boundaries
    Pbnd32 = region3.bnd23P(min(max(T, Tbnd13), 863.15))
    Pbnd4  = satP(min(max(T, Tbnd01), Tbnd13))

    region = 0

    # Only region 1, and 2 via P,T relations implemented
    if (P >= Pbnd0) and (T >= Tbnd01) and (P <= Pbnd1) and (T <= Tbnd25):
        if (T <= Tbnd13) and (P >= Pbnd4):
            region = 1
        elif (T < Tbnd13) or (P <= Pbnd32):
            region = 2
        else:
            region = 0

    assert (region != 0), "Water properties not avalable!"
    return region

#### water properties ####
@_english((unit.P, unit.T), (unit.g, ))
@_region({1: region1.g, 2: region2.g}, idRegion)
def g(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.v, ))
@_region({1: region1.v, 2: region2.v}, idRegion)
def v(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.T), (unit.u, ))
@_region({1: region1.u, 2: region2.u}, idRegion)
def u(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.T), (unit.s, ))
@_region({1: region1.s, 2: region2.s}, idRegion)
def s(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.h, ))
@_region({1: region1.h, 2: region2.h}, idRegion)
def h(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific enthalpy [kJ / kg]"""
    pass

@_english((unit.P, unit.T), (unit.cp, ))
@_region({1: region1.cp, 2: region2.cp}, idRegion)
def cp(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.cv, ))
@_region({1: region1.cv, 2: region2.cv}, idRegion)
def cv(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.w, ))
@_region({1: region1.w, 2: region2.w}, idRegion)
def w(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.T), (unit.a, ))
@_region({1: region1.a, 2: region2.a}, idRegion)
def a(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.T), (unit.k, ))
@_region({1: region1.k, 2: region2.k}, idRegion)
def k(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.T), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP, 2: region2.dgdP}, idRegion)
def dgdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP, 2: region2.dvdP}, idRegion)
def dvdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.du, unit.d_dP))
@_region({1: region1.dudP, 2: region2.dudP}, idRegion)
def dudP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.ds, unit.d_dP))
@_region({1: region1.dsdP, 2: region2.dsdP}, idRegion)
def dsdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dh, unit.d_dP))
@_region({1: region1.dhdP, 2: region2.dhdP}, idRegion)
def dhdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dg, unit.d_dT))
@_region({1: region1.dgdT, 2: region2.dgdT}, idRegion)
def dgdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.dv, unit.d_dT))
@_region({1: region1.dvdT, 2: region2.dvdT}, idRegion)
def dvdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.du, unit.d_dT))
@_region({1: region1.dudT, 2: region2.dudT}, idRegion)
def dudT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.ds, unit.d_dT))
@_region({1: region1.dsdT, 2: region2.dsdT}, idRegion)
def dsdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.dh, unit.d_dT))
@_region({1: region1.dhdT, 2: region2.dhdT}, idRegion)
def dhdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################
@_english((unit.P, unit.h), (_output, ))
def idRegion_h(P: float, h: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # Supporting boundaries
    Tbnd01 = region1.Tbnd01
    Pbnd4  = satP(Tbnd01)
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13
    Tbnd32 = region3.bnd23T(min(max(P, 16.5292), 100.0))
    Tbnd4  = satT(P)

    # Enthalpy- pressure boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1 
    hbnd01 = region1.h(Pbnd4, Tbnd01)
    hbnd25 = region2.h(Pbnd0, Tbnd25)
    Pbndh1 = satP(Tbnd13)
    hbnd13 = region1.h(P, Tbnd13)
    hbnd32 = region2.h(P, Tbnd32)
    hbnd14 = region1.h(P, Tbnd4)
    hbnd42 = region2.h(P, Tbnd4)

    region = 0

    # Only region 1,2,4 via P,h relations implemented
    if (P >= Pbnd0) and (h >= hbnd01) and (P <= Pbnd1) and (h <= hbnd25):
        if (P >= Pbndh1):
            if (h <= hbnd13):
                region = 1
            elif (h >= hbnd32):
                region = 2
            else:
                region = 0
        else:
            if (h <= hbnd14):
                region = 1
            elif (h >= hbnd42):
                region = 2
            else:
                region = 4

    assert (region != 0), "Water properties not avalable!"
    return region

#### water properties ####
@_english((unit.P, unit.h), (unit.g, ))
@_region({1: region1.g_h, 2: region2.g_h, 4: region4.g_h}, idRegion_h)
def g_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg]"""
    pass

@_english((unit.P, unit.h), (unit.v, ))
@_region({1: region1.v_h, 2: region2.v_h, 4: region4.v_h}, idRegion_h)
def v_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.h), (unit.u, ))
@_region({1: region1.u_h, 2: region2.u_h, 4: region4.u_h}, idRegion_h)
def u_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.h), (unit.s, ))
@_region({1: region1.s_h, 2: region2.s_h, 4: region4.s_h}, idRegion_h)
def s_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.T, ))
@_region({1: region1.T_h, 2: region2.T_h, 4: _first(region4.satT)}, idRegion_h)
def T_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Temperature [K]"""
    pass

@_english((unit.P, unit.h), (unit.cp, ))
@_region({1: region1.cp_h, 2: region2.cp_h, 4: region4.cp_h}, idRegion_h)
def cp_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.cv, ))
@_region({1: region1.cv_h, 2: region2.cv_h, 4: region4.cv_h}, idRegion_h)
def cv_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.w, ))
@_region({1: region1.w_h, 2: region2.w_h, 4: region4.w_h}, idRegion_h)
def w_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.h), (unit.a, ))
@_region({1: region1.a_h, 2: region2.a_h, 4: region4.a_h}, idRegion_h)
def a_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.h), (unit.k, ))
@_region({1: region1.k_h, 2: region2.k_h, 4: region4.k_h}, idRegion_h)
def k_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.h), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP_h, 2: region2.dgdP_h, 4: region4.dgdP_h}, idRegion_h)
def dgdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP_h, 2: region2.dvdP_h, 4: region4.dvdP_h}, idRegion_h)
def dvdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.du, unit.d_dP))
@_region({1: region1.dudP_h, 2: region2.dudP_h, 4: region4.dudP_h}, idRegion_h)
def dudP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.ds, unit.d_dP))
@_region({1: region1.dsdP_h, 2: region2.dsdP_h, 4: region4.dsdP_h}, idRegion_h)
def dsdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dh, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: region4.dhdP_h}, idRegion_h)
def dhdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dT, unit.d_dP))
@_region({1: region1.dTdP_h, 2: region2.dTdP_h, 4: _first(region4.dTsdP)}, idRegion_h)
def dTdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dg, unit.d_dh))
@_region({1: region1.dgdh_h, 2: region2.dgdh_h, 4: region4.dgdh_h}, idRegion_h)
def dgdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dv, unit.d_dh))
@_region({1: region1.dvdh_h, 2: region2.dvdh_h, 4: region4.dvdh_h}, idRegion_h)
def dvdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.du, unit.d_dh))
@_region({1: region1.dudh_h, 2: region2.dudh_h, 4: region4.dudh_h}, idRegion_h)
def dudh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.ds, unit.d_dh))
@_region({1: region1.dsdh_h, 2: region2.dsdh_h, 4: region4.dsdh_h}, idRegion_h)
def dsdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ kg / kg K kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dh, unit.d_dh))
@_region({1: _one, 2: _one, 4: _one}, idRegion_h)
def dhdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dT, unit.d_dh))
@_region({1: region1.dTdh_h, 2: region2.dTdh_h, 4: _zero}, idRegion_h)
def dTdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of Temperature [K kg / kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################
@_english((unit.P, unit.s), (_output, ))
def idRegion_s(P: float, s: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # Supporting boundaries
    Tbnd01 = region1.Tbnd01
    Pbnd4  = satP(Tbnd01)
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13
    Tbnd32 = region3.bnd23T(min(max(P, 16.5292), 100.0))
    Tbnd4  = satT(P)

    # Enthalpy-pressure boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1 
    sbnd01 = region1.s(P, Tbnd01)
    sbnd25 = region2.s(P, Tbnd25)
    Pbndh1 = satP(Tbnd13)
    sbnd13 = region1.s(P, Tbnd13)
    sbnd32 = region2.s(P, Tbnd32)
    sbnd14 = region1.s(P, Tbnd4)
    sbnd42 = region2.s(P, Tbnd4)

    region = 0

    # Only region 1,2,4 via P,s relations implemented
    if (P >= Pbnd0) and (s >= sbnd01) and (P <= Pbnd1) and (s <= sbnd25):
        if (P >= Pbndh1):
            if (s <= sbnd13):
                region = 1
            elif (s >= sbnd32):
                region = 2
            else:
                region = 0
        else:
            if (s <= sbnd14):
                region = 1
            elif (s >= sbnd42):
                region = 2
            else:
                region = 4

    assert (region != 0), "Water properties not avalable!"
    return region

#### water properties ####
@_english((unit.P, unit.s), (unit.g, ))
@_region({1: region1.g_s, 2: region2.g_s, 4: region4.g_s}, idRegion_s)
def g_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.v, ))
@_region({1: region1.v_s, 2: region2.v_s, 4: region4.v_s}, idRegion_s)
def v_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.s), (unit.u, ))
@_region({1: region1.u_s, 2: region2.u_s, 4: region4.u_s}, idRegion_s)
def u_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.T, ))
@_region({1: region1.T_s, 2: region2.T_s, 4: _first(region4.satT)}, idRegion_s)
def T_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Temperature [K]"""
    pass

@_english((unit.P, unit.s), (unit.h, ))
@_region({1: region1.h_s, 2: region2.h_s, 4: region4.h_s}, idRegion_s)
def h_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.cp, ))
@_region({1: region1.cp_s, 2: region2.cp_s, 4: region4.cp_s}, idRegion_s)
def cp_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.s), (unit.cv, ))
@_region({1: region1.cv_s, 2: region2.cv_s, 4: region4.cv_s}, idRegion_s)
def cv_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.s), (unit.w, ))
@_region({1: region1.w_s, 2: region2.w_s, 4: region4.w_s}, idRegion_s)
def w_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.s), (unit.a, ))
@_region({1: region1.a_s, 2: region2.a_s, 4: region4.a_s}, idRegion_s)
def a_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.s), (unit.k, ))
@_region({1: region1.k_s, 2: region2.k_s, 4: region4.k_s}, idRegion_s)
def k_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.s), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP_s, 2: region2.dgdP_s, 4: region4.dgdP_s}, idRegion_s)
def dgdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP_s, 2: region2.dvdP_s, 4: region4.dvdP_s}, idRegion_s)
def dvdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.du, unit.d_dP))
@_region({1: region1.dudP_s, 2: region2.dudP_s, 4: region4.dudP_s}, idRegion_s)
def dudP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: region4.dsdP_s}, idRegion_s)
def dsdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific/equilibrium entropy"""
    pass

@_english((unit.P, unit.s), (unit.dh, unit.d_dP))
@_region({1: region1.dhdP_s, 2: region2.dhdP_s, 4: region4.dhdP_s}, idRegion_s)
def dhdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_dP))
@_region({1: region1.dTdP_s, 2: region2.dTdP_s, 4: _first(region4.dTsdP)}, idRegion_s)
def dTdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.dg, unit.d_ds))
@_region({1: region1.dgds_s, 2: region2.dgds_s, 4: region4.dgds_s}, idRegion_s)
def dgds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dv, unit.d_ds))
@_region({1: region1.dvds_s, 2: region2.dvds_s, 4: region4.dvds_s}, idRegion_s)
def dvds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.du, unit.d_ds))
@_region({1: region1.duds_s, 2: region2.duds_s, 4: region4.duds_s}, idRegion_s)
def duds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_ds))
@_region({1: _one, 2: _one, 4: _one}, idRegion_s)
def dsds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ kg K / kg K kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dh, unit.d_ds))
@_region({1: region1.dhds_s, 2: region2.dhds_s, 4: region4.dhds_s}, idRegion_s)
def dhds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dT, unit.d_ds))
@_region({1: region1.dTds_s, 2: region2.dTds_s, 4: _zero}, idRegion_s)
def dTds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of Temperature [K kg K / kJ]
    w.r.t enthalpy at constant pressure"""
    pass

###########################################################
#####     Pressure Only (Saturation) Formulation      #####
###########################################################

#### P-T saturation curves ####
@_english((unit.T, ), (unit.P, ))
def satP(T: float, /, *, english: bool = False) -> float:
    """ Saturation Pressure [Mpa]
    for specified Temperature"""
    return region4.satP(T)

@_english((unit.P, ), (unit.T, ))
def satT(P: float, /, *, english: bool = False) -> float:
    """ Saturation Temperature [K]
    for specified Pressure"""
    return region4.satT(P)

#### Saturated liquid properties ####
@_english((unit.P, ), (unit.g, ))
def gf(P: float, /, *, english: bool = False) -> float:
    """ Specific gibbs free energy [kJ / kg]
    of saturated liquid"""
    return region4.gf(P)

@_english((unit.P, ), (unit.v, ))
def vf(P: float, /, *, english: bool = False) -> float:
    """ Specific volume [m^3 / kg]
    of saturated liquid"""
    return region4.vf(P)

@_english((unit.P, ), (unit.u, ))
def uf(P: float, /, *, english: bool = False) -> float:
    """ Specific internal energy [kJ / kg]
    of saturated liquid"""
    return region4.uf(P)

@_english((unit.P, ), (unit.s, ))
def sf(P: float, /, *, english: bool = False) -> float:
    """ Specific entropy [kJ / kg K]
    of saturated liquid"""
    return region4.sf(P)

@_english((unit.P, ), (unit.h, ))
def hf(P: float, /, *, english: bool = False) -> float:
    """ Specific enthalpy [kJ / kg]
    of saturated liquid"""
    return region4.hf(P)

@_english((unit.P, ), (unit.cp, ))
def cpf(P: float, /, *, english: bool = False) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated liquid"""
    return region4.cpf(P)

@_english((unit.P, ), (unit.cv, ))
def cvf(P: float, /, *, english: bool = False) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated liquid"""
    return region4.cvf(P)

@_english((unit.P, ), (unit.w, ))
def wf(P: float, /, *, english: bool = False) -> float:
    """ Speed of sound [m / s]
    of saturated liquid"""
    return region4.wf(P)

@_english((unit.P, ), (unit.a, ))
def af(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated liquid"""
    return region4.af(P)

@_english((unit.P, ), (unit.k, ))
def kf(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility [kg / kJ]
    of saturated liquid"""
    return region4.kf(P)

#### Saturated vapor properties ####
@_english((unit.P, ), (unit.g, ))
def gg(P: float, /, *, english: bool = False) -> float:
    """ Specific gibbs free energy [kJ / kg]
    of saturated vapor"""
    return region4.gg(P)

@_english((unit.P, ), (unit.v, ))
def vg(P: float, /, *, english: bool = False) -> float:
    """ Specific volume [m^3 / kg]
    of saturated vapor"""
    return region4.vg(P)

@_english((unit.P, ), (unit.u, ))
def ug(P: float, /, *, english: bool = False) -> float:
    """ Specific internal energy [kJ / kg]
    of saturated vapor"""
    return region4.ug(P)

@_english((unit.P, ), (unit.s, ))
def sg(P: float, /, *, english: bool = False) -> float:
    """ Specific entropy [kJ / kg K]
    of saturated vapor"""
    return region4.sg(P)

@_english((unit.P, ), (unit.h, ))
def hg(P: float, /, *, english: bool = False) -> float:
    """ Specific enthalpy [kJ / kg]
    of saturated vapor"""
    return region4.hg(P)

@_english((unit.P, ), (unit.cp, ))
def cpg(P: float, /, *, english: bool = False) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated vapor"""
    return region4.cpg(P)

@_english((unit.P, ), (unit.cv, ))
def cvg(P: float, /, *, english: bool = False) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated vapor"""
    return region4.cvg(P)

@_english((unit.P, ), (unit.w, ))
def wg(P: float, /, *, english: bool = False) -> float:
    """ Speed of sound [m / s]
    of saturated vapor"""
    return region4.wg(P)

@_english((unit.P, ), (unit.a, ))
def ag(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated vapor"""
    return region4.ag(P)

@_english((unit.P, ), (unit.k, ))
def kg(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility [kg / kJ]
    of saturated vapor"""
    return region4.kg(P)

#### delta saturation properties ####
@_english((unit.P, ), (unit.g, ))
def gfg(P: float, /, *, english: bool = False) -> float:
    """ Specific gibbs free energy; [kJ / kg]
    saturation rise of"""
    return region4.gfg(P)

@_english((unit.P, ), (unit.v, ))
def vfg(P: float, /, *, english: bool = False) -> float:
    """ Specific volume; [m^3 / kg]
    saturation rise of"""
    return region4.vfg(P)

@_english((unit.P, ), (unit.u, ))
def ufg(P: float, /, *, english: bool = False) -> float:
    """ Specific internal energy; [kJ / kg]
    saturation rise of"""
    return region4.ufg(P)

@_english((unit.P, ), (unit.s, ))
def sfg(P: float, /, *, english: bool = False) -> float:
    """ Specific entropy; [kJ / kg K]
    saturation rise of"""
    return region4.sfg(P)

@_english((unit.P, ), (unit.h, ))
def hfg(P: float, /, *, english: bool = False) -> float:
    """ Specific enthalpy; [kJ / kg]
    saturation rise of"""
    return region4.hfg(P)

@_english((unit.P, ), (unit.cp, ))
def cpfg(P: float, /, *, english: bool = False) -> float:
    """ Specific isobaric heat capacity; [kJ / kg K]
    saturation rise of"""
    return region4.cpfg(P)

@_english((unit.P, ), (unit.cv, ))
def cvfg(P: float, /, *, english: bool = False) -> float:
    """ Specific isochoric heat capacity; [kJ / kg K]
    saturation rise of"""
    return region4.cvfg(P)

@_english((unit.P, ), (unit.w, ))
def wfg(P: float, /, *, english: bool = False) -> float:
    """ Speed of sound; [m / s]
    saturation rise of"""
    return region4.wfg(P)

@_english((unit.P, ), (unit.a, ))
def afg(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient; [1 / K]
    saturation rise of"""
    return region4.afg(P)

@_english((unit.P, ), (unit.k, ))
def kfg(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility; [kg / kJ]
    saturation rise of"""
    return region4.kfg(P)

#### Saturated liquid derivatives ####
@_english((unit.P, ), (unit.dg, unit.d_dP))
def dgfdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dgfdP(P)

@_english((unit.P, ), (unit.dv, unit.d_dP))
def dvfdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dvfdP(P)

@_english((unit.P, ), (unit.du, unit.d_dP))
def dufdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dufdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsfdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dsfdP(P)

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhfdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dhfdP(P)

#### Saturated vapor derivatives ####
@_english((unit.P, ), (unit.dg, unit.d_dP))
def dggdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dggdP(P)

@_english((unit.P, ), (unit.dv, unit.d_dP))
def dvgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dvgdP(P)

@_english((unit.P, ), (unit.du, unit.d_dP))
def dugdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dugdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dsgdP(P)

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dhgdP(P)

#### Delta saturation derivatives ####
@_english((unit.P, ), (unit.dg, unit.d_dP))
def dgfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dgfgdP(P)

@_english((unit.P, ), (unit.dv, unit.d_dP))
def dvfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dvfgdP(P)

@_english((unit.P, ), (unit.du, unit.d_dP))
def dufgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dufgdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dsfgdP(P)

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dhfgdP(P)
