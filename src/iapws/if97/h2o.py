""" Provides the public interface to the water properties """

# type annotations
from __future__ import annotations

# internal libraries
from . import region1, region2, region3, region4
from . import identify
from . import unit
from .support import _english, _first, _one, _region, _zero

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################

#### water properties ####
@_english
@_region
@_formulation
def g(P: float, T: float, *, h: float, s: float,
      english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg]"""
    pass

@_english((unit.v, ))
@_region({1: region1.v, 2: region2.v}, identify.region)
def v(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.T), (unit.u, ))
@_region({1: region1.u, 2: region2.u}, identify.region)
def u(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.T), (unit.s, ))
@_region({1: region1.s, 2: region2.s}, identify.region)
def s(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.h, ))
@_region({1: region1.h, 2: region2.h}, identify.region)
def h(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific enthalpy [kJ / kg]"""
    pass

@_english((unit.P, unit.T), (unit.cp, ))
@_region({1: region1.cp, 2: region2.cp}, identify.region)
def cp(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.cv, ))
@_region({1: region1.cv, 2: region2.cv}, identify.region)
def cv(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.T), (unit.w, ))
@_region({1: region1.w, 2: region2.w}, identify.region)
def w(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.T), (unit.av, ))
@_region({1: region1.av, 2: region2.av}, identify.region)
def av(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.T), (unit.kT, ))
@_region({1: region1.kT, 2: region2.kT}, identify.region)
def kT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

@_english((unit.P, unit.T), (unit.f, ))
@_region({1: region1.f, 2: region2.f}, identify.region)
def f(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific helmholtz free energy [kJ / kg K]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.T), (unit.dP, unit.d_dP))
@_region({1: _one, 2: _one}, identify.region)
def dPdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of pressure [kJ m^3 / m^3 kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dT, unit.d_dP))
@_region({1: _zero, 2: _zero}, identify.region)
def dTdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of temperature [K m^3 / kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP, 2: region2.dvdP}, identify.region)
def dvdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.du, unit.d_dP))
@_region({1: region1.dudP, 2: region2.dudP}, identify.region)
def dudP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dh, unit.d_dP))
@_region({1: region1.dhdP, 2: region2.dhdP}, identify.region)
def dhdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.ds, unit.d_dP))
@_region({1: region1.dsdP, 2: region2.dsdP}, identify.region)
def dsdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP, 2: region2.dgdP}, identify.region)
def dgdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.df, unit.d_dP))
@_region({1: region1.dfdP, 2: region2.dfdP}, identify.region)
def dfdP(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific helmholtz free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dP, unit.d_dT))
@_region({1: _zero, 2: _zero}, identify.region)
def dPdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of pressure [kJ / m^3 K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.dT, unit.d_dT))
@_region({1: _one, 2: _one}, identify.region)
def dTdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of temperature [K / K]
    w.r.t pressure at constant temperature"""
    pass

@_english((unit.P, unit.T), (unit.dv, unit.d_dT))
@_region({1: region1.dvdT, 2: region2.dvdT}, identify.region)
def dvdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.du, unit.d_dT))
@_region({1: region1.dudT, 2: region2.dudT}, identify.region)
def dudT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.ds, unit.d_dT))
@_region({1: region1.dsdT, 2: region2.dsdT}, identify.region)
def dsdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.dh, unit.d_dT))
@_region({1: region1.dhdT, 2: region2.dhdT}, identify.region)
def dhdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.dg, unit.d_dT))
@_region({1: region1.dgdT, 2: region2.dgdT}, identify.region)
def dgdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    pass

@_english((unit.P, unit.T), (unit.df, unit.d_dT))
@_region({1: region1.dfdT, 2: region2.dfdT}, identify.region)
def dfdT(P: float, T: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of specific helmholtz free energy [kJ / kg]
    w.r.t temperature at constant pressure"""
    pass

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################

#### water properties ####
@_english((unit.P, unit.h), (unit.g, ))
@_region({1: region1.g_h, 2: region2.g_h, 4: region4.g_h}, identify.region_h)
def g_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg]"""
    pass

@_english((unit.P, unit.h), (unit.v, ))
@_region({1: region1.v_h, 2: region2.v_h, 4: region4.v_h}, identify.region_h)
def v_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.h), (unit.u, ))
@_region({1: region1.u_h, 2: region2.u_h, 4: region4.u_h}, identify.region_h)
def u_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.h), (unit.s, ))
@_region({1: region1.s_h, 2: region2.s_h, 4: region4.s_h}, identify.region_h)
def s_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.T, ))
@_region({1: region1.T_h, 2: region2.T_h, 4: _first(region4.satT)}, identify.region_h)
def T_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Temperature [K]"""
    pass

@_english((unit.P, unit.h), (unit.cp, ))
@_region({1: region1.cp_h, 2: region2.cp_h, 4: region4.cp_h}, identify.region_h)
def cp_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.cv, ))
@_region({1: region1.cv_h, 2: region2.cv_h, 4: region4.cv_h}, identify.region_h)
def cv_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.h), (unit.w, ))
@_region({1: region1.w_h, 2: region2.w_h, 4: region4.w_h}, identify.region_h)
def w_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.h), (unit.av, ))
@_region({1: region1.av_h, 2: region2.av_h, 4: region4.av_h}, identify.region_h)
def av_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.h), (unit.kT, ))
@_region({1: region1.kT_h, 2: region2.kT_h, 4: region4.kT_h}, identify.region_h)
def kT_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

@_english((unit.P, unit.h), (unit.f, ))
@_region({1: region1.f_h, 2: region2.f_h, 4: region4.f_h}, identify.region_h)
def f_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific helmholtz free energy [kJ / kg]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.T), (unit.dP, unit.d_dP))
@_region({1: _one, 2: _one, 4: _one}, identify.region_h)
def dPdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of pressure [kJ m^3 / m^3 kJ]
    w.r.t pressure at constant enthalpy"""
    pass

@_english((unit.P, unit.T), (unit.dT, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: _first(region4.dTsdP)}, identify.region_h)
def dTdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of temperature [K m^3 / kJ]
    w.r.t pressure at constant enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP_h, 2: region2.dvdP_h, 4: region4.dvdP_h}, identify.region_h)
def dvdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.du, unit.d_dP))
@_region({1: region1.dudP_h, 2: region2.dudP_h, 4: region4.dudP_h}, identify.region_h)
def dudP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dh, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: region4.dhdP_h}, identify.region_h)
def dhdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.ds, unit.d_dP))
@_region({1: region1.dsdP_h, 2: region2.dsdP_h, 4: region4.dsdP_h}, identify.region_h)
def dsdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP_h, 2: region2.dgdP_h, 4: region4.dgdP_h}, identify.region_h)
def dgdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.h), (unit.df, unit.d_dP))
@_region({1: region1.dfdP_h, 2: region2.dfdP_h, 4: region4.dfdP_h}, identify.region_h)
def dfdP_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific helmholtz free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    pass

@_english((unit.P, unit.T), (unit.dP, unit.d_dh))
@_region({1: _zero, 2: _zero, 4: _zero}, identify.region_h)
def dPdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """Derivative of pressure [kJ / m^3 K]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dT, unit.d_dh))
@_region({1: region1.dTdh_h, 2: region2.dTdh_h, 4: _zero}, identify.region_h)
def dTdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of Temperature [K kg / kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dv, unit.d_dh))
@_region({1: region1.dvdh_h, 2: region2.dvdh_h, 4: region4.dvdh_h}, identify.region_h)
def dvdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.du, unit.d_dh))
@_region({1: region1.dudh_h, 2: region2.dudh_h, 4: region4.dudh_h}, identify.region_h)
def dudh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dh, unit.d_dh))
@_region({1: _one, 2: _one, 4: _one}, identify.region_h)
def dhdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.ds, unit.d_dh))
@_region({1: region1.dsdh_h, 2: region2.dsdh_h, 4: region4.dsdh_h}, identify.region_h)
def dsdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ kg / kg K kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.dg, unit.d_dh))
@_region({1: region1.dgdh_h, 2: region2.dgdh_h, 4: region4.dgdh_h}, identify.region_h)
def dgdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.h), (unit.df, unit.d_dh))
@_region({1: region1.dfdh_h, 2: region2.dfdh_h, 4: region4.dfdh_h}, identify.region_h)
def dfdh_h(P: float, h: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific helmholtz free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    pass

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################

#### water properties ####
@_english((unit.P, unit.s), (unit.g, ))
@_region({1: region1.g_s, 2: region2.g_s, 4: region4.g_s}, identify.region_s)
def g_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific gibbs free energy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.v, ))
@_region({1: region1.v_s, 2: region2.v_s, 4: region4.v_s}, identify.region_s)
def v_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific volume [m^3 / kg]"""
    pass

@_english((unit.P, unit.s), (unit.u, ))
@_region({1: region1.u_s, 2: region2.u_s, 4: region4.u_s}, identify.region_s)
def u_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific internal energy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.T, ))
@_region({1: region1.T_s, 2: region2.T_s, 4: _first(region4.satT)}, identify.region_s)
def T_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Temperature [K]"""
    pass

@_english((unit.P, unit.s), (unit.h, ))
@_region({1: region1.h_s, 2: region2.h_s, 4: region4.h_s}, identify.region_s)
def h_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific entropy [kJ / kg]"""
    pass

@_english((unit.P, unit.s), (unit.cp, ))
@_region({1: region1.cp_s, 2: region2.cp_s, 4: region4.cp_s}, identify.region_s)
def cp_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.s), (unit.cv, ))
@_region({1: region1.cv_s, 2: region2.cv_s, 4: region4.cv_s}, identify.region_s)
def cv_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pass

@_english((unit.P, unit.s), (unit.w, ))
@_region({1: region1.w_s, 2: region2.w_s, 4: region4.w_s}, identify.region_s)
def w_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Speed of sound [m / s]"""
    pass

@_english((unit.P, unit.s), (unit.av, ))
@_region({1: region1.av_s, 2: region2.av_s, 4: region4.av_s}, identify.region_s)
def av_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isobaric cubic expansion coefficient [1 / K]"""
    pass

@_english((unit.P, unit.s), (unit.kT, ))
@_region({1: region1.kT_s, 2: region2.kT_s, 4: region4.kT_s}, identify.region_s)
def kT_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Isothermal compressibility [kg / kJ]"""
    pass

@_english((unit.P, unit.s), (unit.f, ))
@_region({1: region1.f_s, 2: region2.f_s, 4: region4.f_s}, identify.region_s)
def f_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """Specific helmholtz free energy [kJ / kg]"""
    pass

#### water property derivatives ####
@_english((unit.P, unit.s), (unit.dP, unit.d_dP))
@_region({1: _one, 2: _one, 4: _one}, identify.region_s)
def dPdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of pressure [kJ m^3 /m^3 kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: _first(region4.dTsdP)}, identify.region_s)
def dTdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of temperature [K m^3 / kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.dv, unit.d_dP))
@_region({1: region1.dvdP_s, 2: region2.dvdP_s, 4: region4.dvdP_s}, identify.region_s)
def dvdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.du, unit.d_dP))
@_region({1: region1.dudP_s, 2: region2.dudP_s, 4: region4.dudP_s}, identify.region_s)
def dudP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.dh, unit.d_dP))
@_region({1: region1.dhdP_s, 2: region2.dhdP_s, 4: region4.dhdP_s}, identify.region_s)
def dhdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_dP))
@_region({1: _zero, 2: _zero, 4: region4.dsdP_s}, identify.region_s)
def dsdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific/equilibrium entropy"""
    pass

@_english((unit.P, unit.s), (unit.dg, unit.d_dP))
@_region({1: region1.dgdP_s, 2: region2.dgdP_s, 4: region4.dgdP_s}, identify.region_s)
def dgdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.df, unit.d_dP))
@_region({1: region1.dfdP_s, 2: region2.dfdP_s, 4: region4.dfdP_s}, identify.region_s)
def dfdP_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific helmholtz free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    pass

@_english((unit.P, unit.s), (unit.dP, unit.d_ds))
@_region({1: _zero, 2: _zero, 4: _zero}, identify.region_s)
def dPds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of pressure [kJ kg K /m^3 kJ]
    w.r.t enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dT, unit.d_ds))
@_region({1: region1.dTds_s, 2: region2.dTds_s, 4: _zero}, identify.region_s)
def dTds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of temperature [K kg K / kJ]
    w.r.t enthalpy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dv, unit.d_ds))
@_region({1: region1.dvds_s, 2: region2.dvds_s, 4: region4.dvds_s}, identify.region_s)
def dvds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific volume [m^3 kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.du, unit.d_ds))
@_region({1: region1.duds_s, 2: region2.duds_s, 4: region4.duds_s}, identify.region_s)
def duds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific internal energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dh, unit.d_ds))
@_region({1: region1.dhds_s, 2: region2.dhds_s, 4: region4.dhds_s}, identify.region_s)
def dhds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific enthalpy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.ds, unit.d_ds))
@_region({1: _one, 2: _one, 4: _one}, identify.region_s)
def dsds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific entropy [kJ kg K / kg K kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.dg, unit.d_ds))
@_region({1: region1.dgds_s, 2: region2.dgds_s, 4: region4.dgds_s}, identify.region_s)
def dgds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

@_english((unit.P, unit.s), (unit.df, unit.d_ds))
@_region({1: region1.dfds_s, 2: region2.dfds_s, 4: region4.dfds_s}, identify.region_s)
def dfds_s(P: float, s: float, /, *, english: bool = False, region: int = 0) -> float:
    """ Derivative of specific helmholtz free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    pass

###########################################################
#####     Pressure Only (Saturation) Formulation      #####
###########################################################

#### P-T saturation curves ####
@_english((unit.T, ), (unit.P, ))
def satP(T: float, /, *, english: bool = False) -> float:
    """Saturation Pressure [Mpa] for specified Temperature"""
    return region4.satP(T)

@_english((unit.P, ), (unit.T, ))
def satT(P: float, /, *, english: bool = False) -> float:
    """Saturation Temperature [K] for specified Pressure"""
    return region4.satT(P)

#### Saturated liquid properties ####
@_english((unit.P, ), (unit.g, ))
def gf(P: float, /, *, english: bool = False) -> float:
    """Specific gibbs free energy [kJ / kg] of saturated liquid"""
    return region4.gf(P)

@_english((unit.P, ), (unit.v, ))
def vf(P: float, /, *, english: bool = False) -> float:
    """Specific volume [m^3 / kg] of saturated liquid"""
    return region4.vf(P)

@_english((unit.P, ), (unit.u, ))
def uf(P: float, /, *, english: bool = False) -> float:
    """Specific internal energy [kJ / kg] of saturated liquid"""
    return region4.uf(P)

@_english((unit.P, ), (unit.s, ))
def sf(P: float, /, *, english: bool = False) -> float:
    """Specific entropy [kJ / kg K] of saturated liquid"""
    return region4.sf(P)

@_english((unit.P, ), (unit.h, ))
def hf(P: float, /, *, english: bool = False) -> float:
    """Specific enthalpy [kJ / kg] of saturated liquid"""
    return region4.hf(P)

@_english((unit.P, ), (unit.cp, ))
def cpf(P: float, /, *, english: bool = False) -> float:
    """Specific isobaric heat capacity [kJ / kg K] of saturated liquid"""
    return region4.cpf(P)

@_english((unit.P, ), (unit.cv, ))
def cvf(P: float, /, *, english: bool = False) -> float:
    """Specific isochoric heat capacity [kJ / kg K] of saturated liquid"""
    return region4.cvf(P)

@_english((unit.P, ), (unit.w, ))
def wf(P: float, /, *, english: bool = False) -> float:
    """Speed of sound [m / s] of saturated liquid"""
    return region4.wf(P)

@_english((unit.P, ), (unit.av, ))
def avf(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient [1 / K] of saturated liquid"""
    return region4.af(P)

@_english((unit.P, ), (unit.kT, ))
def kTf(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility [kg / kJ] of saturated liquid"""
    return region4.kf(P)

@_english((unit.P, ), (unit.f, ))
def ff(P: float, /, *, english: bool = False) -> float:
    """Specific helmholtz free energy [kJ / kg] of saturated liquid"""
    return region4.ff(P)

#### Saturated vapor properties ####
@_english((unit.P, ), (unit.g, ))
def gg(P: float, /, *, english: bool = False) -> float:
    """Specific gibbs free energy [kJ / kg] of saturated vapor"""
    return region4.gg(P)

@_english((unit.P, ), (unit.v, ))
def vg(P: float, /, *, english: bool = False) -> float:
    """Specific volume [m^3 / kg] of saturated vapor"""
    return region4.vg(P)

@_english((unit.P, ), (unit.u, ))
def ug(P: float, /, *, english: bool = False) -> float:
    """Specific internal energy [kJ / kg] of saturated vapor"""
    return region4.ug(P)

@_english((unit.P, ), (unit.s, ))
def sg(P: float, /, *, english: bool = False) -> float:
    """Specific entropy [kJ / kg K] of saturated vapor"""
    return region4.sg(P)

@_english((unit.P, ), (unit.h, ))
def hg(P: float, /, *, english: bool = False) -> float:
    """Specific enthalpy [kJ / kg] of saturated vapor"""
    return region4.hg(P)

@_english((unit.P, ), (unit.cp, ))
def cpg(P: float, /, *, english: bool = False) -> float:
    """Specific isobaric heat capacity [kJ / kg K] of saturated vapor"""
    return region4.cpg(P)

@_english((unit.P, ), (unit.cv, ))
def cvg(P: float, /, *, english: bool = False) -> float:
    """Specific isochoric heat capacity [kJ / kg K] of saturated vapor"""
    return region4.cvg(P)

@_english((unit.P, ), (unit.w, ))
def wg(P: float, /, *, english: bool = False) -> float:
    """Speed of sound [m / s] of saturated vapor"""
    return region4.wg(P)

@_english((unit.P, ), (unit.av, ))
def avg(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient [1 / K] of saturated vapor"""
    return region4.ag(P)

@_english((unit.P, ), (unit.kT, ))
def kTg(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility [kg / kJ] of saturated vapor"""
    return region4.kg(P)

@_english((unit.P, ), (unit.f, ))
def fg(P: float, /, *, english: bool = False) -> float:
    """Specific helmholtz free energy [kJ / kg] of saturated vapor"""
    return region4.fg(P)

#### delta saturation properties ####
@_english((unit.P, ), (unit.g, ))
def gfg(P: float, /, *, english: bool = False) -> float:
    """Specific gibbs free energy; [kJ / kg] saturation rise of"""
    return region4.gfg(P)

@_english((unit.P, ), (unit.v, ))
def vfg(P: float, /, *, english: bool = False) -> float:
    """Specific volume; [m^3 / kg] saturation rise of"""
    return region4.vfg(P)

@_english((unit.P, ), (unit.u, ))
def ufg(P: float, /, *, english: bool = False) -> float:
    """Specific internal energy; [kJ / kg] saturation rise of"""
    return region4.ufg(P)

@_english((unit.P, ), (unit.s, ))
def sfg(P: float, /, *, english: bool = False) -> float:
    """Specific entropy; [kJ / kg K] saturation rise of"""
    return region4.sfg(P)

@_english((unit.P, ), (unit.h, ))
def hfg(P: float, /, *, english: bool = False) -> float:
    """Specific enthalpy; [kJ / kg] saturation rise of"""
    return region4.hfg(P)

@_english((unit.P, ), (unit.cp, ))
def cpfg(P: float, /, *, english: bool = False) -> float:
    """Specific isobaric heat capacity; [kJ / kg K] saturation rise of"""
    return region4.cpfg(P)

@_english((unit.P, ), (unit.cv, ))
def cvfg(P: float, /, *, english: bool = False) -> float:
    """Specific isochoric heat capacity; [kJ / kg K] saturation rise of"""
    return region4.cvfg(P)

@_english((unit.P, ), (unit.w, ))
def wfg(P: float, /, *, english: bool = False) -> float:
    """Speed of sound; [m / s] saturation rise of"""
    return region4.wfg(P)

@_english((unit.P, ), (unit.av, ))
def avfg(P: float, /, *, english: bool = False) -> float:
    """Isobaric cubic expansion coefficient; [1 / K] saturation rise of"""
    return region4.afg(P)

@_english((unit.P, ), (unit.kT, ))
def kTfg(P: float, /, *, english: bool = False) -> float:
    """Isothermal compressibility; [kg / kJ] saturation rise of"""
    return region4.kfg(P)

@_english((unit.P, ), (unit.f, ))
def ffg(P: float, /, *, english: bool = False) -> float:
    """Specific helmholtz free energy; [kJ / kg] saturation rise of"""
    return region4.ffg(P)

#### Saturated liquid derivatives ####
@_english((unit.P, ), (unit.dv, unit.d_dP))
def dvfdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dvfdP(P)

@_english((unit.P, ), (unit.du, unit.d_dP))
def dufdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dufdP(P)

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhfdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dhfdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsfdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dsfdP(P)

@_english((unit.P, ), (unit.dg, unit.d_dP))
def dgfdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dgfdP(P)

@_english((unit.P, ), (unit.df, unit.d_dP))
def dffdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific helmholtz free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    return region4.dffdP(P)

#### Saturated vapor derivatives ####
@_english((unit.P, ), (unit.dv, unit.d_dP))
def dvgdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dvgdP(P)

@_english((unit.P, ), (unit.du, unit.d_dP))
def dugdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dugdP(P)

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dhgdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsgdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dsgdP(P)

@_english((unit.P, ), (unit.dg, unit.d_dP))
def dggdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dggdP(P)

@_english((unit.P, ), (unit.df, unit.d_dP))
def dfgdP(P: float, /, *, english: bool = False) -> float:
    """Derivative of Specific helmholtz free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    return region4.dfgdP(P)

#### Delta saturation derivatives ####
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

@_english((unit.P, ), (unit.dh, unit.d_dP))
def dhfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dhfgdP(P)

@_english((unit.P, ), (unit.ds, unit.d_dP))
def dsfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dsfgdP(P)

@_english((unit.P, ), (unit.dg, unit.d_dP))
def dgfgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dgfgdP(P)

@_english((unit.P, ), (unit.df, unit.d_dP))
def dffgdP(P: float, /, *, english: bool = False) -> float:
    """ Derivative of Specific helmholtz free energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""
    return region4.dffgdP(P)
