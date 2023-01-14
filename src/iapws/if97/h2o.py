from . import region1, region2, region3, region4

###########################################################
#####          Pressure-Temperature Formulation       #####
###########################################################
def idRegion(P, T):
    """Identification of region from IF97 specification
    using pressure and temperature as primary varibles"""

    # Constant boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1
    Tbnd01 = region1.Tbnd01
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13

    # non-constant boundaries
    Pbnd32 = region3.bnd23P(min(max(T, Tbnd13), 863.15))
    Pbnd4  = satP(min(max(T, Tbnd01), Tbnd13))

    region = 0

    if (P >= Pbnd0) and (T >= Tbnd01) and (P <= Pbnd1) and (T <= Tbnd25):
        if (T <= Tbnd13) and (P >= Pbnd4):
            region = 1
        elif (T < Tbnd13) or (P <= Pbnd32):
            region = 2
        else:
            # region 3 via P,T relations not implemented
            region = 0
    assert (region is not 0), "Water properties not avalable!"
    return region

#### water properties ####
def g(P, T, region = 0):
    """Specific gibbs free energy [kJ / kg K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.g(P, T)
    elif region is 2:
        return region2.g(P, T)
    else:
        return 0.000
def v(P, T, region = 0):
    """Specific volume [m^3 / kg]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.v(P, T)
    elif region is 2:
        return region2.v(P, T)
    else:
        return 0.000
def u(P, T, region = 0):
    """Specific internal energy [kJ / kg]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.u(P, T)
    elif region is 2:
        return region2.u(P, T)
    else:
        return 0.000
def s(P, T, region = 0):
    """Specific entropy [kJ / kg K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.s(P, T)
    elif region is 2:
        return region2.s(P, T)
    else:
        return 0.000
def h(P, T, region = 0):
    """Specific enthalpy [kJ / kg]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.h(P, T)
    elif region is 2:
        return region2.h(P, T)
    else:
        return 0.000
def cp(P, T, region = 0):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.cp(P, T)
    elif region is 2:
        return region2.cp(P, T)
    else:
        return 0.000
def cv(P, T, region = 0):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.cv(P, T)
    elif region is 2:
        return region2.cv(P, T)
    else:
        return 0.000
def w(P, T, region = 0):
    """ Speed of sound [m / s]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.w(P, T)
    elif region is 2:
        return region2.w(P, T)
    else:
        return 0.000
def a(P, T, region = 0):
    """Isobaric cubic expansion coefficient [1 / K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.a(P, T)
    elif region is 2:
        return region2.a(P, T)
    else:
        return 0.000
def k(P, T, region = 0):
    """Isothermal compressibility [kg / kJ]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.k(P, T)
    elif region is 2:
        return region2.k(P, T)
    else:
        return 0.000

#### water property derivatives ####
def dgdP(P, T, region = 0):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dgdP(P, T)
    elif region is 2:
        return region2.dgdP(P, T)
    else:
        return 0.000
def dvdP(P, T, region = 0):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dvdP(P, T)
    elif region is 2:
        return region2.dvdP(P, T)
    else:
        return 0.000
def dudP(P, T, region = 0):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dudP(P, T)
    elif region is 2:
        return region2.dudP(P, T)
    else:
        return 0.000
def dsdP(P, T, region = 0):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant temperature"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dsdP(P, T)
    elif region is 2:
        return region2.dsdP(P, T)
    else:
        return 0.000
def dhdP(P, T, region = 0):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant temperature"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dhdP(P, T)
    elif region is 2:
        return region2.dhdP(P, T)
    else:
        return 0.000

def dgdT(P, T, region = 0):
    """ Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dgdT(P, T)
    elif region is 2:
        return region2.dgdT(P, T)
    else:
        return 0.000
def dvdT(P, T, region = 0):
    """ Derivative of specific volume [m^3 / kg K]
    w.r.t temperature at constant pressure"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dvdT(P, T)
    elif region is 2:
        return region2.dvdT(P, T)
    else:
        return 0.000
def dudT(P, T, region = 0):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dudT(P, T)
    elif region is 2:
        return region2.dudT(P, T)
    else:
        return 0.000
def dsdT(P, T, region = 0):
    """ Derivative of specific entropy [kJ / kg K K]
    w.r.t temperature at constant pressure"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dsdT(P, T)
    elif region is 2:
        return region2.dsdT(P, T)
    else:
        return 0.000
def dhdT(P, T, region = 0):
    """ Derivative of specific enthalpy [kJ / kg K]
    w.r.t temperature at constant pressure"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.dhdT(P, T)
    elif region is 2:
        return region2.dhdT(P, T)
    else:
        return 0.000

###########################################################
#####          Pressure-Enthalpy Formulation          #####
###########################################################
def idRegion_h(P, h):
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # supporting boundaries
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

    if (P >= Pbnd0) and (h >= hbnd01) and (P <= Pbnd1) and (h <= hbnd25):
        if (P >= Pbndh1):
            if (h <= hbnd13):
                region = 1
            elif (h >= hbnd32):
                region = 2
            else:
                # region 3 via P,h relations not implemented
                region = 0
        else:
            if (h <= hbnd14):
                region = 1
            elif (h >= hbnd42):
                region = 2
            else:
                region = 4
    assert (region is not 0), "Water properties not avalable!"
    return region

#### water properties ####
def g_h(P, h, region = 0):
    """Specific gibbs free energy [kJ / kg]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.g_h(P, h)
    elif region is 2:
        return region2.g_h(P, h)
    elif region is 4:
        return region4.g_h(P, h)
    else:
        return 0.000
def v_h(P, h, region = 0):
    """Specific volume [m^3 / kg]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.v_h(P, h)
    elif region is 2:
        return region2.v_h(P, h)
    elif region is 4:
        return region4.v_h(P, h)
    else:
        return 0.000
def u_h(P, h, region = 0):
    """Specific internal energy [kJ / kg]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.u_h(P, h)
    elif region is 2:
        return region2.u_h(P, h)
    elif region is 4:
        return region4.u_h(P, h)
    else:
        return 0.000
def s_h(P, h, region = 0):
    """Specific entropy [kJ / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.s_h(P, h)
    elif region is 2:
        return region2.s_h(P, h)
    elif region is 4:
        return region4.s_h(P, h)
    else:
        return 0.000
def T_h(P, h, region = 0):
    """ Temperature [K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.T_h(P, h)
    elif region is 2:
        return region2.T_h(P, h)
    elif region is 4:
        return region4.satT(P)
    else:
        return 0.000
def cp_h(P, h, region = 0):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.cp_h(P, h)
    elif region is 2:
        return region2.cp_h(P, h)
    elif region is 4:
        return region4.cp_h(P, h)
    else:
        return 0.000
def cv_h(P, h, region = 0):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.cv_h(P, h)
    elif region is 2:
        return region2.cv_h(P, h)
    elif region is 4:
        return region4.cv_h(P, h)
    else:
        return 0.000
def w_h(P, h, region = 0):
    """ Speed of sound [m / s]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.w_h(P, h)
    elif region is 2:
        return region2.w_h(P, h)
    elif region is 4:
        return region4.w_h(P, h)
    else:
        return 0.000
def a_h(P, h, region = 0):
    """Isobaric cubic expansion coefficient [1 / K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.a_h(P, h)
    elif region is 2:
        return region2.a_h(P, h)
    elif region is 4:
        return region4.a_h(P, h)
    else:
        return 0.000
def k_h(P, h, region = 0):
    """Isothermal compressibility [kg / kJ]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.k_h(P, h)
    elif region is 2:
        return region2.k_h(P, h)
    elif region is 4:
        return region4.k_h(P, h)
    else:
        return 0.000

#### water property derivatives ####
def dgdP_h(P, h, region = 0):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dgdP_h(P, h)
    elif region is 2:
        return region2.dgdP_h(P, h)
    elif region is 4:
        return region4.dgdP_h(P, h)
    else:
        return 0.000
def dvdP_h(P, h, region = 0):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dvdp_h(P, h)
    elif region is 2:
        return region2.dvdP_h(P, h)
    elif region is 4:
        return region4.dvdP_h(P, h)
    else:
        return 0.000
def dudP_h(P, h, region = 0):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dudP_h(P, h)
    elif region is 2:
        return region2.dudP_h(P, h)
    elif region is 4:
        return region4.dudP_h(P, h)
    else:
        return 0.000
def dsdP_h(P, h, region = 0):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dsdP_h(P, h)
    elif region is 2:
        return region2.dsdP_h(P, h)
    elif region is 4:
        return region4.dsdP_h(P, h)
    else:
        return 0.000
def dhdP_h(P, h, region = 0):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return 0.000
    elif region is 2:
        return 0.000
    elif region is 4:
        return region4.dhdP_h(P, h)
    else:
        return 0.000
def dTdP_h(P, h, region = 0):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dTdP_h(P, h)
    elif region is 2:
        return region2.dTdP_h(P, h)
    elif region is 4:
        return region4.dTsdP(P)
    else:
        return 0.000

def dgdh_h(P, h, region = 0):
    """ Derivative of specific gibbs free energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dgdh_h(P, h)
    elif region is 2:
        return region2.dgdh_h(P, h)
    elif region is 4:
        return region4.dgdh_h(P, h)
    else:
        return 0.000
def dvdh_h(P, h, region = 0):
    """ Derivative of specific volume [m^3 kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dvdh_h(P, h)
    elif region is 2:
        return region2.dvdh_h(P, h)
    elif region is 4:
        return region4.dvdh_h(P, h)
    else:
        return 0.000
def dudh_h(P, h, region = 0):
    """ Derivative of specific internal energy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dudh_h(P, h)
    elif region is 2:
        return region2.dudh_h(P, h)
    elif region is 4:
        return region4.dudh_h(P, h)
    else:
        return 0.000
def dsdh_h(P, h, region = 0):
    """ Derivative of specific entropy [kJ kg / kg K kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dsdh_h(P, h)
    elif region is 2:
        return region2.dsdh_h(P, h)
    elif region is 4:
        return region4.dsdh_h(P, h)
    else:
        return 0.000
def dhdh_h(P, h, region = 0):
    """ Derivative of specific enthalpy [kJ kg / kg kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return 1.000
    elif region is 2:
        return 1.000
    elif region is 4:
        return 1.000
    else:
        return 0.000
def dTdh_h(P, h, region = 0):
    """ Derivative of Temperature [K kg / kJ]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dTdh_h(P, h)
    elif region is 2:
        return region2.dTdh_h(P, h)
    elif region is 4:
        return 0.000
    else:
        return 0.000

###########################################################
#####           Pressure-Entropy Formulation          #####
###########################################################
def idRegion_s(P, s):
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # supporting boundaries
    Tbnd01 = region1.Tbnd01
    Pbnd4  = satP(Tbnd01)
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13
    Tbnd32 = region3.bnd23T(min(max(P, 16.5292), 100.0))
    Tbnd4  = satT(P)

    # Enthalpy- pressure boundaries
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

    if (P >= Pbnd0) and (s >= sbnd01) and (P <= Pbnd1) and (s <= sbnd25):
        if (P >= Pbndh1):
            if (s <= sbnd13):
                region = 1
            elif (s >= sbnd32):
                region = 2
            else:
                # region 3 via P,h relations not implemented
                region = 0
        else:
            if (s <= sbnd14):
                region = 1
            elif (s >= sbnd42):
                region = 2
            else:
                region = 4
    assert (region is not 0), "Water properties not avalable!"
    return region

#### water properties ####
def g_s(P, s, region = 0):
    """Specific gibbs free energy [kJ / kg]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.g_s(P, s)
    elif region is 2:
        return region2.g_s(P, s)
    elif region is 4:
        return region4.g_s(P, s)
    else:
        return 0.000
def v_s(P, s, region = 0):
    """Specific volume [m^3 / kg]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.v_s(P, s)
    elif region is 2:
        return region2.v_s(P, s)
    elif region is 4:
        return region4.v_s(P, s)
    else:
        return 0.000
def u_s(P, s, region = 0):
    """Specific internal energy [kJ / kg]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.u_s(P, s)
    elif region is 2:
        return region2.u_s(P, s)
    elif region is 4:
        return region4.u_s(P, s)
    else:
        return 0.000
def T_s(P, s, region = 0):
    """ Temperature [K]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.T_s(P, s)
    elif region is 2:
        return region2.T_s(P, s)
    elif region is 4:
        return region4.satT(P)
    else:
        return 0.000
def h_s(P, s, region = 0):
    """Specific entropy [kJ / kg]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.h_s(P, s)
    elif region is 2:
        return region2.h_s(P, s)
    elif region is 4:
        return region4.h_s(P, s)
    else:
        return 0.000
def cp_s(P, s, region = 0):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.cp_s(P, s)
    elif region is 2:
        return region2.cp_s(P, s)
    elif region is 4:
        return region4.cp_s(P, s)
    else:
        return 0.000
def cv_s(P, s, region = 0):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.cv_s(P, s)
    elif region is 2:
        return region2.cv_s(P, s)
    elif region is 4:
        return region4.cv_s(P, s)
    else:
        return 0.000
def w_s(P, s, region = 0):
    """ Speed of sound [m / s]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.w_s(P, s)
    elif region is 2:
        return region2.w_s(P, s)
    elif region is 4:
        return region4.w_s(P, s)
    else:
        return 0.000
def a_s(P, s, region = 0):
    """Isobaric cubic expansion coefficient [1 / K]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.a_s(P, s)
    elif region is 2:
        return region2.a_s(P, s)
    elif region is 4:
        return region4.a_s(P, s)
    else:
        return 0.000
def k_s(P, s, region = 0):
    """Isothermal compressibility [kg / kJ]"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.k_s(P, s)
    elif region is 2:
        return region2.k_s(P, s)
    elif region is 4:
        return region4.k_s(P, s)
    else:
        return 0.000

#### water property derivatives ####
def dgdP_s(P, s, region = 0):
    """ Derivative of specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dgdP_s(P, s)
    elif region is 2:
        return region2.dgdP_s(P, s)
    elif region is 4:
        return region4.dgdP_s(P, s)
    else:
        return 0.000
def dvdP_s(P, s, region = 0):
    """ Derivative of specific volume [m^3 m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dvdP_s(P, s)
    elif region is 2:
        return region2.dvdP_s(P, s)
    elif region is 4:
        return region4.dvdP_s(P, s)
    else:
        return 0.000
def dudP_s(P, s, region = 0):
    """ Derivative of specific internal energy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dudP_s(P, s)
    elif region is 2:
        return region2.dudP_s(P, s)
    elif region is 4:
        return region4.dudP_s(P, s)
    else:
        return 0.000
def dsdP_s(P, s, region = 0):
    """ Derivative of specific entropy [kJ m^3 / kg K kJ]
    w.r.t pressure at constant specific/equilibrium entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return 0.000
    elif region is 2:
        return 0.000
    elif region is 4:
        return region4.dsdP_s(P, s)
    else:
        return 0.000
def dhdP_s(P, s, region = 0):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dhdP_s(P, s)
    elif region is 2:
        return region2.dhdP_s(P, s)
    elif region is 4:
        return region4.dhdP_s(P, s)
    else:
        return 0.000
def dTdP_s(P, s, region = 0):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant specific entropy"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dTdP_s(P, s)
    elif region is 2:
        return region2.dTdP_s(P, s)
    elif region is 4:
        return region4.dTsdP(P)
    else:
        return 0.000

def dgds_s(P, s, region = 0):
    """ Derivative of specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dgds_s(P, s)
    elif region is 2:
        return region2.dgds_s(P, s)
    elif region is 4:
        return region4.dgds_s(P, s)
    else:
        return 0.000
def dvds_s(P, s, region = 0):
    """ Derivative of specific volume [m^3 kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dvds_s(P, s)
    elif region is 2:
        return region2.dvds_s(P, s)
    elif region is 4:
        return region4.dvds_s(P, s)
    else:
        return 0.000
def duds_s(P, s, region = 0):
    """ Derivative of specific internal energy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.duds_s(P, s)
    elif region is 2:
        return region2.duds_s(P, s)
    elif region is 4:
        return region4.duds_s(P, s)
    else:
        return 0.000
def dsds_s(P, s, region = 0):
    """ Derivative of specific entropy [kJ kg K / kg K kJ]
    w.r.t specific entropy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return 1.000
    elif region is 2:
        return 1.000
    elif region is 4:
        return 1.000
    else:
        return 0.000
def dhds_s(P, s, region = 0):
    """ Derivative of specific enthalpy [kJ kg K / kg kJ]
    w.r.t specific entropy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dhds_s(P, s)
    elif region is 2:
        return region2.dhds_s(P, s)
    elif region is 4:
        return region4.dhds_s(P, s)
    else:
        return 0.000
def dTds_s(P, s, region = 0):
    """ Derivative of Temperature [K kg K / kJ]
    w.r.t enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_s(P, s)

    if region is 1:
        return region1.dTds_s(P, s)
    elif region is 2:
        return region2.dTds_s(P, s)
    elif region is 4:
        return 0.000
    else:
        return 0.000

###########################################################
#####     Pressure Only (Saturation) Formulation      #####
###########################################################

#### P-T saturation curves ####
def satP(T):
    """ Saturation Pressure [Mpa]
    for specified Temperature"""

    return region4.satP(T)
def satT(P):
    """ Saturation Temperature [K]
    for specified Pressure"""

    return region4.satT(P)

#### Saturated liquid properties ####
def gf(P):
    """ Specific gibbs free energy [kJ / kg]
    of saturated liquid"""

    return region4.gf(P)
def vf(P):
    """ Specific volume [m^3 / kg]
    of saturated liquid"""

    return region4.vf(P)
def uf(P):
    """ Specific internal energy [kJ / kg]
    of saturated liquid"""

    return region4.uf(P)
def sf(P):
    """ Specific entropy [kJ / kg K]
    of saturated liquid"""

    return region4.sf(P)
def hf(P):
    """ Specific enthalpy [kJ / kg]
    of saturated liquid"""

    return region4.hf(P)
def cpf(P):
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated liquid"""

    return region4.cpf(P)
def cvf(P):
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated liquid"""

    return region4.cvf(P)
def wf(P):
    """ Speed of sound [m / s]
    of saturated liquid"""

    return region4.wf(P)
def af(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated liquid"""

    return region4.af(P)
def kf(P):
    """Isothermal compressibility [kg / kJ]
    of saturated liquid"""

    return region4.kf(P)

#### Saturated vapor properties ####
def gg(P):
    """ Specific gibbs free energy [kJ / kg]
    of saturated vapor"""

    return region4.gg(P)
def vg(P):
    """ Specific volume [m^3 / kg]
    of saturated vapor"""

    return region4.vg(P)
def ug(P):
    """ Specific internal energy [kJ / kg]
    of saturated vapor"""
    
    return region4.ug(P)
def sg(P):
    """ Specific entropy [kJ / kg K]
    of saturated vapor"""

    return region4.sg(P)
def hg(P):
    """ Specific enthalpy [kJ / kg]
    of saturated vapor"""
    
    return region4.hg(P)
def cpg(P):
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated vapor"""

    return region4.cpg(P)
def cvg(P):
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated vapor"""

    return region4.cvg(P)
def wg(P):
    """ Speed of sound [m / s]
    of saturated vapor"""

    return region4.wg(P)
def ag(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated vapor"""

    return region4.ag(P)
def kg(P):
    """Isothermal compressibility [kg / kJ]
    of saturated vapor"""

    return region4.kg(P)

#### delta saturation properties ####
def gfg(P):
    """ Specific gibbs free energy; [kJ / kg]
    saturation rise of"""

    return region4.gfg(P)
def vfg(P):
    """ Specific volume; [m^3 / kg]
    saturation rise of"""

    return region4.vfg(P)
def ufg(P):
    """ Specific internal energy; [kJ / kg]
    saturation rise of"""

    return region4.ufg(P)
def sfg(P):
    """ Specific entropy; [kJ / kg K]
    saturation rise of"""

    return region4.sfg(P)
def hfg(P):
    """ Specific enthalpy; [kJ / kg]
    saturation rise of"""

    return region4.hfg(P)
def cpfg(P):
    """ Specific isobaric heat capacity; [kJ / kg K]
    saturation rise of"""

    return region4.cpfg(P)
def cvfg(P):
    """ Specific isochoric heat capacity; [kJ / kg K]
    saturation rise of"""

    return region4.cvfg(P)
def wfg(P):
    """ Speed of sound; [m / s]
    saturation rise of"""

    return region4.wfg(P)
def afg(P):
    """Isobaric cubic expansion coefficient; [1 / K]
    saturation rise of"""

    return region4.afg(P)
def kfg(P):
    """Isothermal compressibility; [kg / kJ]
    saturation rise of"""

    return region4.kfg(P)

#### Saturated liquid derivatives ####
def dgfdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""

    return region4.dgfdP(P)
def dvfdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""

    return region4.dvfdP(P)
def dufdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""

    return region4.dufdP(P)
def dsfdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated liquid w.r.t. pressure"""

    return region4.dsfdP(P)
def dhfdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""

    return region4.dhfdP(P)

#### Saturated vapor derivatives ####
def dggdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""

    return region4.dggdP(P)
def dvgdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""

    return region4.dvgdP(P)
def dugdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""

    return region4.dugdP(P)
def dsgdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated vapor w.r.t. pressure"""

    return region4.dsgdP(P)
def dhgdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""

    return region4.dhgdP(P)

#### Delta saturation derivatives ####
def dgfgdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return region4.dgfgdP(P)
def dvfgdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return region4.dvfgdP(P)
def dufgdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return region4.dufgdP(P)
def dsfgdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure; saturation rise of"""

    return region4.dsfgdP(P)
def dhfgdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return region4.dhfgdP(P)
