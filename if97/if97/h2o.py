from if97 import region1, region2, region3, region4

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
    Tbnd13 = region1.Tbnd13
    Tbnd25 = region2.Tbnd25

    # non-constant boundaries
    Pbnd32 = region3.bnd23P(min(max(T, Tbnd13), 863.15))
    Pbnd4  = satP(T)

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
    """Specific volume [m^3 / kg K]"""
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.v(P, T)
    elif region is 2:
        return region2.v(P, T)
    else:
        return 0.000
def u(P, T, region = 0):
    """Specific internal energy [kJ / kg K]"""
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
    """Specific enthalpy [kJ / kg K]"""
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
    return -v(P, T) * a(P, T)
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
    """ Derivative of specific entropy [kJ / kg K]
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

    # Constant boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1 
    Tbnd01 = region1.Tbnd01
    Tbnd13 = region1.Tbnd13
    Tbnd25 = region2.Tbnd25

    # non-constant boundaries
    Pbnd4  = satP
    Tbnd32 = region3.bnd23T
    Tbnd4  = satT

    # Enthalpy- pressure boundaries
    Pbndh1 = satP(Tbnd13)
    hbnd01 = region1.h(Pbnd4(Tbnd01), Tbnd01)
    hbnd40 = region2.h(Pbnd1, Tbnd25)
    hbnd13 = region1.h(P, Tbnd13)
    hbnd14 = region1.h(P, Tbnd4(P))
    hbnd42 = region2.h(P, Tbnd4(P))
    hbnd32 = region2.h(P, Tbnd32(min(max(P, 16.5292), 100.0)))

    region = 0

    if (P >= Pbnd0) and (h >= hbnd01) and (P <= Pbnd1) and (h <= hbnd40):
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
    """Specific gibbs free energy [kJ / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.g_h(P, h)
    elif region is 2:
        return region2.g_h(P, h)
    elif region is 4:
        return g_e(P, h)
    else:
        return 0.000
def v_h(P, h, region = 0):
    """Specific volume [m^3 / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.v_h(P, h)
    elif region is 2:
        return region2.v_h(P, h)
    elif region is 4:
        return v_e(P, h)
    else:
        return 0.000
def u_h(P, h, region = 0):
    """Specific internal energy [kJ / kg K]"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.u_h(P, h)
    elif region is 2:
        return region2.u_h(P, h)
    elif region is 4:
        return u_e(P, h)
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
        return s_e(P, h)
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
        return satT(P)
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
        return cp_e(P, h)
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
        return cv_e(P, h)
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
        return w_e(P, h)
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
        return a_e(P, h)
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
        return k_e(P, h)
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
        return dgdP_e(P, h)
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
        return dvdP_e(P, h)
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
        return dudP_e(P, h)
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
        return dsdP_e(P, h)
    else:
        return 0.000
    return -v(P, h) * a(P, h)
def dhdP_h(P, h, region = 0):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
    w.r.t pressure at constant specific/equilibrium enthalpy"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return 0.000
    elif region is 2:
        return 0.000
    elif region is 4:
        return dhdP_e(P, h)
    else:
        return 0.000
def dTdP_h(P, h, region = 0):
    """ Derivative of Temperature [K m^3 / kJ]
    w.r.t pressure at constant enthalpy"""
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
    """ Derivative of specific gibbs free energy [kJ / kg K]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dgdh_h(P, h)
    elif region is 2:
        return region2.dgdh_h(P, h)
    elif region is 4:
        return dgdh_e(P, h)
    else:
        return 0.000
def dvdh_h(P, h, region = 0):
    """ Derivative of specific volume [m^3 / kg K]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dvdh_h(P, h)
    elif region is 2:
        return region2.dvdh_h(P, h)
    elif region is 4:
        return dvdh_e(P, h)
    else:
        return 0.000
def dudh_h(P, h, region = 0):
    """ Derivative of specific internal energy [kJ / kg K]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dudh_h(P, h)
    elif region is 2:
        return region2.dudh_h(P, h)
    elif region is 4:
        return dudh_e(P, h)
    else:
        return 0.000
def dsdh_h(P, h, region = 0):
    """ Derivative of specific entropy [kJ / kg K]
    w.r.t specific enthalpy at constant pressure"""
    if region is 0:
        region = idRegion_h(P, h)

    if region is 1:
        return region1.dsdh_h(P, h)
    elif region is 2:
        return region2.dsdh_h(P, h)
    elif region is 4:
        return dsdh_e(P, h)
    else:
        return 0.000
def dhdh_h(P, h, region = 0):
    """ Derivative of specific enthalpy [kJ m^3 / kg kJ]
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
    w.r.t enthalpy at constant pressure"""
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

    return region1.g(P, satT(P))
def vf(P):
    """ Specific volume [m^3 / kg]
    of saturated liquid"""

    return region1.v(P, satT(P))
def uf(P):
    """ Specific internal energy [kJ / kg]
    of saturated liquid"""

    return region1.u(P, satT(P))
def sf(P):
    """ Specific entropy [kJ / kg K]
    of saturated liquid"""

    return region1.s(P, satT(P))
def hf(P):
    """ Specific enthalpy [kJ / kg]
    of saturated liquid"""

    return region1.h(P, satT(P))
def cpf(P):
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated liquid"""

    return region1.cp(P, satT(P))
def cvf(P):
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated liquid"""

    return region1.cv(P, satT(P))
def wf(P):
    """ Speed of sound [m / s]
    of saturated liquid"""

    return region1.w(P, satT(P))
def af(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated liquid"""

    return region1.a(P, satT(P))
def kf(P):
    """Isothermal compressibility [kg / kJ]
    of saturated liquid"""

    return region1.k(P, satT(P))

#### Saturated vapor properties ####
def gg(P):
    """ Specific gibbs free energy [kJ / kg]
    of saturated vapor"""

    return region2.g(P, satT(P))
def vg(P):
    """ Specific volume [m^3 / kg]
    of saturated vapor"""

    return region2.v(P, satT(P))
def ug(P):
    """ Specific internal energy [kJ / kg]
    of saturated vapor"""
    
    return region2.u(P, satT(P))
def sg(P):
    """ Specific entropy [kJ / kg K]
    of saturated vapor"""

    return region2.s(P, satT(P))
def hg(P):
    """ Specific enthalpy [kJ / kg]
    of saturated vapor"""
    
    return region2.h(P, satT(P))
def cpg(P):
    """ Specific isobaric heat capacity [kJ / kg K]
    of saturated vapor"""

    return region2.cp(P, satT(P))
def cvg(P):
    """ Specific isochoric heat capacity [kJ / kg K]
    of saturated vapor"""

    return region2.cv(P, satT(P))
def wg(P):
    """ Speed of sound [m / s]
    of saturated vapor"""

    return region2.w(P, satT(P))
def ag(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated vapor"""

    return region2.a(P, satT(P))
def kg(P):
    """Isothermal compressibility [kg / kJ]
    of saturated vapor"""

    return region2.k(P, satT(P))

#### delta saturation properties ####
def gfg(P):
    """ Specific gibbs free energy; [kJ / kg]
    saturation rise of"""

    return gg(P) - gf(P)
def vfg(P):
    """ Specific volume; [m^3 / kg]
    saturation rise of"""

    return vg(P) - vf(P)
def ufg(P):
    """ Specific internal energy; [kJ / kg]
    saturation rise of"""

    return ug(P) - uf(P)
def sfg(P):
    """ Specific entropy; [kJ / kg K]
    saturation rise of"""

    return sg(P) - sf(P)
def hfg(P):
    """ Specific enthalpy; [kJ / kg]
    saturation rise of"""

    return hg(P) - hf(P)
def cpfg(P):
    """ Specific isobaric heat capacity; [kJ / kg K]
    saturation rise of"""

    return cpg(P) - cpf(P)
def cvfg(P):
    """ Specific isochoric heat capacity; [kJ / kg K]
    saturation rise of"""

    return cvg(P) - cvf(P)
def wfg(P):
    """ Speed of sound; [m / s]
    saturation rise of"""

    return wg(P) - wf(P)
def afg(P):
    """Isobaric cubic expansion coefficient; [1 / K]
    saturation rise of"""

    return ag(P) - af(P)
def kfg(P):
    """Isothermal compressibility; [kg / kJ]
    saturation rise of"""

    return kg(P) - kf(P)

#### Saturated liquid derivatives ####
def dgfdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)

    return region1.dgdP(P, T) + region1.dgdT(P, T) * region4.dTsdP(P)
def dvfdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)

    return region1.dvdP(P, T) + region1.dvdT(P, T) * region4.dTsdP(P)
def dufdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)

    return region1.dudP(P, T) + region1.dudT(P, T) * region4.dTsdP(P)
def dsfdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)

    return region1.dsdP(P, T) + region1.dsdT(P, T) * region4.dTsdP(P)
def dhfdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)

    return region1.dhdP(P, T) + region1.dhdT(P, T) * region4.dTsdP(P)

#### Saturated vapor derivatives ####
def dggdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)

    return region2.dgdP(P, T) + region2.dgdT(P, T) * region4.dTsdP(P)
def dvgdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)

    return region2.dvdP(P, T) + region2.dvdT(P, T) * region4.dTsdP(P)
def dugdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)

    return region2.dudP(P, T) + region2.dudT(P, T) * region4.dTsdP(P)
def dsgdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)

    return region2.dsdP(P, T) + region2.dsdT(P, T) * region4.dTsdP(P)
def dhgdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)

    return region2.dhdP(P, T) + region2.dhdT(P, T) * region4.dTsdP(P)

#### Delta saturation derivatives ####
def dgfgdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return dggdP(P) - dgfdP(P)
def dvfgdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return dvgdP(P) - dvfdP(P)
def dufgdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return dugdP(P) - dufdP(P)
def dsfgdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure; saturation rise of"""

    return dsgdP(P) - dsfdP(P)
def dhfgdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure; saturation rise of"""

    return dhgdP(P) - dhfdP(P)

###########################################################
#####    Two-phase Mixture (Saturation)   Properties  #####
###########################################################

#### equilibrium quantities ####
def g_e(P, h):
    """ Equilibrium specific gibbs free energy [kJ / kg]"""
    x = x_e(P, h)

    return x * gg(P) + (1 - x) * gf(P)
def v_e(P, h):
    """ Equilibrium specific volume [m^3 / kg]"""
    x = x_e(P, h)

    return x * vg(P) + (1 - x) * vf(P)
def u_e(P, h):
    """ Equilibrium specific internal energy [kJ / kg]"""
    x = x_e(P, h)

    return x * ug(P) + (1 - x) * uf(P)
def s_e(P, h):
    """ Equilibrium specific entropy [kJ / kg K]"""
    x = x_e(P, h)

    return x * sg(P) + (1 - x) * sf(P)
def x_e(P, h):
    """ Equilibrium quality [-]
    as a function of P, h"""

    return (h - hf(P)) / hfg(P)
def cp_e(P, h):
    """ Equilibrium specific isobaric heat capacity [kJ / kg K]"""
    x = x_e(P, h)

    return x * cpg(P) + (1 - x) * cpf(P)
def cv_e(P, h):
    """ Equilibrium specific isochoric heat capacity [kJ / kg K]"""
    x = x_e(P, h)

    return x * cvg(P) + (1 - x) * cvf(P)
def w_e(P, h):
    """ Equilibrium speed of sound [m / s]"""
    x = x_e(P, h)

    return x * wg(P) + (1 - x) * wf(P)
def a_e(P, h):
    """ Equilibrium Isobaric cubic expansion coefficient [1 / K]"""
    x = x_e(P, h)

    return x * ag(P) + (1 - x) * af(P)
def k_e(P, h):
    """ Equilibrium Isothermal compressibility [kg / kJ]"""
    x = x_e(P, h)

    return x * kg(P) + (1 - x) * kf(P)
def h_e(P, x):
    """ Equilibrium specific enthaply [kJ / kg]"""

    return x * hg(P) + (1 - x) * hf(P)

#### equilibrium derivatives ####
def dgdP_e(P, h):
    """ Derivative of equilibrium specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_e(P, h)

    return x * dggdP(P) + (1 - x) * dgfdP(P) + gfg(P) * dxdP_e(P, h)
def dvdP_e(P, h):
    """ Derivative of equilibrium specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_e(P, h)

    return x * dvgdP(P) + (1 - x) * dvfdP(P) + vfg(P) * dxdP_e(P, h)
def dudP_e(P, h):
    """ Equilibrium specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_e(P, h)

    return x * dugdP(P) + (1 - x) * dufdP(P) + ufg(P) * dxdP_e(P, h)
def dsdP_e(P, h):
    """ Derivative of equilibrium specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_e(P, h)

    return x * dsgdP(P) + (1 - x) * dsfdP(P) + sfg(P) * dxdP_e(P, h)
def dxdP_e(P, h):
    """ Derivative of equilibrium quality [m^3 / kJ]
    w.r.t. pressure"""

    return (dhfgdP(P) * (hf(P) - h) - dhfdP(P) * hfg(P)) / hfg(P)**2
def dhdP_e(P, h):
    """ Derivative of equilibrium specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_e(P, h)

    return x * dhgdP(P) + (1 - x) * dhfdP(P) + dxdP_e(P, h) * dhfgdP(P)

def dgdh_e(P, h):
    """ Derivative of equilibrium specific gibbs free energy [kJ kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""

    return gfg(P) * dxdh_e(P)
def dvdh_e(P, h):
    """ Derivative of equilibrium specific volume [m^3 kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""

    return vfg(P) * dxdh_e(P)
def dudh_e(P, h):
    """ Equilibrium specific internal energy [kJ kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""

    return ufg(P) * dxdh_e(P)
def dsdh_e(P, h):
    """ Derivative of equilibrium specific entropy [kJ kg / kg K kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""

    return sfg(P) * dxdh_e(P)
def dxdh_e(P):
    """ Derivative of equilibrium quality [kg / kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""

    return 1 / hfg(P)