from . import region1, region2

###########################################################
#####       Constants and Dimensionless Functions     #####
###########################################################

# constants for Region 4
n = [ 0.11670521452767e+04, -0.72421316703206e+06, -0.17073846940092e+02, 0.12020824702470e+05, -0.32325550322333e+07, 0.14915108613530e+02, 
     -0.48232657361591e+04,  0.40511340542057e+06, -0.23855557567849e+00, 0.65017534844798e+03]

# non-dimenionalization for Region 4
Pb = 1.0 # [Mpa]
Tb = 1.0 # [K  ]

###########################################################
#####                Saturation Curves                #####
###########################################################

#### saturation curves ####
def satP(T):
    """ Saturation pressure as a function of temperature [MPa]"""
    v = T / Tb + n[8] / ((T / Tb) - n[9])

    A = v**2 + n[0] * v + n[1]
    B = n[2] * v**2 + n[3] * v + n[4]
    C = n[5] * v**2 + n[6] * v + n[7]

    return Pb * (2 * C / (-B + (B**2 - 4 * A * C)**0.5))**4

def satT(P):
    """ Saturation Temperature as a function of pressure [K]"""
    b = (P / Pb)**0.25

    E = b**2 + n[2] * b + n[5]
    F = n[0] * b**2 + n[3] * b + n[6]
    G = n[1] * b**2 + n[4] * b + n[7]
    D = 2 * G / (-F - (F**2 - 4 * E * G)**0.5)

    return Tb * (n[9] + D - ((n[9] + D)**2 - 4 * (n[8] + n[9] * D))**0.5) / 2

#### saturation curves derivatives ####
def dTsdP(P):
    """ Derivative of saturation temperature [K m^3 / kJ]
    w.r.t pressure """
    b = (P / Pb)**0.25

    E = b**2 + n[2] * b + n[5]
    F = n[0] * b**2 + n[3] * b + n[6]
    G = n[1] * b**2 + n[4] * b + n[7]
    D = 2 * G / (-F - (F**2 - 4 * E * G)**0.5)

    bp = 1 / (4 * Pb * (P / Pb)**0.75)

    Ep = (2 * b + n[2]) * bp
    Fp = (2 * n[0] * b + n[3]) * bp
    Gp = (2 * n[1] * b + n[4]) * bp
    Dp = (-F**2 * Ep + F * (Ep * (F**2 - 4 * E * G)**0.5 + E * Fp) - E * (-2 * G * Ep + Fp * (F**2 - 4 * E * G)**0.5 + 2 * E * Gp)) / (2 * E**2 * (F**2 - 4 * E * G)**0.5)
    
    return (Tb / 2) * (1 + (n[9] - D) / (n[9]**2 - 4 * n[8] - 2 * n[9] * D + D**2)**(1/2)) * Dp / (10**6 / 1000)

def dPsdT(T):
    """ Derivative of saturation pressure [P / K]
    w.r.t temperature """
    v = T / Tb + n[8] / ((T / Tb) - n[9])

    A = v**2 + n[0] * v + n[1]
    B = n[2] * v**2 + n[3] * v + n[4]
    C = n[5] * v**2 + n[6] * v + n[7]

    vp = 1 / Tb - n[8] * Tb / (T - n[9] * Tb)**2

    Ap = (2 * v + n[0]) * vp
    Bp = (2 * n[2] * v + n[3]) * vp
    Cp = (2 * n[5] * v + n[6]) * vp
    
    return 64 * Pb * C**3 * (C * Bp + Cp * (-B + (B**2 - 4 * A * C)**0.5) + C * (2 * C * Ap - B * Bp + 2 * A * Cp) / (B**2 - 4 * A * C)**0.5) / (-B + (B**2 - 4 * A * C)**0.5)**5

###########################################################
#####     Pressure Only (Saturation) Formulation      #####
###########################################################

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

def avf(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated liquid"""
    return region1.av(P, satT(P))

def kTf(P):
    """Isothermal compressibility [m^3 / kJ]
    of saturated liquid"""
    return region1.kT(P, satT(P))

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

def avg(P):
    """Isobaric cubic expansion coefficient [1 / K]
    of saturated vapor"""
    return region2.av(P, satT(P))

def kTg(P):
    """Isothermal compressibility [m^3 / kJ]
    of saturated vapor"""
    return region2.kT(P, satT(P))

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

def avfg(P):
    """Isobaric cubic expansion coefficient; [1 / K]
    saturation rise of"""
    return avg(P) - avf(P)

def kTfg(P):
    """Isothermal compressibility; [m^3 / kJ]
    saturation rise of"""
    return kTg(P) - kTf(P)

#### Saturated liquid derivatives ####
def dgfdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)
    return region1.dgdP(P, T) + region1.dgdT(P, T) * dTsdP(P)

def dvfdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)
    return region1.dvdP(P, T) + region1.dvdT(P, T) * dTsdP(P)

def dufdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)
    return region1.dudP(P, T) + region1.dudT(P, T) * dTsdP(P)

def dsfdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)
    return region1.dsdP(P, T) + region1.dsdT(P, T) * dTsdP(P)

def dhfdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated liquid w.r.t. pressure"""
    T = satT(P)
    return region1.dhdP(P, T) + region1.dhdT(P, T) * dTsdP(P)

#### Saturated vapor derivatives ####
def dggdP(P):
    """ Derivative of Specific gibbs free energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)
    return region2.dgdP(P, T) + region2.dgdT(P, T) * dTsdP(P)

def dvgdP(P):
    """ Derivative of Specific volume [m^3 m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)
    return region2.dvdP(P, T) + region2.dvdT(P, T) * dTsdP(P)

def dugdP(P):
    """ Derivative of Specific internal energy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)
    return region2.dudP(P, T) + region2.dudT(P, T) * dTsdP(P)

def dsgdP(P):
    """ Derivative of Specific entropy [kJ m^3 / kg K kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)
    return region2.dsdP(P, T) + region2.dsdT(P, T) * dTsdP(P)

def dhgdP(P):
    """ Derivative of Specific enthalpy [kJ m^3 / kg kJ]
    of saturated vapor w.r.t. pressure"""
    T = satT(P)
    return region2.dhdP(P, T) + region2.dhdT(P, T) * dTsdP(P)

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
###    Two-phase Mixture (Saturation) Formulation(P,h)  ###
###########################################################

#### equilibrium quantities ####
def g_h(P, h):
    """ Equilibrium specific gibbs free energy [kJ / kg]"""
    x = x_h(P, h)
    return x * gg(P) + (1 - x) * gf(P)

def v_h(P, h):
    """ Equilibrium specific volume [m^3 / kg]"""
    x = x_h(P, h)
    return x * vg(P) + (1 - x) * vf(P)

def u_h(P, h):
    """ Equilibrium specific internal energy [kJ / kg]"""
    x = x_h(P, h)
    return x * ug(P) + (1 - x) * uf(P)

def s_h(P, h):
    """ Equilibrium specific entropy [kJ / kg K]"""
    x = x_h(P, h)
    return x * sg(P) + (1 - x) * sf(P)

def h_h(P, x):
    """ Equilibrium specific enthaply [kJ / kg]"""
    return x * hg(P) + (1 - x) * hf(P)

def cp_h(P, h):
    """ Equilibrium specific isobaric heat capacity [kJ / kg K]"""
    x = x_h(P, h)
    return x * cpg(P) + (1 - x) * cpf(P)

def cv_h(P, h):
    """ Equilibrium specific isochoric heat capacity [kJ / kg K]"""
    x = x_h(P, h)
    return x * cvg(P) + (1 - x) * cvf(P)

def w_h(P, h):
    """ Equilibrium speed of sound [m / s]"""
    x = x_h(P, h)
    return x * wg(P) + (1 - x) * wf(P)

def av_h(P, h):
    """ Equilibrium Isobaric cubic expansion coefficient [1 / K]"""
    x = x_h(P, h)
    return x * avg(P) + (1 - x) * avf(P)

def kT_h(P, h):
    """ Equilibrium Isothermal compressibility [m^3 / kJ]"""
    x = x_h(P, h)
    return x * kTg(P) + (1 - x) * kTf(P)

def x_h(P, h):
    """ Equilibrium quality [-]
    as a function of P, h"""
    return (h - hf(P)) / hfg(P)

#### equilibrium derivatives ####
def dgdP_h(P, h):
    """ Derivative of equilibrium specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_h(P, h)
    return x * dggdP(P) + (1 - x) * dgfdP(P) + gfg(P) * dxdP_h(P, h)

def dvdP_h(P, h):
    """ Derivative of equilibrium specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_h(P, h)
    return x * dvgdP(P) + (1 - x) * dvfdP(P) + vfg(P) * dxdP_h(P, h)

def dudP_h(P, h):
    """ Equilibrium specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_h(P, h)
    return x * dugdP(P) + (1 - x) * dufdP(P) + ufg(P) * dxdP_h(P, h)

def dsdP_h(P, h):
    """ Derivative of equilibrium specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_h(P, h)
    return x * dsgdP(P) + (1 - x) * dsfdP(P) + sfg(P) * dxdP_h(P, h)

def dhdP_h(P, h):
    """ Derivative of equilibrium specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium enthalpy"""
    x = x_h(P, h)
    return x * dhgdP(P) + (1 - x) * dhfdP(P) + dxdP_h(P, h) * dhfgdP(P)

def dxdP_h(P, h):
    """ Derivative of equilibrium quality [m^3 / kJ]
    w.r.t. pressure"""
    return (dhfgdP(P) * (hf(P) - h) - dhfdP(P) * hfg(P)) / hfg(P)**2

def dgdh_h(P, h):
    """ Derivative of equilibrium specific gibbs free energy [kJ kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""
    return gfg(P) * dxdh_h(P)

def dvdh_h(P, h):
    """ Derivative of equilibrium specific volume [m^3 kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""
    return vfg(P) * dxdh_h(P)

def dudh_h(P, h):
    """ Equilibrium specific internal energy [kJ kg / kg kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""
    return ufg(P) * dxdh_h(P)

def dsdh_h(P, h):
    """ Derivative of equilibrium specific entropy [kJ kg / kg K kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""
    return sfg(P) * dxdh_h(P)

def dxdh_h(P):
    """ Derivative of equilibrium quality [kg / kJ]
    w.r.t. equilibrium enthalpy @ a given pressure"""
    return 1 / hfg(P)

###########################################################
###    Two-phase Mixture (Saturation) Formulation(P,s)  ###
###########################################################

#### equilibrium quantities ####
def g_s(P, s):
    """ Equilibrium specific gibbs free energy [kJ / kg]"""
    x = x_s(P, s)
    return x * gg(P) + (1 - x) * gf(P)

def v_s(P, s):
    """ Equilibrium specific volume [m^3 / kg]"""
    x = x_s(P, s)
    return x * vg(P) + (1 - x) * vf(P)

def u_s(P, s):
    """ Equilibrium specific internal energy [kJ / kg]"""
    x = x_s(P, s)
    return x * ug(P) + (1 - x) * uf(P)

def s_s(P, x):
    """ Equilibrium specific entropy [kJ / kg K]"""
    return x * sg(P) + (1 - x) * sf(P)

def h_s(P, s):
    """ Equilibrium specific enthaply [kJ / kg]"""
    x = x_s(P, s)
    return x * hg(P) + (1 - x) * hf(P)

def cp_s(P, s):
    """ Equilibrium specific isobaric heat capacity [kJ / kg K]"""
    x = x_s(P, s)
    return x * cpg(P) + (1 - x) * cpf(P)

def cv_s(P, s):
    """ Equilibrium specific isochoric heat capacity [kJ / kg K]"""
    x = x_s(P, s)
    return x * cvg(P) + (1 - x) * cvf(P)

def w_s(P, s):
    """ Equilibrium speed of sound [m / s]"""
    x = x_s(P, s)
    return x * wg(P) + (1 - x) * wf(P)

def av_s(P, s):
    """ Equilibrium Isobaric cubic expansion coefficient [1 / K]"""
    x = x_s(P, s)
    return x * avg(P) + (1 - x) * avf(P)

def kT_s(P, s):
    """ Equilibrium Isothermal compressibility [m^3 / kJ]"""
    x = x_s(P, s)
    return x * kTg(P) + (1 - x) * kTf(P)

def x_s(P, s):
    """ Equilibrium quality [-]"""
    return (s - sf(P)) / sfg(P)

#### equilibrium derivatives ####
def dgdP_s(P, s):
    """ Derivative of equilibrium specific gibbs free energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium entropy"""
    x = x_s(P, s)
    return x * dggdP(P) + (1 - x) * dgfdP(P) + gfg(P) * dxdP_s(P, s)

def dvdP_s(P, s):
    """ Derivative of equilibrium specific volume [m^3 m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium entropy"""
    x = x_s(P, s)
    return x * dvgdP(P) + (1 - x) * dvfdP(P) + vfg(P) * dxdP_s(P, s)

def dudP_s(P, s):
    """ Equilibrium specific internal energy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium entropy"""
    x = x_s(P, s)
    return x * dugdP(P) + (1 - x) * dufdP(P) + ufg(P) * dxdP_s(P, s)

def dsdP_s(P, s):
    """ Derivative of equilibrium specific entropy [kJ m^3 / kg K kJ]
    w.r.t. pressure @ a given equilibrium entropy"""
    x = x_s(P, s)
    return x * dsgdP(P) + (1 - x) * dsfdP(P) + dxdP_s(P, s) * dsfgdP(P)

def dhdP_s(P, s):
    """ Derivative of equilibrium specific enthalpy [kJ m^3 / kg kJ]
    w.r.t. pressure @ a given equilibrium entropy"""
    x = x_s(P, s)
    return x * dhgdP(P) + (1 - x) * dhfdP(P) + hfg(P) * dxdP_s(P, s)

def dxdP_s(P, s):
    """ Derivative of equilibrium quality [m^3 / kJ]
    w.r.t. pressure"""
    return (dsfgdP(P) * (sf(P) - s) - dsfdP(P) * sfg(P)) / sfg(P)**2

def dgds_s(P, s):
    """ Derivative of equilibrium specific gibbs free energy [kJ kg K / kg kJ]
    w.r.t. equilibrium entropy @ a given pressure"""
    return gfg(P) * dxds_s(P)

def dvds_s(P, s):
    """ Derivative of equilibrium specific volume [m^3 kg K / kg kJ]
    w.r.t. equilibrium entropy @ a given pressure"""
    return vfg(P) * dxds_s(P)

def duds_s(P, s):
    """ Equilibrium specific internal energy [kJ kg K / kg kJ]
    w.r.t. equilibrium entropy @ a given pressure"""
    return ufg(P) * dxds_s(P)

def dhds_s(P, s):
    """ Derivative of equilibrium specific enthalpy [kJ kg K / kg kJ]
    w.r.t. equilibrium entropy @ a given pressure"""
    return hfg(P) * dxds_s(P)

def dxds_s(P):
    """ Derivative of equilibrium quality [kg K / kJ]
    w.r.t. equilibrium entropy @ a given pressure"""
    return 1 / sfg(P)
