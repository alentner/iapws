# Constants for region 4
n = [ 0.11670521452767e+04, -0.72421316703206e+06, -0.17073846940092e+02, 0.12020824702470e+05, -0.32325550322333e+07, 0.14915108613530e+02, 
     -0.48232657361591e+04,  0.40511340542057e+06, -0.23855557567849e+00, 0.65017534844798e+03]

# Non-dimenionalization for region 4
Pb = 1.0    #[Mpa]
Tb = 1.0    #[K]

#### Saturation curves ####
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

#### Saturation curves derivatives ####
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
    Dp = (-F**2 * Ep + F * (Ep * (F**2 - 4 * E * G)**0.5 + E * Fp) - \
                       E * (-2 * G * Ep + Fp * (F**2 - 4 * E * G)**0.5 + 2 * E * Gp)) / \
         (2 * E**2 * (F**2 - 4 * E * G)**0.5)
    
    return (Tb / 2) * (1 + (n[9] - D) / (n[9]**2 - 4 * n[8] - 2 * n[9] * D + D**2)**(1/2)) * Dp / (10**6 / 1000)
def dPsdT(T):
    """ Derivative of saturation pressure [P / K]
    w.r.t temperature """
    v = T / Tb + n[8] / ((T / Tb) - n[9])

    A = v**2 + n[0] * v + n[1]
    B = n[2] * v**2 + n[3] * v + n[4]
    C = n[5] * v**2 + n[6] * v + n[7]

    #return Pb * (2 * C / (-B + (B**2 - 4 * A * C)**0.5))**4

    vp = 1 / Tb - n[8] * Tb / (T - n[9] * Tb)**2

    Ap = (2 * v + n[0]) * vp
    Bp = (2 * n[2] * v + n[3]) * vp
    Cp = (2 * n[5] * v + n[6]) * vp
    
    return 64 * Pb * C**3 * (C * Bp + Cp * (-B + (B**2 - 4 * A * C)**0.5) + C * (2 * C * Ap - B * Bp + 2 * A * Cp) / (B**2 - 4 * A * C)**0.5) / \
                (-B + (B**2 - 4 * A * C)**0.5)**5