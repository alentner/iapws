# Constants for region 4
n = [ 0.11670521452767e+04, -0.72421316703206e+06, -0.17073846940092e+02, 0.12020824702470e+05, -0.32325550322333e+07, 0.14915108613530e+02, 
     -0.48232657361591e+04,  0.40511340542057e+06, -0.23855557567849e+00, 0.65017534844798e+03]

# Non-dimenionalization for region 4
Pb = 1.0    #[Mpa]
Tb = 1.0    #[K]

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
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    return Tb * (n[9] + D - ((n[9] + D)**2 - 4 * (n[8] + n[9] * D))**0.5) / 2