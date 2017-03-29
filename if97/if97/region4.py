# Constants for region 4
n = [0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02, 0.12020824702470E+05, -0.32325550322333E+07, 0.14915108613530E+02, -0.48232657361591E+04, 0.40511340542057E+06, -0.23855557567849E+00, 0.65017534844798E+03]

# Non-dimenionalization for region 1
Ps = 1.0    #[Mpa]
Ts = 1.0    #[K]

def satP(T):
    v = T / Ts + n[8] / ((T / Ts) - n[9])

    A = v**2 + n[0] * v + n[1]
    B = n[2] * v**2 + n[3] * v + n[4]
    C = n[5] * v**2 + n[6] * v + n[7]

    return Ps * (2 * C / (-B + (B**2 - 4 * A * C)**0.5))**4
def satT(P):
    b = (P / Ps)**0.25

    E = b**2 + n[2] * b + n[5]
    F = n[0] * b**2 + n[3] * b + n[6]
    G = n[1] * b**2 + n[4] * b + n[7]
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    return Ts * (n[9] + D - ((n[9] + D)**2 - 4 * (n[8] + n[9] * D))**0.5) / 2