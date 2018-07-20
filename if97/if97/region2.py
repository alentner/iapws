from math import log

# Constants for region 2
J0 = [0, 1, -5, -4, -3, -2, -1, 2, 3]
n0 = [-0.96927686500217e1,  0.10086655968018e2, -0.56087911283020e-2, 0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1,
      -0.43839511319450e1, -0.28408632460772,    0.21268463753307e-1]
Ir = [1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  5,  6,  6, 6, 
      7,  7,  7,  8,  8,  9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
Jr = [0,  1,  2,  3,  6,  1,  2,  4,  7, 36,  0,  1,  3,  6, 35,  1,  2,  3,  7,  3, 16, 35,
      0, 11, 25,  8, 36, 13,  4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
nr = [-0.0017731742473212999,  -0.017834862292357999,   -0.045996013696365003,   -0.057581259083432,     -0.050325278727930002,   -3.3032641670203e-05,
      -0.00018948987516315,    -0.0039392777243355001,  -0.043797295650572998,   -2.6674547914087001e-05, 2.0481737692308999e-08,  4.3870667284435001e-07, 
      -3.2277677238570002e-05, -0.0015033924542148,     -0.040668253562648998,   -7.8847309559367001e-10, 1.2790717852285001e-08,  4.8225372718507002e-07,
       2.2922076337661001e-06, -1.6714766451061001e-11, -0.0021171472321354998, -23.895741934103999,     -5.9059564324270004e-18, -1.2621808899101e-06, 
      -0.038946842435739003,    1.1256211360459e-11,    -8.2311340897998004,      1.9809712802088e-08,    1.0406965210174e-19,    -1.0234747095929e-13, 
      -1.0018179379511e-09,    -8.0882908646984998e-11,  0.10693031879409,       -0.33662250574170999,    8.9185845355420999e-25,  3.0629316876231997e-13, 
      -4.2002467698208001e-06, -5.9056029685639003e-26,  3.7826947613457002e-06, -1.2768608934681e-15,    7.3087610595061e-29,     5.5414715350778001e-17,
      -9.4369707241209998e-07]

# Non-dimenionalization for region 2
Ps = 1.0        #[Mpa]
Ts = 540.0      #[K]
R  = 0.461526   #[kJ / kg K]

def gamma(pi, tau):
    """ Dimensionless form for the specific Gibbs free energy"""
    sum = log(pi)
    for Ji, ni in zip(J0, n0):
        sum += ni * tau**Ji

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * (tau - 0.5)**Ji
    return sum

def gamma_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""
    sum = 1 / pi

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gamma_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi)"""   
    sum = -1 / pi**2

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gamma_tau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * tau**(Ji - 1)

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gamma_tautau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau)"""
    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * (Ji - 1) * tau**(Ji - 2)

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * (Ji - 1) * (tau - 0.5)**(Ji - 2)
    return sum

def gamma_pitau(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi) and temperature (tau)"""
    sum = 0

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def gammaR_pi(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""
    sum = 0

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * (tau - 0.5)**Ji
    return sum

def gammaR_pipi(pi, tau):
    """ Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""   
    sum = 0

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * (Ii - 1) * pi**(Ii - 2) * (tau - 0.5)**Ji
    return sum

def gammaR_pitau(pi, tau):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); residual part"""
    sum = 0

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * Ji * (tau - 0.5)**(Ji - 1)
    return sum

def g(P, T):
    """ Specific Gibbs free energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T

    return gamma(pi, tau) * R * T

def v(P, T):
    """ Specific volume [m^3 / kg]"""
    pi = P / Ps
    tau = Ts / T

    return pi * gamma_pi(pi, tau) * R * T / (P * 10**6 / 1000)

def u(P, T):
    """ Specific internal energy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T

    return (tau * gamma_tau(pi, tau) - pi * gamma_pi(pi, tau)) * R * T    

def s(P, T):
    """ Specific entropy [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T

    return (tau * gamma_tau(pi, tau) - gamma(pi, tau)) * R

def h(P, T):
    """ Specific enthalpy [kJ / kg]"""
    pi = P / Ps
    tau = Ts / T

    return tau * gamma_tau(pi, tau) * R * T

def cp(P, T):
    """ Specific isobaric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T

    return -tau**2 * gamma_tautau(pi, tau) * R

def cv(P, T):
    """ Specific isochoric heat capacity [kJ / kg K]"""
    pi = P / Ps
    tau = Ts / T

    return (-tau**2 * gamma_tautau(pi, tau) - (1 + pi * gammaR_pi(pi, tau) - tau * gammaR_pitau(pi, tau))**2 / (1 - pi**2 * gammaR_pipi(pi, tau))) * R

def w(P, T):
    """ Speed of sound [m / s]"""
    pi = P / Ps
    tau = Ts / T

    return (R * T * 1000 * (1 + 2 * pi * gammaR_pi(pi, tau) + pi**2 * gammaR_pi(pi, tau)**2) / ((1 - pi**2 * gammaR_pipi(pi, tau)) + (1 + pi * gammaR_pi(pi, tau) - tau * pi * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau))))**0.5