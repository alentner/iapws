from if97 import region1, region2, region4

# Constants for region 4
n = region4.n

# Non-dimenionalization for region 4
Pb = region4.Pb #[Mpa]
Tb = region4.Tb #[K]

#### Support for region 1 ####
Ps1 = region1.Ps    #[MPa]
Ts1 = region1.Ts    #[K]
R   = region1.R     #[kJ / kg K]

def theta1(pi):
    """ Dimensionless form for the temperature 
        as a function of pressure"""
    b = (pi * Ps1 / Pb)**0.25

    E = b**2 + n4[2] * b + n4[5]
    F = n4[0] * b**2 + n4[3] * b + n4[6]
    G = n4[1] * b**2 + n4[4] * b + n4[7]
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    return (n4[9] + D - ((n4[9] + D)**2 - 4 * (n4[8] + n4[9] * D))**0.5) / 2

def theta1_pi(pi):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature 
        as a function of pressure w.r.t. pressure (pi)"""

    pib = pi * Ps1 / Pb

    E = pib**2 + n4[2] * pib + n4[5]
    F = n4[0] * pib**2 + n4[3] * pib + n4[6]
    G = n4[1] * pib**2 + n4[4] * pib + n4[7]
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    E_pib = (1 / 2) * pib**(-1/2) + (n4[2] / 4) * pib**(-3/4)
    F_pib = (n4[0] / 2) * pib**(-1/2) + (n4[3] / 4) * pib**(-3/4)
    G_pib = (n4[1] / 2) * pib**(-1/2) + (n4[4] / 4) * pib**(-3/4)    
    D_pib = 2 * G_pib / (-F - (F**2 - 4 * E * G)**0.5) - \
        2 * G / (-F - (F**2 - 4 * E * G)**0.5)**2 * \
        (-F_pib - (2 * F * F_pib - 4 * E_pib * G - 4 * E * G_pib) / (2 * (F**2 - 4 * E * G)**0.5))

    return  (D_pib - (2 * (n4[9] + D) * D_pib - 4 * n4[9] * D_pib) / ((n4[9] + D)**2 - 4 * (n4[8] + n4[9] * D))**0.5) * (Ps1 / (2 * Pb))

def gamma1_pi(pi):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (Ts1 / (Tb * theta1(pi)) - 1.222)**Ji
    return sum

def gamma1_pi_pi(pi):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += (ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (Ts1 / (Tb * thetaR(pi)) - 1.222)**Ji) + \
               (ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * ((Ts1 * theta1_pi(pi)) / (Tb * theta1(pi)**2)) * (Ts1 / (Tb * theta1(pi)) - 1.222)**(Ji - 1))
    return sum

def gamma1_tau(pi):
    """ --- Not Tested ---
        Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau); using backward relation for 
        temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ts1 / (Tb * theta1(pi)) - 1.222)**(Ji - 1)
    return sum

def gamma1_tau_pi(pi):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau) and pressure (pi); using backward 
        relation for temperature"""
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += (-ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (Ts1 / (Tb * theta1(pi)) - 1.222)**(Ji - 1)) + \
               (-ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * ((Ts1 * theta1_pi(pi)) / (Tb * theta1(pi)**2)) * (Ts1 / (Tb * theta1(pi)) - 1.222)**(Ji - 2))
    return sum

#### Insert the remaining backward relavent relations ####

def vf(P):
    """ Specific volume [m^3 / kg]"""
    pi = P / Ps1

    return theta1(pi) * gamma1_pi(pi) * R * Tb / (Ps1 * 10**6 / 1000)

def dvfdp(P):
    """ --- Not Tested ---
        Derivative of specific volume [m^3 m^3 / kg kJ]
        w.r.t. pressure"""
    pi = P / Ps1

    return (theta1_pi(pi) * gamma1_pi(pi) + theta1(pi) * gamma1_pi_pi(pi)) * R * Tb / (Ps1 * 10**6 / 1000)**2

def hf(P):
    """ --- Not Tested ---
        specific enthalpy [kJ / kg]
        """
    pi = P / Ps1

    return gamma1_tau(pi) * R * Ts1

def dhfdp(P):
    """ --- Not Tested ---
        Derivative of specific enthalpy [kJ m^3 / kg kJ]
        w.r.t. pressure"""
    pi = P / Ps1

    return gamma1_tau_pi(pi) * R * Ts1 / (Ps1 * 10**6 / 1000)

#### Support for region 2 ####
Ps2 = region2.Ps    #[MPa]
Ts2 = region2.Ts    #[K]

def theta2(pi):
    """ Dimensionless form for the temperature 
        as a function of pressure"""
    b = (pi * Ps2 / Pb)**0.25

    E = b**2 + n4[2] * b + n4[5]
    F = n4[0] * b**2 + n4[3] * b + n4[6]
    G = n4[1] * b**2 + n4[4] * b + n4[7]
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    return (n4[9] + D - ((n4[9] + D)**2 - 4 * (n4[8] + n4[9] * D))**0.5) / 2

def theta2_pi(pi):  
    """ --- Not Tested ---
        Derivative of dimensionless form for the temperature 
        as a function of pressure w.r.t. pressure (pi)"""

    pib = pi * Ps2 / Pb

    E = pib**2 + n4[2] * pib + n4[5]
    F = n4[0] * pib**2 + n4[3] * pib + n4[6]
    G = n4[1] * pib**2 + n4[4] * pib + n4[7]
    D = 2*G / (-F - (F**2 - 4 * E * G)**0.5)

    E_pib = (1 / 2) * pib**(-1/2) + (n4[2] / 4) * pib**(-3/4)
    F_pib = (n4[0] / 2) * pib**(-1/2) + (n4[3] / 4) * pib**(-3/4)
    G_pib = (n4[1] / 2) * pib**(-1/2) + (n4[4] / 4) * pib**(-3/4)    
    D_pib = 2 * G_pib / (-F - (F**2 - 4 * E * G)**0.5) - \
        2 * G / (-F - (F**2 - 4 * E * G)**0.5)**2 * \
        (-F_pib - (2 * F * F_pib - 4 * E_pib * G - 4 * E * G_pib) / (2 * (F**2 - 4 * E * G)**0.5))

    return  (D_pib - (2 * (n4[9] + D) * D_pib - 4 * n4[9] * D_pib) / ((n4[9] + D)**2 - 4 * (n4[8] + n4[9] * D))**0.5) * (Ps2 / (2 * Pb))

def gamma2_pi(pi):
    """ Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = 1 / pi

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * Ii * pi**(Ii - 1) * ((Ts2 / (Tb * theta2(pi))) - 0.5)**Ji
    return sum

def gamma2_pi_pi(pi):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        pressure (pi); using backward relation for 
        temperature"""
    sum = -1 / pi**2

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += (ni * Ii * (Ii - 1) * pi**(Ii - 2) * ((Ts2 / (Tb * theta2(pi))) - 0.5)**Ji) + \
              (-ni * Ii * pi**(Ii - 1) * Ji * ((Ts2 * theta2_pi(pi)) / (Tb * theta2(pi)**2)) * ((Ts2 / (Tb * theta2(pi))) - 0.5)**(Ji - 1))
    return sum

def gamma2_tau(pi):
    """ --- Not Tested ---
        Derivative of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau); using backward relation for 
        temperature"""
    sum = 0 
    for Ji, ni in zip(J0, n0):
        sum += ni * Ji * (Ts2 / (Tb * theta2(pi)))**(Ji - 1)

    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += ni * pi**Ii * Ji * ((Ts2 / (Tb * theta2(pi))) - 0.5)**(Ji - 1)
    return sum

def gamma2_tau_pi(pi):
    """ --- Not Tested ---
        Derivative (second) of dimensionless form for the 
        specific Gibbs free energy w.r.t. dimensionless 
        temperature (tau) and pressure (pi); using backward 
        relation for temperature"""
    sum = 0 
    for Ii, Ji, ni in zip(Ir, Jr, nr):
        sum += (ni * Ii * pi**(Ii - 1) * Ji * ((Ts2 / (Tb * theta2(pi))) - 0.5)**(Ji - 1)) + \
              (-ni * pi**Ii * Ji * (Ji - 1) * ((Ts2 * theta2_pi(pi)) / (Tb * theta2(pi)**2)) * ((Ts2 / (Tb * theta2(pi))) - 0.5)**(Ji - 2))
    return sum

#### Insert the remaining backward relavent relations ####

def vg(P):
    """ Specific volume [m^3 / kg]"""
    pi = P / Ps2

    return theta2(pi) * gamma2_pi(pi) * R * Tb / (Ps2 * 10**6 / 1000)

def dvgdp(P):
    """ --- Not Tested ---
        Derivative of specific volume [m^3 m^3 / kg kJ]
        w.r.t. pressure"""
    pi = P / Ps2

    return (theta2_pi(pi) * gamma2_pi(pi) + theta2(pi) * gamma2_pi_pi(pi)) * R * Tb / (Ps2 * 10**6 / 1000)**2

def hg(P):
    """ --- Not Tested ---
        specific enthalpy [kJ / kg]
        """
    pi = P / Ps2

    return gamma2_tau(pi) * R * Ts2

def dhgdp(P):
    """ --- Not Tested ---
        Derivative of specific enthalpy [kJ m^3 / kg kJ]
        w.r.t. pressure"""
    pi = P / Ps2

    return gamma2_tau_pi(pi) * R * Ts2 / (Ps2 * 10**6 / 1000)