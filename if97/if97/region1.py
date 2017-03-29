# Constants for region 1
I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32]
J = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
n = [0.14632971213167, -0.84548187169114, -0.37563603672040E1, 0.33855169168385E1, -0.95791963387872, 0.15772038513228, -0.16616417199501E-1,  0.81214629983568E-3, 0.28319080123804E-3, -0.60706301565874E-3, -0.18990068218419E-1,
     -0.32529748770505E-1, -0.21841717175414E-1, -0.52838357969930E-4, -0.47184321073267E-3, -0.30001780793026E-3, 0.47661393906987E-4, -0.44141845330846E-5, -0.72694996297594E-15, -0.31679644845054E-4, -0.28270797985312E-5,
     -0.85205128120103E-9, -0.22425281908000E-5, -0.65171222895601E-6, -0.14341729937924E-12, -0.40516996860117E-6, -0.12734301741641E-8, -0.17424871230634E-9, -0.68762131295531E-18, 0.14478307828521E-19, 0.26335781662795E-22,
     -0.11947622640071E-22, 0.18228094581404E-23, -0.93537087292458E-25]

# Non-dimenionalization for region 1
Ps = 16.53      #[Mpa]
Ts = 1386       #[K]
R = 0.461526    #[KJ / Kg K]

def gamma(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * (tau - 1.222)**Ji
    return sum
def gamma_pi(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * (tau - 1.222)**Ji
    return sum
def gamma_pipi(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * Ii * (Ii - 1) * (7.1 - pi)**(Ii - 2) * (tau - 1.222)**Ji
    return sum
def gamma_tau(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (tau - 1.222)**(Ji - 1)
    return sum
def gamma_tautau(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += ni * (7.1 - pi)**Ii * Ji * (Ji - 1) * (tau - 1.222)**(Ji - 2)
    return sum
def gamma_pitau(pi, tau):
    sum = 0
    for Ii, Ji, ni in zip(I, J, n):
        sum += -ni * Ii * (7.1 - pi)**(Ii - 1) * Ji * (tau - 1.222)**(Ji - 1)
    return sum
def gibbs(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return gamma(pi, tau) * R * T
def v(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return pi * gamma_pi(pi, tau) * R * T / (P * 10**6 / 1000)
def h(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return tau * gamma_tau(pi, tau) * R * T
def u(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return (tau * gamma_tau(pi, tau) - pi * gamma_pi(pi, tau)) * R * T    
def s(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return (tau * gamma_tau(pi, tau) - gamma(pi, tau)) * R
def cp(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return -tau**2 * gamma_tautau(pi, tau) * R
def cv(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return (-tau**2 * gamma_tautau(pi, tau) + (gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / gamma_pipi(pi, tau)) * R
def w(P, T):
    pi = P / Ps
    tau = Ts / T
    R = 0.461526

    return (gamma_pi(pi, tau)**2 * R * T * 1000 / ((gamma_pi(pi, tau) - tau * gamma_pitau(pi, tau))**2 / (tau**2 * gamma_tautau(pi, tau)) - gamma_pipi(pi, tau)))**0.5