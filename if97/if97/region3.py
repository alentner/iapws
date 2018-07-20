# Constants for region 3
n = [ 0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2]

# Non-dimenionalization for region 3
Ps = 1.00       #[Mpa]
Ts = 1.00       #[K]

def bnd23P(T):
    theta = T / Ts

    return Ps * (n[0] + n[1] * theta + n[2] * theta**2)

def bnd23T(P):
    pi = P / Ps

    return Ts * (n[3] + ((pi - n[4]) / n[2])**0.5)
