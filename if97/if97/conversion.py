# Temperature conversions
def temp_C_F(T):

    return 180 / 100 * T + 32
def temp_F_C(T):

    return (T - 32) * 100 / 180
def temp_C_K(T):

    return T + 273.15
def temp_K_C(T):

    return T - 273.15
def temp_F_R(T):

    return T + 459.67
def temp_R_F(T):

    return T - 459.67

# Pressure conversions
def press_MPa_psia(P):

    return P * 145.0377439
def press_psia_MPa(P):

    return P / 145.0377439