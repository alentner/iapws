from if97 import region1, region2, region3, region4

#### identify if97 region ####
def idRegion(P, T):
    if (P >= 0) and (T >= 0) and (P <= region1.Pbnd1):
        if (T <= region1.Tbnd13) and (P >= region4.satP(T)):
            return 1
        elif (T < region1.Tbnd13) or (P <= region3.bnd23P(T)):
            return 2
        else:
            return 3
    return 0

#### water properties ####
def g(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.g(P, T)
    elif region is 2:
        return region2.g(P, T)
    else:
        return 0.000

def v(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.v(P, T)
    elif region is 2:
        return region2.v(P, T)
    else:
        return 0.000

def u(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.u(P, T)
    elif region is 2:
        return region2.u(P, T)
    else:
        return 0.000

def s(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.s(P, T)
    elif region is 2:
        return region2.s(P, T)
    else:
        return 0.000

def h(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.h(P, T)
    elif region is 2:
        return region2.h(P, T)
    else:
        return 0.000

def cp(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.cp(P, T)
    elif region is 2:
        return region2.cp(P, T)
    else:
        return 0.000

def cv(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.cv(P, T)
    elif region is 2:
        return region2.cv(P, T)
    else:
        return 0.000

def w(P, T, region = 0):
    if region is 0:
        region = idRegion(P, T)

    if region is 1:
        return region1.w(P, T)
    elif region is 2:
        return region2.w(P, T)
    else:
        return 0.000

#### Saturated liquid properties ####
def gf(P):
    return region2.g(P, region4.satT(P))

def vf(P):
    return region2.v(P, region4.satT(P))

def uf(P):
    return region2.u(P, region4.satT(P))

def sf(P):
    return region2.s(P, region4.satT(P))

def hf(P):
    return region2.h(P, region4.satT(P))

def cpf(P):
    return region2.cp(P, region4.satT(P))

def cvf(P):
    return region2.cv(P, region4.satT(P))

def wf(P):
    return region2.w(P, region4.satT(P))

#### Saturated vapor properties ####
def gg(P):
    return region2.g(P, region4.satT(P))

def vg(P):
    return region2.v(P, region4.satT(P))

def ug(P):
    return region2.u(P, region4.satT(P))

def sg(P):
    return region2.s(P, region4.satT(P))

def hg(P):
    return region2.h(P, region4.satT(P))

def cpg(P):
    return region2.cp(P, region4.satT(P))

def cvg(P):
    return region2.cv(P, region4.satT(P))

def wg(P):
    return region2.w(P, region4.satT(P))

#### delta saturation properties ####
def gfg(P):
    return gg(P) - gf(P)

def vfg(P):
    return vg(P) - vf(P)

def ufg(P):
    return ug(P) - uf(P)

def sfg(P):
    return sg(P) - sf(P)

def hfg(P):
    return hg(P) - hf(P)

def cpfg(P):
    return cpg(P) - cpf(P)

def cvfg(P):
    return cvg(P) - cvf(P)

def wfg(P):
    return wg(P) - wf(P)