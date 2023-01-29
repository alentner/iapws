

# type annotations
from __future__ import annotations

# internal libraries
from . import region1, region2, region3, region4, region5
from . import unit
from .support import _english, _output

###########################################################
#####          Region Identification Functions        #####
###########################################################

@_english((unit.P, unit.T), (_output, ))
def region(P: float, T: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and temperature as primary varibles"""

    # Constant boundaries
    Pbnd0  = region2.Pbnd0
    Pbnd1  = region1.Pbnd1
    Tbnd0  = region1.Tbnd0
    Tbnd1  = region5.Tbnd1
    Pbnd50 = region5.Pbnd1
    Tbnd25 = region2.Tbnd1
    Tbnd13 = region1.Tbnd1
    Tbnd32 = auxiliary.Tbnd1

    # Functional boundaries, as constants
    Pbnd32 = region3.bnd23P(min(max(T, Tbnd13), Tbnd32))
    Pbnd12 = satP(min(max(T, Tbnd0), Tbnd13))

    # Determine region based on pressure and temperature
    if (Pbnd0 <= P <= Pbnd1) and (Tbnd0 <= T <= Tbnd1):
        if T <= Tbnd13:
            region = 1 if P >= Pbnd12 else 2
        elif T <= Tbnd25:
            region = 5 if P >= Pbnd32 else 2
        elif P <= Pbnd50:
            region = 5
        else:
            region = 0
    else:
        region = 0

    # Validate proper region identification
    assert (region != 0), "Water properties not avalable!"
    return region

@_english((unit.P, unit.h), (_output, ))
def region_h(P: float, h: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # Supporting boundaries
    Tbnd01 = region1.Tbnd01
    Pbnd4  = satP(Tbnd01)
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13
    Tbnd32 = region3.bnd23T(min(max(P, 16.5292), 100.0))
    Tbnd4  = satT(P)

    # Enthalpy- pressure boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1 
    hbnd01 = region1.h(Pbnd4, Tbnd01)
    hbnd25 = region2.h(Pbnd0, Tbnd25)
    Pbndh1 = satP(Tbnd13)
    hbnd13 = region1.h(P, Tbnd13)
    hbnd32 = region2.h(P, Tbnd32)
    hbnd14 = region1.h(P, Tbnd4)
    hbnd42 = region2.h(P, Tbnd4)

    region = 0

    # Only region 1,2,4 via P,h relations implemented
    if (P >= Pbnd0) and (h >= hbnd01) and (P <= Pbnd1) and (h <= hbnd25):
        if (P >= Pbndh1):
            if (h <= hbnd13):
                region = 1
            elif (h >= hbnd32):
                region = 2
            else:
                region = 0
        else:
            if (h <= hbnd14):
                region = 1
            elif (h >= hbnd42):
                region = 2
            else:
                region = 4

    assert (region != 0), "Water properties not avalable!"
    return region

@_english((unit.P, unit.s), (_output, ))
def idRegion_s(P: float, s: float, /, *, english: bool = False) -> int:
    """Identification of region from IF97 specification
    using pressure and enthalpy as primary variables"""

    # Supporting boundaries
    Tbnd01 = region1.Tbnd01
    Pbnd4  = satP(Tbnd01)
    Tbnd25 = region2.Tbnd25
    Tbnd13 = region1.Tbnd13
    Tbnd32 = region3.bnd23T(min(max(P, 16.5292), 100.0))
    Tbnd4  = satT(P)

    # Enthalpy-pressure boundaries
    Pbnd0  = region1.Pbnd0
    Pbnd1  = region1.Pbnd1 
    sbnd01 = region1.s(P, Tbnd01)
    sbnd25 = region2.s(P, Tbnd25)
    Pbndh1 = satP(Tbnd13)
    sbnd13 = region1.s(P, Tbnd13)
    sbnd32 = region2.s(P, Tbnd32)
    sbnd14 = region1.s(P, Tbnd4)
    sbnd42 = region2.s(P, Tbnd4)

    region = 0

    # Only region 1,2,4 via P,s relations implemented
    if (P >= Pbnd0) and (s >= sbnd01) and (P <= Pbnd1) and (s <= sbnd25):
        if (P >= Pbndh1):
            if (s <= sbnd13):
                region = 1
            elif (s >= sbnd32):
                region = 2
            else:
                region = 0
        else:
            if (s <= sbnd14):
                region = 1
            elif (s >= sbnd42):
                region = 2
            else:
                region = 4

    assert (region != 0), "Water properties not avalable!"
    return region

