"""Provides supporting routines for region identification in IF97 specification"""

# type annotations
from __future__ import annotations

# internal libraries
from . import region1, region2, region3, region4, region5, auxiliary

###########################################################
#####          Region Identification Functions        #####
###########################################################

def region_P(P: float) -> int:
    """Identification (validate) of region from IF97 specification,
    using pressure only and assuming two-phase mixture.
    Reference: Figure (1) from R7-97(2012)"""
    Pbnd0 = region4.Pbnd0
    Pbnd1 = region4.Pnbd1
    assert Pbnd0 <= P <= Pbnd1, "Invalid pressure for two-phase state!" 
    return 4

def region_T(T: float) -> int:
    """Identification (validate) of region from IF97 specification,
    using temperature only and assuming two-phase mixture.
    Reference: Figure (1) from R7-97(2012)"""
    Tbnd0 = region4.Tbnd0
    Tbnd1 = region4.Tbnd1
    assert Tbnd0 <= T <= Tbnd1, "Invalid temperature for two-phase state!" 
    return 4

def region_PT(P: float, T: float) -> int:
    """Identification of region from IF97 specification,
    using pressure and temperature as primary varibles.
    Reference: Figure (1) from R7-97(2012)"""

    # Constant boundaries
    Pbnd0  = region2.Pbnd0
    Pbnd1  = region2.Pbnd1
    Tbnd0  = region2.Tbnd0
    Tbnd1  = region5.Tbnd1
    Pbnd50 = region5.Pbnd1
    Tbnd25 = region2.Tbnd1
    Tbnd13 = region1.Tbnd1
    Tbnd32 = auxiliary.Tbnd1

    # Functional boundaries, as constants
    Pbnd32 = region3.bnd23P(min(max(T, Tbnd13), Tbnd32))
    Pbnd12 = region4.satP(min(max(T, Tbnd0), Tbnd13))

    # Determine region based on pressure and temperature
    if (Pbnd0 <= P <= Pbnd1) and (Tbnd0 <= T <= Tbnd1):
        if T <= Tbnd13:
            region = 1 if P >= Pbnd12 else 2
        elif T <= Tbnd25:
            region = 3 if P >= Pbnd32 else 2
        else:
            region = 5 if P <= Pbnd50 else 0
    else:
        region = 0

    # Validate proper region identification
    assert (region != 0), "Water properties not available!"
    assert (region != 3), "Water properties (r3) not implemented!"
    assert (region != 5), "Water properties (r5) not implemented!"
    return region

def region_Ph(P: float, h: float) -> int:
    """Identification of region from IF97 specification,
    using pressure and enthalpy as primary variables.
    Reference: Figure (1) from R7-97(2012)"""

    # Supporting constants
    Tbnd0  = region2.Tbnd0
    Tbnd13 = region1.Tbnd13
    Tbnd25 = region2.Tbnd1

    # Constant boundaries
    Pbnd0  = region2.Pbnd0
    Pbnd13 = region4.satP(Tbnd13)
    Pbnd50 = region5.Pbnd1
    Pbnd1  = region2.Pbnd1 
    hbnd0  = region2.h(Pbnd0, Tbnd0)
    hbnd1  = region2.h(Pbnd1, Tbnd25) # region 5 not implemented
    
    # Functional boundaries, as constants
    Tbnd32 = auxiliary.bnd23T(min(max(P, auxiliary.Pbnd0), Pbnd1))
    hbnd13 = region1.h(P, Tbnd13)
    hbnd32 = region2.h(P, Tbnd32)
    hbnd25 = region2.h(P, Tbnd25)
    Tbnd4  = region4.satT(P)
    hbnd14 = region1.h(P, Tbnd4)
    hbnd42 = region2.h(P, Tbnd4)

    # Determine region based on pressure and enthalpy
    if (Pbnd0 <= P <= Pbnd1) and (hbnd0 <= h <= hbnd1):
        if P >= Pbnd13:
            if h <= hbnd32:
                region = 1 if h <= hbnd13 else 3
            elif h <= hbnd25
                region = 2
            else:
                region = 5 if P <= Pbnd50 else 0
        else:
            if h <= hbnd14:
                region = 1
            elif h <= hbnd25:
                region = 2 if h >= hbnd42 else 4
            else:
                region = 5 if P <= Pbnd50 else 0
    else:
        region = 0

    # Validate proper region identification
    assert (region != 0), "Water properties not avalable!"
    assert (region != 3), "Water properties not avalable!"
    assert (region != 5), "Water properties not avalable!"
    return region

def region_Ps(P: float, s: float) -> int:
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
