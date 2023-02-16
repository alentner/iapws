
# type annotations
from __future__ import annotations

# system libraries
from typing import NamedTuple
from functools import partial

# internal libraries
from . import region1, region2, region3, region4, region5
from . import backward1, backward2, backward3, backward5, iterate
from . import gibbs, helmholtz
from . import identify

# static analysis
from typing import Optional

class State:
    P: float
    T: float
    x: Optional[float]
    english: bool
    form: str
    region: int

    def __init__(self, P: Optional[float] = None, T: Optional[float] = None, *,
                 h: Optional[float] = None, s: Optional[float] = None,
                 form: Optional[str] = None, region: Optional[int] = None, english: bool = False):
        """Initialize a state based on provided properties, units, and possible known region.""" 
        form = get_formulation(P, T, h=h, s=s, form=form)
        metric = get_metric(P, T, h=h, s=s, form=form, english=english)
        region = get_region(form=form, region=region, **metric)
        self._P, self._T = get_forward(form=form, region=region, **metric)
        seld._h, self._s = unit.h(h) if english else h, unit.s(s) if english else s
        self._equation = {1: region1, 2: region2, 3: region3, 4: region4, 5: region5}[region]
        self._property = {None: self._z,
                          'fg': self._zfg,
                          'f' : self._zf,
                          'g' : self._zg}
        self._partials = {None: {1: partial(gibbs.dzdx_y, equation=self._equation),
                                 2: partial(gibbs.dzdx_y, equation=self._equation),
                                 3: partial(helmholtz.dzdx_y, equation=self._equation),
                                 4: self._dzdxy,
                                 5: partial(gibbs.dzdxy, equation=self._equation)}[region],
                          'fg': self._dzfgdxy,
                          'f' : self._dzfdxy,
                          'g' : self._dzgdxy} 
        self.english = english
        self.form = form
        self.region = region

    @property
    def P(self) -> float:
        """Return the Pressure of the thermodynamic state."""
        return unit.P(self._P, english=False) if self.english else self._P

    @property
    def T(self) -> float:
        """Return the Temperature of the thermodynamic state."""
        return unit.T(self._T, english=False) if self.english else self._T

    @property
    def x(self) -> float:
        """Return the Equilibrium Quality of the thermodynamic state."""
        if self.region == 4:
            if self._h is not None:
                return (self._h - self._zf('h')) / self._zfg('h')
            elif self._h is not None:
                if self.english: s = unit.s(s)
                return (self._s - self._zf('s')) / self._zfg('s')
            assert False, "Two properties required to uniquely define state!"
        assert False, "Quality undefined for single phase fluid!"
        return None

    def __getattr__(self, name: str) -> float:
        """Determine the specified fluid property based on the name provided."""
        base, suffix, partial, constant = get_details(name)
        if constant is None:
            constant = self.form.strip(partial)
        if partial:
            value = self._partials[suffix](base, partial, constant)
            if self.english:
                value = getattr(unit, f"d{base}")(value, english=False)
                value = getattr(unit, f"d_d{partial}")(value, english=False)
        else:
            value = self._property[suffix](base) 
            if self.english:
                value = getattr(unit, base)(value, english=False)
        return value
    
    def _zf(self, name: str) -> Optional[float]:
        """Fluid property of saturated liquid."""
        if self.region == 4:
            equation = region1 if self.T <= auxiliary.Tbnd0 else region3
            return getattr(equation, name)(self.P, self.T)
        assert False, "Saturated liquid properties undefined for single phase fluid!"
        return None
    
    def _zg(self, name: str) -> Optional[float]:
        """Fluid property of saturated vapor."""
        if self.region == 4:
            equation = region2 if self.T <= auxiliary.Tbnd0 else region3
            return getattr(equation, name)(self.P, self.T)
        assert False, "Saturated vapor properties undefined for single phase fluid!"
        return None

    def _zfg(self, name: str) -> Optional[float]:
        """Fluid property rise between saturated liquid and vapor."""
        if self.region == 4:
            return self._zg(name) - self._zf(name)
        assert False, "Saturation rise properties undefined for single phase fluid!"
        return None

    def _z(self, name: str) -> float:
        """Fluid property of single phase fluid."""
        if self.region == 4:
            x = self.x
            return x * self._zg(name) + (1 - x) * self._zf(name)
        else:
            return getattr(self._equation, name)(self.P, self.T)

def get_details(function: str) -> tuple[str, Optional[str], Optional[str]]:
    """Determine key information about the desired property function."""
    phases = ('fg', 'f', 'g')
    basic, constant = basic.split('_') if '_' in function else function, None
    _, basic, partial = basic.split('d') if 'd' in function else None, basic, None
    suffix = next(filter(basic[1:].endswith, phases), None)
    basic = basic[:-len(suffix)] if suffix is not None else basic
    return basic, suffix, partial, constant

def get_formulation(P: Optional[float] = None, T: Optional[float] = None, *,
                    h: Optional[float] = None, s: Optional[float] = None,
                    form: Optional[str] = None) -> str:
    """Identify (validate) the formulation based on the provided properties."""
    if form is not None:
        assert form in {'P', 'T', 'PT', 'Ph', 'Ps', 'hs', 'Th', 'Ts'
                }, f"Provided formulation, {form}, is not available!"
        return form
    if P is not None and T is not None:
        form = 'PT'
    elif P is not None:
        valid = (h is not None) or (s is not None)
        form = 'P' if not valid else 'Ph' if h is not None else 'Ps'
    elif T is not None:
        valid = (h is not None) or (s is not None)
        form = 'T' if not valid else 'Th' if h is not None else 'Ts'
    else:
        valid = (h is not None) and (s is not None)
        assert valid, "Two properties required to uniquely define state!" 
        form = 'hs'
    return form

def get_forward(P: Optional[float] = None, T: Optional[float] = None, *,
                h: Optional[float] = None, s: Optional[float] = None,
                form: Optional[str] = None, region: Optional[int] = None) -> tuple[float, float]:
    """Determine (ensure) the forward formulation using backward functions;
    this method requires provided properties to be in metric units."""
    form = get_formulation(P, T, h=h, s=s, form=form)
    region = get_region(P, T, h=h, s=s, form=form, region=region)
    backward = {1: backward1, 2: backward2, 3: backward3, 4: backward4, 5: backward5}[region]
    return {'P' : (P, region4.satT(P)),
            'T' : (region4.satP(T), T),
            'PT': (P, T),
            'Ph': (P, backward.T_Ph(P, h)),
            'Ps': (P, backward.T_Ps(P, s)),
            'Th': (iterate.P_Th(T, h), T),
            'Ts': (iterate.P_Ts(T, s), T),
            'hs': (backward.P_hs(h, s), backward.T_hs(h, s))
            }[form]

def get_metric(P: Optional[float] = None, T: Optional[float] = None, *,
               h: Optional[float] = None, s: Optional[float] = None,
               form: Optional[str] = None, english: bool = False) -> dict[str, float]:
    """Determine (ensure) the metric units based formulation using appropriate conversions."""
    form = get_formulation(P, T, h=h, s=s, form=form)
    prop = {'P' : (unit.P(P),          ) if english else (P,  ),
            'T' : (unit.T(T),          ) if english else (T,  ),
            'PT': (unit.P(P), unit.T(T)) if english else (P, T),
            'Ph': (unit.P(P), unit.h(h)) if english else (P, h),
            'Ps': (unit.P(P), unit.s(s)) if english else (P, s),
            'Th': (unit.T(T), unit.h(h)) if english else (T, h),
            'Ts': (unit.T(T), unit.s(s)) if english else (T, s),
            'hs': (unit.h(h), unit.s(s)) if english else (h, s)
            }[form]
    return {f: p for f, p in zip(form, prop)}

def get_region(P: Optional[float] = None, T: Optional[float] = None, *
               h: Optional[float] = None, s: Optional[float] = None,
               form: Optional[str] = None, region: Optional[int] = None) -> int:
    """Identification of region (IF97 specification) based on provided properties;
    this method requires provided properties to be in metric units."""
    if region is not None:
        assert region in {1, 2, 3, 4, 5}, f"Region, {region}, is not defined in IF97!"
        return region
    form = get_formulation(P, T, h=h, s=s, form=form)
    return {'P' : identify.region_P(P),
            'T' : identify.region_T(T),
            'PT': identify.region_PT(P, T), 
            'Ph': identify.region_Ph(P, h),
            'Ps': identify.region_Ps(P, s),
            'hs': identify.region_hs(h, s),
            'Th': identify.region_Th(T, h),
            'Ts': identify.region_Ts(T, s),
            }[form]
