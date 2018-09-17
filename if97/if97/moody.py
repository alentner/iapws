from if97 import h2o
from scipy.optimize import brentq as root
from matplotlib import pyplot

def dGdPc(P, h0, s0):
    h = h2o.h_s(P, s0)

    if (h >= h0):
        return 100.0

    dhdp = h2o.dhdP_s(P, s0)
    v = h2o.v_s(P, s0)
    dvdp = h2o.dvdP_s(P, s0)
    return (dvdp / v**2) * (2 * (h0 - h))**0.5 + dhdp / (v * (2 * (h0 - h))**0.5)
def Gc(P0, h0):
    
    s0 = h2o.s_h(P0, h0)
    Pc = root(dGdPc, 611.213e-6, P0 - 1e-6, (h0, s0))
    
    h = h2o.h_s(Pc, s0)
    v = h2o.v_s(Pc, s0)
    
    return 1000**0.5 * (2 * (h0 - h))**0.5 / v


import numpy
P0 = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 
      10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
h0 = numpy.linspace(250, 2750, 50)

g = [[Gc(j, i) for i in h0] for j in P0[:15]]

#pyplot.xkcd()
pyplot.figure(figsize=(11, 7))
[pyplot.semilogy(h0, g[j], markersize=5, marker='.', label=str(P0[j])) for j in range(len(P0[:15]))]
pyplot.ylim(10, 10**6)
pyplot.ylabel('Critical Mass Flux [kg / s m**2]')
pyplot.xlim(0, 3000)
pyplot.xlabel('Stagnation Enthalpy [kJ / kg]')
pyplot.legend(loc='lower left', ncol=5)
pyplot.text(750, 8e1, 'Stagnation pressures [Mpa]', backgroundcolor="white", ha='center', va='top', color='black')
pyplot.grid(linestyle="--", linewidth=0.5, color='.25', zorder=-10)
pyplot.title('Critical Flow')
pyplot.savefig("if97/__testout__/test_moody", dpi=300)
pyplot.show()