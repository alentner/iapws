import unittest
from iapws.if97 import region2, region3, region4, h2o
from iapws.critical import moody
import numpy
import matplotlib
from matplotlib import pyplot, cm, patches, lines

class test_ThermodynamicPlots(unittest.TestCase):

    def test_ThermodynamicDerivative_Plot(self):

        Pfg = numpy.logspace(-3, numpy.log10(h2o.satP(623.14)), 50)
        P   = numpy.logspace(-3, numpy.log10(50), 100)
        Ts  = numpy.linspace(273.16, 623.14, 50)
        Tl  = numpy.linspace(273.16, 1073.14, 50)
        Tn = [10, 40, 60, 100, 200]

        hf = [h2o.hf(n) for n in Pfg]
        hg = [h2o.hg(n) for n in Pfg]
        h  = [[h2o.h(n, 273.16 + i) for n in P] for i in Tn]

        vf = [h2o.vf(n) for n in Pfg]
        vg = [h2o.vg(n) for n in Pfg]
        v  = [[h2o.v(n, 273.16 + i) for n in P] for i in Tn]

        uf = [h2o.uf(n) for n in Pfg]
        ug = [h2o.ug(n) for n in Pfg]
        u  = [[h2o.u(n, 273.16 + i) for n in P] for i in Tn]

        Tsat = [h2o.satT(n) for n in Pfg]
        P1 = numpy.array([numpy.linspace(h2o.satP(n), 50, len(Ts)) for n in Ts])
        T1 = numpy.array([numpy.ones(len(Ts)) * Ts[n] for n in range(len(Ts))])
        h1 = numpy.array([[h2o.h(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])
        v1 = numpy.array([[h2o.v(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])
        u1 = numpy.array([[h2o.u(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])
        dhdp1 = numpy.array([[h2o.dhdP(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])
        dvdp1 = numpy.array([[h2o.dvdP(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])
        dudp1 = numpy.array([[h2o.dudP(P1[j, i], Ts[j], region=1) for i in range(len(P1))] for j in range(len(Ts))])

        P2 = numpy.array([numpy.linspace(10**-3, (h2o.satP(n) if n <= 623.15 else region3.bnd23P(n)), len(Ts)) for n in Tl])
        T2 = numpy.array([numpy.ones(len(Ts)) * n for n in Tl])
        h2 = numpy.array([[h2o.h(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        v2 = numpy.array([[h2o.v(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        u2 = numpy.array([[h2o.u(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dhdp2 = numpy.array([[h2o.dhdP(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dvdp2 = numpy.array([[h2o.dvdP(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dudp2 = numpy.array([[h2o.dudP(P2[j, i], Tl[j], region=2) for i in range(len(P2[0]))] for j in range(len(Tl))])

        P4 = numpy.array([numpy.ones(len(Pfg)) * Pfg[n] for n in range(len(Pfg))])
        h4 = numpy.array([numpy.linspace(h2o.hf(n), h2o.hg(n), len(Ts)) for n in Pfg])
        dhdp4 = numpy.array([[region4.dhdP_h(Pfg[j], h4[j, i]) for i in range(len(h4))] for j in range(len(Pfg))])
        v4 = numpy.array([numpy.linspace(h2o.vf(n), h2o.vg(n), len(Ts)) for n in Pfg])
        dvdp4 = numpy.array([[region4.dvdP_h(Pfg[j], v4[j, i]) for i in range(len(v4))] for j in range(len(Pfg))])
        u4 = numpy.array([numpy.linspace(h2o.uf(n), h2o.ug(n), len(Ts)) for n in Pfg])
        dudp4 = numpy.array([[region4.dudP_h(Pfg[j], u4[j, i]) for i in range(len(u4))] for j in range(len(Pfg))])

        lvlh = numpy.linspace(-0.150, 0.015, 22+1)
        lvlh4 = numpy.linspace(-2.5, 50, 22)
        lvlv = numpy.linspace(-500, 0, 21)
        lvlu = numpy.linspace(-0.150, 0.00, 16)

        matplotlib.rcParams.update({'font.size': 8})
        clrmp = cm.viridis.copy()
        clrmp.set_under(cm.viridis.colors[0])
        clrmp.set_over(cm.viridis.colors[-1])
        clrmp.set_bad(color='black')
        fig = pyplot.figure(figsize=(13,11))

        pyplot.subplot(321, facecolor='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dhdp1, levels=lvlh/15, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dhdp2[:], levels=lvlh, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific enthalpy w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/15th', style='italic')

        pyplot.subplot(322, facecolor='darkgrey')
        pyplot.loglog(hf, Pfg, 'k', hg, Pfg, 'k')
        [pyplot.loglog(h[:][i], P, 'b', linewidth=0.70) for i in range(1,5)]
        pyplot.contourf(h1, P1, dhdp1, levels=lvlh/15, cmap=clrmp, extend="both")
        pyplot.contourf(h2, P2, dhdp2, levels=lvlh, cmap=clrmp, extend="both")
        pyplot.colorbar(ticks=lvlh[::4])
        pyplot.contourf(h4[:], P4[:], dhdp4[:], levels=lvlh4, cmap=clrmp, extend="both")
        pyplot.colorbar(ticks=lvlh4[1::4])
        pyplot.title('                    -- (w/ enthalpy @ constant T lines)', fontsize=8)
        pyplot.xlabel('Specific Enthalpy [kJ / kg]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.xlim(5*10**1, 4*10**3)
        pyplot.ylim(10**-3, 50)
        pyplot.text(60, 10, 'scale 1/15th', style='italic')
        pyplot.text(150, 2e-3, 'scale left only', style='italic')

        pyplot.subplot(323, facecolor='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dvdp1, levels=lvlv/15e9, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dvdp2[:], levels=lvlv, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific volume w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/15e9', style='italic')

        pyplot.subplot(324, facecolor='darkgrey')
        pyplot.loglog(vf, Pfg, 'k', vg, Pfg, 'k')
        [pyplot.loglog(v[:][i], P, 'b', linewidth=0.70) for i in range(5)]
        pyplot.contourf(v1, P1, dvdp1, levels=lvlv/15e9, cmap=clrmp, extend="both")
        pyplot.contourf(v4[:], P4[:], dvdp4[:], levels=lvlv/50, cmap=clrmp, extend="both")
        pyplot.contourf(v2[:], P2[:], dvdp2[:], levels=lvlv, cmap=clrmp, extend="both")
        pyplot.colorbar(ticks=lvlv[::4])
        pyplot.title('                    -- (w/ volume @ constant T lines)', fontsize=8)
        pyplot.xlabel('Specific volume [m^3 / kg]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.xlim(1*10**-3, 5*10**2)
        pyplot.ylim(10**-3, 50)
        pyplot.annotate('scale 1/15e9', xy=(1.25e-3, 15), xycoords='data', xytext=(2.0e-3, 2.5),
                    textcoords='data', arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=0,angleB=90"))
        pyplot.text(2.0e-3, 3.0e-3, 'scale 1/50th', style='italic')

        pyplot.subplot(325, facecolor='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dudp1, levels=lvlu/20, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dudp2[:], levels=lvlu, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific internal energy w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/20th', style='italic')

        pyplot.subplot(326, facecolor='darkgrey')
        pyplot.loglog(uf, Pfg, 'k', ug, Pfg, 'k')
        [pyplot.loglog(u[:][i], P, 'b', linewidth=0.70) for i in range(1,5)]
        pyplot.contourf(u1, P1, dudp1, levels=lvlu/20, cmap=clrmp, extend="both")
        pyplot.contourf(u4[:], P4[:], dudp4[:], levels=lvlu*30, cmap=clrmp, extend="both")
        pyplot.contourf(u2[:], P2[:], dudp2[:], levels=lvlu, cmap=clrmp, extend="both")
        pyplot.colorbar(ticks=lvlu[::2])
        pyplot.title('                    --  (w/ int. energy @ constant T lines)', fontsize=8)
        pyplot.xlabel('Specific Internal Energy [kJ / kg]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.xlim(2*10**1, 3.5*10**3)
        pyplot.ylim(10**-3, 50)
        pyplot.text(23, 10, 'scale 1/20th', style='italic')
        pyplot.text(150, 2e-3, 'scale x30', style='italic')

        fig.tight_layout()
        pyplot.savefig("tests/PartialP", dpi=300)

    def test_ThermodynamicProperty_Plot(self):
        P0 = -3
        np = 500
        Pfg = numpy.logspace(P0, numpy.log10(h2o.satP(623.14)), np)
        Pb3 = numpy.logspace(numpy.log10(h2o.satP(623.14)), 2, np)
        Pb2 = numpy.logspace(P0, 2, np)
        Pb5 = numpy.logspace(P0, numpy.log10(50), np)
        PTn = numpy.logspace(P0, 2, np)
        Tn = [10, 50, 100, 150, 200, 250, 300]
        Phn = numpy.logspace(P0, 2, 2*np)
        hn = [250, 500, 750, 1000, 1250, 1500]
        Psn = numpy.logspace(P0, 2, np)
        sn = [5.5, 6.0, 6.5, 7, 7.5, 8, 8.5]
        Pvp = [3, 80, 3, 
               0.0035, 0.0035, 30, 
               0.1, 1, 10, 0.1, 1, 10]
        Pvh = [3, 80, 80, 
               0.001, 3, 3, 
               5, 5, 25, 
               40, 60, 60]
        Pvs = [3, 80, 80, 
               0.1, 0.1, 2.5, 
               8, 8, 90, 
               20, 80, 80]

        hf = [h2o.hf(i) for i in Pfg]
        hg = [h2o.hg(i) for i in Pfg]
        h13 = [h2o.h(i, 623.15) for i in Pb3]
        h23 = [h2o.h(i, h2o.region3.bnd23T(i)) for i in Pb3]
        h25 = [h2o.h(i, 1073.15) for i in Pb5]
        h20  = [h2o.h(i, 1073.15) for i in Pb2]
        h20a = [h2o.h(i, 1073.15) for i in Pb3]
        h20b = [h2o.h(i, 1073.15) for i in Pfg]
        hTn = [[h2o.h(i, j + 273.15) for i in PTn] for j in Tn]
        hsn = [[h2o.h_s(i, j) if j <= h2o.region2.s(i, 1073.15) else h2o.h_s(i, h2o.region2.s(i, 1073.15)) for i in Psn] for j in sn]
        hvp = [h2o.h(3, 300),          h2o.h(80, 300),          h2o.h(3, 500),
               h2o.h(0.0035, 300),     h2o.h(0.0035, 700),      h2o.h(30, 700),
               h2o.hf(0.1), h2o.hf(1), h2o.hf(10), h2o.hg(0.1), h2o.hg(1), h2o.hg(10)]
        hvh = [500,  500,  1500, 
               3000, 3000, 4000, 
               3500, 4000, 3500, 
               2700, 2700, 3200]
        hvs = [h2o.h_s(3, 0.5),   h2o.h_s(80, 0.5),  h2o.h_s(80, 3),
               h2o.h_s(0.1, 7.5), h2o.h_s(0.1, 8),   h2o.h_s(2.5, 8),
               h2o.h_s(8, 6),     h2o.h_s(8, 7.5),   h2o.h_s(90, 6),
               h2o.h_s(20, 5.75), h2o.h_s(80, 5.25), h2o.h_s(80, 5.75)]

        sf = [h2o.sf(i) for i in Pfg]
        sg = [h2o.sg(i) for i in Pfg]
        s13 = [h2o.s(i, 623.15) for i in Pb3]
        s23 = [h2o.s(i, h2o.region3.bnd23T(i)) for i in Pb3]
        s25 = [h2o.s(i, 1073.15) for i in Pb5]
        s20  = [h2o.s(i, 1073.15) for i in Pb2]
        s20a = [h2o.s(i, 1073.15) for i in Pb3]
        s20b = [h2o.s(i, 1073.15) for i in Pfg]
        sTn = [[h2o.s(i, j + 273.15) for i in PTn] for j in Tn]
        shn = [[h2o.s_h(i, j) for i in Phn] for j in hn]
        svp = [h2o.s(3, 300),          h2o.s(80, 300),          h2o.s(3, 500),
               h2o.s(0.0035, 300),     h2o.s(0.0035, 700),      h2o.s(30, 700),
               h2o.sf(0.1), h2o.hf(1), h2o.sf(10), h2o.sg(0.1), h2o.sg(1), h2o.sg(10)]
        svh = [h2o.s_h(3, 500),      h2o.s_h(80, 500),  h2o.s_h(80, 1500), 
               h2o.s_h(0.001, 3000), h2o.s_h(3, 3000),  h2o.s_h(3, 4000), 
               h2o.s_h(5, 3500),     h2o.s_h(5, 4000),  h2o.s_h(25, 3500), 
               h2o.s_h(40, 2700),    h2o.s_h(60, 2700), h2o.s_h(60, 3200)]
        svs = [0.5,  0.5,  3,
               7.5,  8,    8,
               6,    7.5,  6,
               5.75, 5.25, 5.75]

        vf = [h2o.vf(i) for i in Pfg]
        vg = [h2o.vg(i) for i in Pfg]
        v13 = [h2o.v(i, 623.15) for i in Pb3]
        v23 = [h2o.v(i, h2o.region3.bnd23T(i)) for i in Pb3]
        v25 = [h2o.v(i, 1073.15) for i in Pb5]
        v20  = [h2o.v(i, 1073.15) for i in Pb2]
        v20a = [h2o.v(i, 1073.15) for i in Pb3]
        v20b = [h2o.v(i, 1073.15) for i in Pfg]
        vTn = [[h2o.v(i, j + 273.15) for i in PTn] for j in Tn]
        vhn = [[h2o.v_h(i, j) for i in Phn] for j in hn]
        vsn = [[h2o.v_s(i, j) if j <= h2o.region2.s(i, 1073.15) else h2o.v_s(i, h2o.region2.s(i, 1073.15)) for i in Psn] for j in sn]
        vvp = [h2o.v(3, 300),          h2o.v(80, 300),          h2o.v(3, 500),
               h2o.v(0.0035, 300),     h2o.v(0.0035, 700),      h2o.v(30, 700),
               h2o.vf(0.1), h2o.vf(1), h2o.vf(10), h2o.vg(0.1), h2o.vg(1), h2o.vg(10)]
        vvh = [h2o.v_h(3, 500),      h2o.v_h(80, 500),  h2o.v_h(80, 1500), 
               h2o.v_h(0.001, 3000), h2o.v_h(3, 3000),  h2o.v_h(3, 4000), 
               h2o.v_h(5, 3500),     h2o.v_h(5, 4000),  h2o.v_h(25, 3500), 
               h2o.v_h(40, 2700),    h2o.v_h(60, 2700), h2o.v_h(60, 3200)]
        vvs = [h2o.v_s(3, 0.5),   h2o.v_s(80, 0.5),  h2o.v_s(80, 3),
               h2o.v_s(0.1, 7.5), h2o.v_s(0.1, 8),   h2o.v_s(2.5, 8),
               h2o.v_s(8, 6),     h2o.v_s(8, 7.5),   h2o.v_s(90, 6),
               h2o.v_s(20, 5.75), h2o.v_s(80, 5.25), h2o.v_s(80, 5.75)]

        uf = [h2o.uf(i) for i in Pfg]
        ug = [h2o.ug(i) for i in Pfg]
        u13 = [h2o.u(i, 623.15) for i in Pb3]
        u23 = [h2o.u(i, h2o.region3.bnd23T(i)) for i in Pb3]
        u25 = [h2o.u(i, 1073.15) for i in Pb5]
        u20  = [h2o.u(i, 1073.15) for i in Pb2]
        u20a = [h2o.u(i, 1073.15) for i in Pb3]
        u20b = [h2o.u(i, 1073.15) for i in Pfg]
        uTn = [[h2o.u(i, j + 273.15) for i in PTn] for j in Tn]
        uhn = [[h2o.u_h(i, j) for i in Phn] for j in hn]
        usn = [[h2o.u_s(i, j) if j <= h2o.region2.s(i, 1073.15) else h2o.u_s(i, h2o.region2.s(i, 1073.15)) for i in Psn] for j in sn]
        uvp = [h2o.u(3, 300),          h2o.u(80, 300),          h2o.u(3, 500),
               h2o.u(0.0035, 300),     h2o.u(0.0035, 700),      h2o.u(30, 700),
               h2o.uf(0.1), h2o.uf(1), h2o.uf(10), h2o.ug(0.1), h2o.ug(1), h2o.ug(10)]
        uvh = [h2o.u_h(3, 500),      h2o.u_h(80, 500),  h2o.u_h(80, 1500), 
               h2o.u_h(0.001, 3000), h2o.u_h(3, 3000),  h2o.u_h(3, 4000), 
               h2o.u_h(5, 3500),     h2o.u_h(5, 4000),  h2o.u_h(25, 3500), 
               h2o.u_h(40, 2700),    h2o.u_h(60, 2700), h2o.u_h(60, 3200)]
        uvs = [h2o.u_s(3, 0.5),   h2o.u_s(80, 0.5),  h2o.u_s(80, 3),
               h2o.u_s(0.1, 7.5), h2o.u_s(0.1, 8),   h2o.u_s(2.5, 8),
               h2o.u_s(8, 6),     h2o.u_s(8, 7.5),   h2o.u_s(90, 6),
               h2o.u_s(20, 5.75), h2o.u_s(80, 5.25), h2o.u_s(80, 5.75)]

        colormap = cm.Set1.colors
        clr = [colormap[i] for i in [0, 1, 8, 2, 8, 3, 3, 6]]
        matplotlib.rcParams.update({'font.size': 8})
        matplotlib.rcParams.update({'lines.linewidth': 0.7})
        fig = pyplot.figure(figsize=(15, 11))
        pyplot.suptitle('Thermodynamic Properties of Water and Steam [1997 release of IAPWS]')

        ax1 = pyplot.subplot(221)
        pyplot.ylim(10**P0, 100)
        pyplot.xlim(0, 4500)
        pyplot.ylabel('Pressure [MPa]')
        pyplot.xlabel('Specific Enthalpy [kJ / kg]')
        pyplot.grid(linestyle="-", linewidth=0.4, color='.25', zorder=-10, alpha=0.2)
        pyplot.fill_betweenx(Pfg, 0.0, hf, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.fill_betweenx(Pb3, 0.0, h13, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.semilogy(h20, Pb2, color=clr[1], zorder=4)
        pyplot.fill_betweenx(Pfg, hg, h20b, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.fill_betweenx(Pb3, h23, h20a, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.semilogy(h13, Pb3, color=clr[2], zorder=2)
        pyplot.semilogy(h23, Pb3, color=clr[2], zorder=2)
        pyplot.fill_betweenx(Pb3, h13, h23, alpha=0.2, facecolor=clr[2], zorder=1)
        pyplot.semilogy(hf, Pfg, color=clr[3], zorder=2)
        pyplot.semilogy(hg, Pfg, color=clr[3], zorder=2)
        pyplot.fill_betweenx(Pfg, hf, hg, alpha=0.2, facecolor=clr[3], zorder=1)
        pyplot.semilogy(h25, Pb5, color=clr[4], zorder=2)
        pyplot.fill_betweenx(Pb5, h25, 5000, alpha=0.2, facecolor=clr[4], zorder=1)
        [pyplot.semilogy(hTn[j], PTn, color=clr[5], linestyle='--', zorder=3) for j in range(len(Tn))]
        [pyplot.semilogy(hsn[j], Psn, color=clr[6], linestyle=':', zorder=3) for j in range(len(sn))]
        pyplot.scatter(hvp[0:3],  Pvp[0:3],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(hvh[0:3],  Pvh[0:3],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(hvs[0:3],  Pvs[0:3],  facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(hvp[3:6],  Pvp[3:6],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(hvh[3:6],  Pvh[3:6],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(hvp[6:9],  Pvp[6:9],  s=5, facecolors=clr[7], marker='>', zorder=5)
        pyplot.scatter(hvp[9:12], Pvp[9:12], s=5, facecolors=clr[7], marker='<', zorder=5)

        ax2 = pyplot.subplot(222, sharey=ax1)
        pyplot.ylim(10**P0, 100)
        pyplot.xlim(0, 12)
        pyplot.xlabel('Specific Entropy [kJ / kg K]')
        pyplot.grid(linestyle="-", linewidth=0.4, color='.25', zorder=-10, alpha=0.2)
        pyplot.fill_betweenx(Pfg, 0.0, sf, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.fill_betweenx(Pb3, 0.0, s13, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.semilogy(s20, Pb2, color=clr[1], zorder=4)
        pyplot.fill_betweenx(Pfg, sg, s20b, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.fill_betweenx(Pb3, s23, s20a, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.semilogy(s13, Pb3, color=clr[2], zorder=2)
        pyplot.semilogy(s23, Pb3, color=clr[2], zorder=2)
        pyplot.fill_betweenx(Pb3, s13, s23, alpha=0.2, facecolor=clr[2], zorder=1)
        pyplot.semilogy(sf, Pfg, color=clr[3], zorder=2)
        pyplot.semilogy(sg, Pfg, color=clr[3], zorder=2)
        pyplot.fill_betweenx(Pfg, sf, sg, alpha=0.2, facecolor=clr[3], zorder=1)
        pyplot.semilogy(s25, Pb5, color=clr[4], zorder=2)
        pyplot.fill_betweenx(Pb5, s25, 15, alpha=0.2, facecolor=clr[4], zorder=1)
        [pyplot.semilogy(sTn[j], PTn, color=clr[5], linestyle='--', zorder=3) for j in range(len(Tn))]
        [pyplot.semilogy(shn[j], Phn, color=clr[6], linestyle='-.', zorder=3) for j in range(len(hn))]
        pyplot.scatter(svp[0:3],  Pvp[0:3],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(svh[0:3],  Pvh[0:3],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(svs[0:3],  Pvs[0:3],  facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(svp[3:6],  Pvp[3:6],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(svh[3:6],  Pvh[3:6],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(svp[6:9],  Pvp[6:9],  s=5, facecolors=clr[7], marker='>', zorder=5)
        pyplot.scatter(svp[9:12], Pvp[9:12], s=5, facecolors=clr[7], marker='<', zorder=5)
        pyplot.setp(ax2.get_yticklabels(), visible=False)

        ax3 = pyplot.subplot(223)
        pyplot.ylim(10**P0, 100)
        pyplot.xlim(20, 5000)
        pyplot.ylabel('Pressure [MPa]')
        pyplot.xlabel('Specific internal energy [kJ / kg]')
        pyplot.grid(linestyle="-", linewidth=0.4, color='.25', zorder=-10, alpha=0.2)
        pyplot.fill_betweenx(Pfg, 0.0, uf, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.fill_betweenx(Pb3, 0.0, u13, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.loglog(u20, Pb2, color=clr[1], zorder=4)
        pyplot.fill_betweenx(Pfg, ug, u20b, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.fill_betweenx(Pb3, u23, u20a, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.loglog(u13, Pb3, color=clr[2], zorder=2)
        pyplot.loglog(u23, Pb3, color=clr[2], zorder=2)
        pyplot.fill_betweenx(Pb3, u13, u23, alpha=0.2, facecolor=clr[2], zorder=1)
        pyplot.loglog(uf, Pfg, color=clr[3], zorder=2)
        pyplot.loglog(ug, Pfg, color=clr[3], zorder=2)
        pyplot.fill_betweenx(Pfg, uf, ug, alpha=0.2, facecolor=clr[3], zorder=1)
        pyplot.loglog(u25, Pb5, color=clr[4], zorder=2)
        pyplot.fill_betweenx(Pb5, u25, 10000, alpha=0.2, facecolor=clr[4], zorder=1)
        [pyplot.loglog(uTn[j], PTn, color=clr[5], linestyle='--', zorder=3) for j in range(len(Tn))]
        [pyplot.loglog(uhn[j], Phn, color=clr[5], linestyle='-.', zorder=3) for j in range(len(hn))]
        [pyplot.loglog(usn[j], Psn, color=clr[6], linestyle=':', zorder=3) for j in range(len(sn))]
        pyplot.scatter(uvp[0:3],  Pvp[0:3],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(uvh[0:3],  Pvh[0:3],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(uvs[0:3],  Pvs[0:3],  facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(uvp[3:6],  Pvp[3:6],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(uvh[3:6],  Pvh[3:6],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(uvs[3:12], Pvs[3:12], facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(uvp[6:9],  Pvp[6:9],  s=5, facecolors=clr[7], marker='>', zorder=5)
        pyplot.scatter(uvp[9:12], Pvp[9:12], s=5, facecolors=clr[7], marker='<', zorder=5)

        ax4 = pyplot.subplot(224, sharey=ax3)
        pyplot.ylim(10**P0, 100)
        pyplot.xlim(0.001, 1000)
        pyplot.xlabel('Specific volume [m**3 / kg]')
        pyplot.grid(linestyle="-", linewidth=0.4, color='.25', zorder=-10, alpha=0.2)
        pyplot.fill_betweenx(Pfg, 0.0, vf, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.fill_betweenx(Pb3, 0.0, v13, alpha=0.2, facecolor=clr[0], zorder=1)
        pyplot.loglog(v20, Pb2, color=clr[1], zorder=4)
        pyplot.fill_betweenx(Pfg, vg, v20b, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.fill_betweenx(Pb3, v23, v20a, alpha=0.2, facecolor=clr[1], zorder=1)
        pyplot.loglog(v13, Pb3, color=clr[2], zorder=2)
        pyplot.loglog(v23, Pb3, color=clr[2], zorder=2)
        pyplot.fill_betweenx(Pb3, v13, v23, alpha=0.2, facecolor=clr[2], zorder=1)
        pyplot.loglog(vf, Pfg, color=clr[3], zorder=2)
        pyplot.loglog(vg, Pfg, color=clr[3], zorder=2)
        pyplot.fill_betweenx(Pfg, vf, vg, alpha=0.2, facecolor=clr[3], zorder=1)
        pyplot.loglog(v25, Pb5, color=clr[4], zorder=2)
        pyplot.fill_betweenx(Pb5, v25, 1000, alpha=0.2, facecolor=clr[4], zorder=1)
        [pyplot.loglog(vTn[j], PTn, color=clr[5], linestyle='--', zorder=3) for j in range(len(Tn))]
        [pyplot.loglog(vhn[j], Phn, color=clr[5], linestyle='-.', zorder=3) for j in range(len(hn))]
        [pyplot.loglog(vsn[j], Psn, color=clr[6], linestyle=':', zorder=3) for j in range(len(sn))]
        pyplot.scatter(vvp[0:3],  Pvp[0:3],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(vvh[0:3],  Pvh[0:3],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(vvs[0:3],  Pvs[0:3],  facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(vvp[3:6],  Pvp[3:6],  s=5, facecolors=clr[7], marker='o', zorder=5)
        pyplot.scatter(vvh[3:6],  Pvh[3:6],  facecolors=clr[7], marker='+', zorder=5)
        pyplot.scatter(vvs[3:12], Pvs[3:12], facecolors=clr[7], marker='2', zorder=5)
        pyplot.scatter(vvp[6:9],  Pvp[6:9],  s=5, facecolors=clr[7], marker='>', zorder=5)
        pyplot.scatter(vvp[9:12], Pvp[9:12], s=5, facecolors=clr[7], marker='<', zorder=5)
        pyplot.setp(ax4.get_yticklabels(), visible=False)

        r1 = patches.Patch(color=clr[0], alpha=0.2, label='Liquid Region         (1)')
        r2 = patches.Patch(color=clr[1], alpha=0.2, label='Vapor Region         (2)')
        r3 = patches.Patch(color=clr[2], alpha=0.2, label='Not Implemented   (3, 5)')
        r4 = patches.Patch(color=clr[3], alpha=0.2, label='Two-phase Region  (4)')
        cT = lines.Line2D([], [], color=clr[5], linestyle='--', label='Lines of Constant Temperature')
        ch = lines.Line2D([], [], color=clr[5], linestyle='-.', label='Lines of Constant Enthalpy')
        cs = lines.Line2D([], [], color=clr[6], linestyle=':',  label='Lines of Constant Entropy')
        vT = lines.Line2D([], [], color=clr[7], linestyle='none', marker='o', markersize=5, label='f(P, T) validation Points')
        vh = lines.Line2D([], [], color=clr[7], linestyle='none', marker='+', markersize=5, label='f(P, h) validation Points')
        vs = lines.Line2D([], [], color=clr[7], linestyle='none', marker='2', markersize=5, label='f(P, s) validation Points')
        vp = lines.Line2D([], [], color=clr[7], linestyle='none', marker='>', markersize=5, label='Saturation line validation Points')
        pyplot.legend(handles=[r1, r2, r3, r4, cT, ch, cs, vT, vh, vs, vp])

        fig.tight_layout()
        pyplot.subplots_adjust(top=0.95)
        pyplot.savefig("tests/Properties", dpi=300)

    def test_CriticalFlowModuleHEM_Plot(self):
        P0 = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 
              10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        h0 = numpy.linspace(250, 2750, 50)

        g = [[moody.Gcritical(j, i) for i in h0] for j in P0[:15]]

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
        pyplot.savefig("tests/Moody", dpi=300)

if __name__ == '__main__':
    unittest.main()
