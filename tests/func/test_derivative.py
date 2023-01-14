import unittest
from iapws.if97 import region1, region2, region3, region4, h2o
import numpy

class test_ThermodynamicDerivative(unittest.TestCase):

    def test_ThermodynamicDerivative_Region1(self):
        n = 100

        ## partial with respect to T
        Pn = [3, 80]
        Tn = [numpy.linspace(273.15 + 0.5, region4.satT(Pn[0]) - 0.5, n),
              numpy.linspace(273.15 + 0.5, 623.15 - 0.5, n)]
        for k, p in enumerate(Pn):
            T = Tn[k]

            dvdtn = [(region1.v(p, i + 0.1) - region1.v(p, i - 0.1)) / 0.2 for i in T]
            dudtn = [(region1.u(p, i + 0.1) - region1.u(p, i - 0.1)) / 0.2 for i in T]
            dhdtn = [(region1.h(p, i + 0.1) - region1.h(p, i - 0.1)) / 0.2 for i in T]
            dsdtn = [(region1.s(p, i + 0.1) - region1.s(p, i - 0.1)) / 0.2 for i in T]
            dgdtn = [(region1.g(p, i + 0.1) - region1.g(p, i - 0.1)) / 0.2 for i in T]
        
            dvdt  = [region1.dvdT(p, i) for i in T]
            dudt  = [region1.dudT(p, i) for i in T]
            dhdt  = [region1.dhdT(p, i) for i in T]
            dsdt  = [region1.dsdT(p, i) for i in T]
            dgdt  = [region1.dgdT(p, i) for i in T]

            self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state '+str(k+1)+', region 1!')
            
        ## partial with respect to P
        Tn = [300, 500]
        Pn = [numpy.linspace(region4.satP(Tn[0]) + 0.5, 100 - 0.5, n),
              numpy.linspace(region4.satP(Tn[1]) + 0.5, 100 - 0.5, n)]
        for k, T in enumerate(Tn):
            p = Pn[k]

            dvdpn = [(region1.v(i + 0.1, T) - region1.v(i - 0.1, T)) / (0.2 * 1e3) for i in p]
            dudpn = [(region1.u(i + 0.1, T) - region1.u(i - 0.1, T)) / (0.2 * 1e3) for i in p]
            dhdpn = [(region1.h(i + 0.1, T) - region1.h(i - 0.1, T)) / (0.2 * 1e3) for i in p]
            dsdpn = [(region1.s(i + 0.1, T) - region1.s(i - 0.1, T)) / (0.2 * 1e3) for i in p]
            dgdpn = [(region1.g(i + 0.1, T) - region1.g(i - 0.1, T)) / (0.2 * 1e3) for i in p]
        
            dvdp  = [region1.dvdP(i, T) for i in p]
            dudp  = [region1.dudP(i, T) for i in p]
            dhdp  = [region1.dhdP(i, T) for i in p]
            dsdp  = [region1.dsdP(i, T) for i in p]
            dgdp  = [region1.dgdP(i, T) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state '+str(k+1)+', region 1!')

    def test_ThermodynamicDerivative_Region2(self):
        n = 100

        ## partial with respect to T
        Pn = [0.0035, 10]
        Tn = [numpy.linspace(region4.satT(Pn[0]) + 0.5, 1073.15 - 0.5, n),
              numpy.linspace(region4.satT(Pn[1]) + 0.5, 1073.15 - 0.5, n)]
        for k, p in enumerate(Pn):
            T = Tn[k]

            dvdtn = [(region2.v(p, i + 0.1) - region2.v(p, i - 0.1)) / 0.2 for i in T]
            dudtn = [(region2.u(p, i + 0.1) - region2.u(p, i - 0.1)) / 0.2 for i in T]
            dhdtn = [(region2.h(p, i + 0.1) - region2.h(p, i - 0.1)) / 0.2 for i in T]
            dsdtn = [(region2.s(p, i + 0.1) - region2.s(p, i - 0.1)) / 0.2 for i in T]
            dgdtn = [(region2.g(p, i + 0.1) - region2.g(p, i - 0.1)) / 0.2 for i in T]

            dvdt  = [region2.dvdT(p, i) for i in T]
            dudt  = [region2.dudT(p, i) for i in T]
            dhdt  = [region2.dhdT(p, i) for i in T]
            dsdt  = [region2.dsdT(p, i) for i in T]
            dgdt  = [region2.dgdT(p, i) for i in T]

            self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state '+str(k+1)+', region 2!')

        ## partial with respect to P
        Tn = [700.0, 554.0]
        Pn = [numpy.linspace(0.0 + 0.5, region3.bnd23P(Tn[0]) - 0.5, n),
              numpy.linspace(0.0 + 0.5, region4.satP(Tn[1]) - 0.5, n)]
        for k, T in enumerate(Tn):
            p = Pn[k]

            dvdpn = [(region2.v(i + 0.001, T) - region2.v(i - 0.001, T)) / (0.002 * 1e3) for i in p]
            dudpn = [(region2.u(i + 0.01, T)  - region2.u(i - 0.01, T))  / (0.02 * 1e3)  for i in p]
            dhdpn = [(region2.h(i + 0.01, T)  - region2.h(i - 0.01, T))  / (0.02 * 1e3)  for i in p]
            dsdpn = [(region2.s(i + 0.001, T) - region2.s(i - 0.001, T)) / (0.002 * 1e3) for i in p]
            dgdpn = [(region2.g(i + 0.001, T) - region2.g(i - 0.001, T)) / (0.002 * 1e3) for i in p]
        
            dvdp  = [region2.dvdP(i, T) for i in p]
            dudp  = [region2.dudP(i, T) for i in p]
            dhdp  = [region2.dhdP(i, T) for i in p]
            dsdp  = [region2.dsdP(i, T) for i in p]
            dgdp  = [region2.dgdP(i, T) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state '+str(k+1)+', region 2!')

    def test_ThermodynamicDerivative_Region3(self):
        n = 100

        ## partial with respect to T
        nun = [0.00144225, 0.01147]
        Tn = [numpy.linspace(623.15 + 0.5, 760.688440, n),
              numpy.linspace(623.15 + 0.5, 539.975, n)]
        for k, nu in enumerate(nun):
            T = Tn[k]

            dpdtn = [(region3.P(nu, i + 0.1) - region3.P(nu, i - 0.1)) / 0.2 for i in T]
            dudtn = [(region3.u(nu, i + 0.1) - region3.u(nu, i - 0.1)) / 0.2 for i in T]
            dhdtn = [(region3.h(nu, i + 0.1) - region3.h(nu, i - 0.1)) / 0.2 for i in T]
            dsdtn = [(region3.s(nu, i + 0.1) - region3.s(nu, i - 0.1)) / 0.2 for i in T]
            dfdtn = [(region3.f(nu, i + 0.1) - region3.f(nu, i - 0.1)) / 0.2 for i in T]
        
            dpdt  = [region3.dPdT(nu, i) for i in T]
            dudt  = [region3.dudT(nu, i) for i in T]
            dhdt  = [region3.dhdT(nu, i) for i in T]
            dsdt  = [region3.dsdT(nu, i) for i in T]
            dfdt  = [region3.dfdT(nu, i) for i in T]

            self.assertLessEqual(abs(sum([(abs(dpdt[i] - dpdtn[i]) / dpdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dpdt, state '+str(k+1)+', region 3!')
            self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state '+str(k+1)+', region 3!')
            self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state '+str(k+1)+', region 3!')
            self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state '+str(k+1)+', region 3!')
            self.assertLessEqual(abs(sum([(abs(dfdt[i] - dfdtn[i]) / dfdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dfdt, state '+str(k+1)+', region 3!')

        ## partial with respect to nu
        Tn = [649.79, 800.0]
        nun = [numpy.linspace(0.00156888 + 5e-8, 0.00787900 - 5e-8, n),
               numpy.linspace(0.00459770 + 5e-8, 0.00220000 - 5e-8, n)]
        for k, T in enumerate(Tn):
            nu = nun[k]

            dpdvn = [(region3.P(i + 1e-8, T) - region3.P(i - 1e-8, T)) / 2e-8 for i in nu]
            dudvn = [(region3.u(i + 1e-8, T) - region3.u(i - 1e-8, T)) / 2e-8 for i in nu]
            dhdvn = [(region3.h(i + 1e-8, T) - region3.h(i - 1e-8, T)) / 2e-8 for i in nu]
            dsdvn = [(region3.s(i + 1e-8, T) - region3.s(i - 1e-8, T)) / 2e-8 for i in nu]
            dfdvn = [(region3.f(i + 1e-8, T) - region3.f(i - 1e-8, T)) / 2e-8 for i in nu]
        
            dpdv  = [region3.dPdv(i, T) for i in nu]
            dudv  = [region3.dudv(i, T) for i in nu]
            dhdv  = [region3.dhdv(i, T) for i in nu]
            dsdv  = [region3.dsdv(i, T) for i in nu]
            dfdv  = [region3.dfdv(i, T) for i in nu]

            self.assertLessEqual(abs(sum([(abs(dpdv[i] - dpdvn[i]) / dpdvn[i]) for i in range(n)]) * 100), 0.05,  'Failed dpdv, state 1, region 3!')
            self.assertLessEqual(abs(sum([(abs(dudv[i] - dudvn[i]) / dudvn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudv, state 1, region 3!')
            self.assertLessEqual(abs(sum([(abs(dhdv[i] - dhdvn[i]) / dhdvn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdv, state 1, region 3!')
            self.assertLessEqual(abs(sum([(abs(dsdv[i] - dsdvn[i]) / dsdvn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdv, state 1, region 3!')
            self.assertLessEqual(abs(sum([(abs(dfdv[i] - dfdvn[i]) / dfdvn[i]) for i in range(n)]) * 100), 0.05,  'Failed dfdv, state 1, region 3!')

    def test_ThermodynamicDerivative_Region4(self):
        n = 10

        ## partial with respect to P @ f
        Pn = [numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n),
              numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n)]
        for k, P in enumerate(Pn):
            dvfdpn = [(h2o.vf(i + i/100) - h2o.vf(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dufdpn = [(h2o.uf(i + i/100) - h2o.uf(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dhfdpn = [(h2o.hf(i + i/100) - h2o.hf(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dsfdpn = [(h2o.sf(i + i/100) - h2o.sf(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dgfdpn = [(h2o.gf(i + i/100) - h2o.gf(i - i/100)) / (2*i/100 * 1e3) for i in P]
        
            dvfdp  = [h2o.dvfdP(i) for i in P]
            dufdp  = [h2o.dufdP(i) for i in P]
            dhfdp  = [h2o.dhfdP(i) for i in P]
            dsfdp  = [h2o.dsfdP(i) for i in P]
            dgfdp  = [h2o.dgfdP(i) for i in P]

            self.assertLessEqual(abs(sum([(abs(dvfdp[i] - dvfdpn[i]) / dvfdpn[i]) for i in range(n)]) * 100), 0.10,  'Failed dvfdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dufdp[i] - dufdpn[i]) / dufdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dufdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dhfdp[i] - dhfdpn[i]) / dhfdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhfdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dsfdp[i] - dsfdpn[i]) / dsfdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsfdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dgfdp[i] - dgfdpn[i]) / dgfdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgfdP, state '+str(k+1)+', region 4!')

        ## partial with respect to P @ g
        Pn = [numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n),
              numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n)]
        for k, P in enumerate(Pn):
            dvgdpn = [(h2o.vg(i + i/100) - h2o.vg(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dugdpn = [(h2o.ug(i + i/100) - h2o.ug(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dhgdpn = [(h2o.hg(i + i/100) - h2o.hg(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dsgdpn = [(h2o.sg(i + i/100) - h2o.sg(i - i/100)) / (2*i/100 * 1e3) for i in P]
            dggdpn = [(h2o.gg(i + i/100) - h2o.gg(i - i/100)) / (2*i/100 * 1e3) for i in P]
        
            dvgdp  = [h2o.dvgdP(i) for i in P]
            dugdp  = [h2o.dugdP(i) for i in P]
            dhgdp  = [h2o.dhgdP(i) for i in P]
            dsgdp  = [h2o.dsgdP(i) for i in P]
            dggdp  = [h2o.dggdP(i) for i in P]

            self.assertLessEqual(abs(sum([(abs(dvgdp[i] - dvgdpn[i]) / dvgdpn[i]) for i in range(n)]) * 100), 0.10,  'Failed dvgdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dugdp[i] - dugdpn[i]) / dugdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dugdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dhgdp[i] - dhgdpn[i]) / dhgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhgdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dsgdp[i] - dsgdpn[i]) / dsgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsgdP, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dggdp[i] - dggdpn[i]) / dggdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dggdP, state '+str(k+1)+', region 4!')

class test_ThermodynamicDerivativeBackwards(unittest.TestCase):

    def test_ThermodynamicDerivative_Region1_Ph(self):
        n = 100

        ## partial with respect to P
        hn = [1500, 50]
        Pn = [numpy.linspace(20 + 0.5, 80 - 0.5, n),
              numpy.linspace(1.0 + 0.5, 50 - 0.5, n)]
        for k, h in enumerate(hn):
            p = Pn[k]
            dvdpn = [(region1.v_h(i + 0.1, h) - region1.v_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dudpn = [(region1.u_h(i + 0.1, h) - region1.u_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dsdpn = [(region1.s_h(i + 0.1, h) - region1.s_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dgdpn = [(region1.g_h(i + 0.1, h) - region1.g_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dTdpn = [(region1.T_h(i + 0.1, h) - region1.T_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
        
            dvdp  = [region1.dvdP_h(i, h) for i in p]
            dudp  = [region1.dudP_h(i, h) for i in p]
            dsdp  = [region1.dsdP_h(i, h) for i in p]
            dgdp  = [region1.dgdP_h(i, h) for i in p]
            dTdp  = [region1.dTdP_h(i, h) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdp_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dudp_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dsdp_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdp_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdp_h, state '+str(k+1)+', region 1!')

        ## partial with respect to h
        Pn = [3, 80]
        hn = [numpy.linspace(50, region4.hf(Pn[0]), n),
              numpy.linspace(50, region4.hf(Pn[1]), n)]
        for k, p in enumerate(Pn):
            h = hn[k]
            dvdhn = [(region1.v_h(p, i + 0.1) - region1.v_h(p, i - 0.1)) / 0.2 for i in h]
            dudhn = [(region1.u_h(p, i + 0.1) - region1.u_h(p, i - 0.1)) / 0.2 for i in h]
            dsdhn = [(region1.s_h(p, i + 0.1) - region1.s_h(p, i - 0.1)) / 0.2 for i in h]
            dgdhn = [(region1.g_h(p, i + 0.1) - region1.g_h(p, i - 0.1)) / 0.2 for i in h]
            dTdhn = [(region1.T_h(p, i + 0.1) - region1.T_h(p, i - 0.1)) / 0.2 for i in h]
        
            dvdh  = [region1.dvdh_h(p, i) for i in h]
            dudh  = [region1.dudh_h(p, i) for i in h]
            dsdh  = [region1.dsdh_h(p, i) for i in h]
            dgdh  = [region1.dgdh_h(p, i) for i in h]
            dTdh  = [region1.dTdh_h(p, i) for i in h]

            self.assertLessEqual(abs(sum([(abs(dvdh[i] - dvdhn[i]) / dvdhn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdh_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dudh[i] - dudhn[i]) / dudhn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dudh_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dsdh[i] - dsdhn[i]) / dsdhn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dsdh_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgdh[i] - dgdhn[i]) / dgdhn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdh_h, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dTdh[i] - dTdhn[i]) / dTdhn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdh_h, state '+str(k+1)+', region 1!')

    def test_ThermodynamicDerivative_Region1_Ps(self):
        n = 100

        ## partial with respect to P
        sn = [1, 3]
        Pn = [numpy.linspace(0.0 + 0.5, 80 - 0.5, n),
              numpy.linspace(20 + 0.5,  80 - 0.5, n)]
        for k, s in enumerate(sn):
            p = Pn[k]
            dvdpn = [(region1.v_s(i + 0.1, s) - region1.v_s(i - 0.1, s)) / (0.2 * 1e3) for i in p]
            dudpn = [(region1.u_s(i + 0.1, s) - region1.u_s(i - 0.1, s)) / (0.2 * 1e3) for i in p]
            dhdpn = [(region1.h_s(i + 0.1, s) - region1.h_s(i - 0.1, s)) / (0.2 * 1e3) for i in p]
            dgdpn = [(region1.g_s(i + 0.1, s) - region1.g_s(i - 0.1, s)) / (0.2 * 1e3) for i in p]
            dTdpn = [(region1.T_s(i + 0.1, s) - region1.T_s(i - 0.1, s)) / (0.2 * 1e3) for i in p]
        
            dvdp  = [region1.dvdP_s(i, s) for i in p]
            dudp  = [region1.dudP_s(i, s) for i in p]
            dhdp  = [region1.dhdP_s(i, s) for i in p]
            dgdp  = [region1.dgdP_s(i, s) for i in p]
            dTdp  = [region1.dTdP_s(i, s) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdp_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 7.5 * n,  'Failed dudp_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dhdp_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdp_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdp_s, state '+str(k+1)+', region 1!')

        ## partial with respect to s
        Pn = [5, 80]
        s = numpy.linspace(0.1, 2.0, n)
        for k, p in enumerate(Pn):
            dvdsn = [(region1.v_s(p, i + 0.1) - region1.v_s(p, i - 0.1)) / 0.2 for i in s]
            dudsn = [(region1.u_s(p, i + 0.1) - region1.u_s(p, i - 0.1)) / 0.2 for i in s]
            dhdsn = [(region1.h_s(p, i + 0.1) - region1.h_s(p, i - 0.1)) / 0.2 for i in s]
            dgdsn = [(region1.g_s(p, i + 0.1) - region1.g_s(p, i - 0.1)) / 0.2 for i in s]
            dTdsn = [(region1.T_s(p, i + 0.1) - region1.T_s(p, i - 0.1)) / 0.2 for i in s]
        
            dvds  = [region1.dvds_s(p, i) for i in s]
            duds  = [region1.duds_s(p, i) for i in s]
            dhds  = [region1.dhds_s(p, i) for i in s]
            dgds  = [region1.dgds_s(p, i) for i in s]
            dTds  = [region1.dTds_s(p, i) for i in s]

            self.assertLessEqual(abs(sum([(abs(dvds[i] - dvdsn[i]) / dvdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvds_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(duds[i] - dudsn[i]) / dudsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed duds_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dhds[i] - dhdsn[i]) / dhdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dhds_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dgds[i] - dgdsn[i]) / dgdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgds_s, state '+str(k+1)+', region 1!')
            self.assertLessEqual(abs(sum([(abs(dTds[i] - dTdsn[i]) / dTdsn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTds_s, state '+str(k+1)+', region 1!')

    def test_ThermodynamicDerivative_Region2_Ph(self):
        n = 100

        ## partial with respect to P
        hn = [2800, 3500]
        p = numpy.linspace(1.0 + 0.5, 80 - 0.5, n)
        for k, h in enumerate(hn):
            dvdpn = [(region2.v_h(i + 0.1, h) - region2.v_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dudpn = [(region2.u_h(i + 0.1, h) - region2.u_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dsdpn = [(region2.s_h(i + 0.1, h) - region2.s_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dgdpn = [(region2.g_h(i + 0.1, h) - region2.g_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
            dTdpn = [(region2.T_h(i + 0.1, h) - region2.T_h(i - 0.1, h)) / (0.2 * 1e3) for i in p]
        
            dvdp  = [region2.dvdP_h(i, h) for i in p]
            dudp  = [region2.dudP_h(i, h) for i in p]
            dsdp  = [region2.dsdP_h(i, h) for i in p]
            dgdp  = [region2.dgdP_h(i, h) for i in p]
            dTdp  = [region2.dTdP_h(i, h) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dvdp_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dudp_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dsdp_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dgdp_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 2.90**2 * n,  'Failed dTdp_h, state '+str(k+1)+', region 2!')

        ## partial with respect to h
        Pn = [0.0035, 10]
        hn = [numpy.linspace(region4.hg(Pn[0]), region2.h(Pn[0], 1073.15), n),
              numpy.linspace(region4.hg(Pn[1]), region2.h(Pn[1], 1073.15), n)]
        for k, p in enumerate(Pn):
            h = hn[k]
            dvdhn = [(region2.v_h(p, i + 0.1) - region2.v_h(p, i - 0.1)) / 0.2 for i in h]
            dudhn = [(region2.u_h(p, i + 0.1) - region2.u_h(p, i - 0.1)) / 0.2 for i in h]
            dsdhn = [(region2.s_h(p, i + 0.1) - region2.s_h(p, i - 0.1)) / 0.2 for i in h]
            dgdhn = [(region2.g_h(p, i + 0.1) - region2.g_h(p, i - 0.1)) / 0.2 for i in h]
            dTdhn = [(region2.T_h(p, i + 0.1) - region2.T_h(p, i - 0.1)) / 0.2 for i in h]
        
            dvdh  = [region2.dvdh_h(p, i) for i in h]
            dudh  = [region2.dudh_h(p, i) for i in h]
            dsdh  = [region2.dsdh_h(p, i) for i in h]
            dgdh  = [region2.dgdh_h(p, i) for i in h]
            dTdh  = [region2.dTdh_h(p, i) for i in h]

            self.assertLessEqual(abs(sum([(abs(dvdh[i] - dvdhn[i]) / dvdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dvdh_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dudh[i] - dudhn[i]) / dudhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dudh_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dsdh[i] - dsdhn[i]) / dsdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dsdh_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgdh[i] - dgdhn[i]) / dgdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dgdh_h, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dTdh[i] - dTdhn[i]) / dTdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dTdh_h, state '+str(k+1)+', region 2!')

    def test_ThermodynamicDerivative_Region2_Ps(self):
        n = 100

        ## partial with respect to P
        sn = [8, 6]
        Pn = [numpy.linspace(0.1, 1.0, n),
              numpy.linspace(10, 80, n)]
        for k, s in enumerate(sn):
            p = Pn[k]
            dvdpn = [(region2.v_s(i + 0.01, s) - region2.v_s(i - 0.01, s)) / (0.02 * 1e3) for i in p]
            dudpn = [(region2.u_s(i + 0.01, s) - region2.u_s(i - 0.01, s)) / (0.02 * 1e3) for i in p]
            dhdpn = [(region2.h_s(i + 0.01, s) - region2.h_s(i - 0.01, s)) / (0.02 * 1e3) for i in p]
            dgdpn = [(region2.g_s(i + 0.01, s) - region2.g_s(i - 0.01, s)) / (0.02 * 1e3) for i in p]
            dTdpn = [(region2.T_s(i + 0.01, s) - region2.T_s(i - 0.01, s)) / (0.02 * 1e3) for i in p]
        
            dvdp  = [region2.dvdP_s(i, s) for i in p]
            dudp  = [region2.dudP_s(i, s) for i in p]
            dhdp  = [region2.dhdP_s(i, s) for i in p]
            dgdp  = [region2.dgdP_s(i, s) for i in p]
            dTdp  = [region2.dTdP_s(i, s) for i in p]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdp_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dudp_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dhdp_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdp_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdp_s, state '+str(k+1)+', region 2!')

        ## partial with respect to s
        Pn = [0.1, 10]
        sn = [numpy.linspace(8, 9, n),
              numpy.linspace(6, 7, n)]
        for k, p in enumerate(Pn):
            s = sn[k]
            dvdsn = [(region2.v_s(p, i + 0.01) - region2.v_s(p, i - 0.01)) / 0.02 for i in s]
            dudsn = [(region2.u_s(p, i + 0.01) - region2.u_s(p, i - 0.01)) / 0.02 for i in s]
            dhdsn = [(region2.h_s(p, i + 0.01) - region2.h_s(p, i - 0.01)) / 0.02 for i in s]
            dgdsn = [(region2.g_s(p, i + 0.01) - region2.g_s(p, i - 0.01)) / 0.02 for i in s]
            dTdsn = [(region2.T_s(p, i + 0.01) - region2.T_s(p, i - 0.01)) / 0.02 for i in s]
        
            dvds  = [region2.dvds_s(p, i) for i in s]
            duds  = [region2.duds_s(p, i) for i in s]
            dhds  = [region2.dhds_s(p, i) for i in s]
            dgds  = [region2.dgds_s(p, i) for i in s]
            dTds  = [region2.dTds_s(p, i) for i in s]

            self.assertLessEqual(abs(sum([(abs(dvds[i] - dvdsn[i]) / dvdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvds_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(duds[i] - dudsn[i]) / dudsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed duds_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dhds[i] - dhdsn[i]) / dhdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dhds_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dgds[i] - dgdsn[i]) / dgdsn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgds_s, state '+str(k+1)+', region 2!')
            self.assertLessEqual(abs(sum([(abs(dTds[i] - dTdsn[i]) / dTdsn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTds_s, state '+str(k+1)+', region 2!')

    def test_ThermodynamicDerivative_Region4_Ph(self):
        n = 10

        ## partial with respect to p
        hn = [region4.h_h(h2o.satP(623.15), 0.4),
              region4.h_h(h2o.satP(623.15), 0.8)]
        P = numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n)
        for k, h in enumerate(hn):
            dvdpn = [(region4.v_h(i + i/100, h) - region4.v_h(i - i/100, h)) / (2*i/100 * 1e3) for i in P]
            dudpn = [(region4.u_h(i + i/100, h) - region4.u_h(i - i/100, h)) / (2*i/100 * 1e3) for i in P]
            dsdpn = [(region4.s_h(i + i/100, h) - region4.s_h(i - i/100, h)) / (2*i/100 * 1e3) for i in P]
            dgdpn = [(region4.g_h(i + i/100, h) - region4.g_h(i - i/100, h)) / (2*i/100 * 1e3) for i in P]
        
            dvdp  = [region4.dvdP_h(i, h) for i in P]
            dudp  = [region4.dudP_h(i, h) for i in P]
            dsdp  = [region4.dsdP_h(i, h) for i in P]
            dgdp  = [region4.dgdP_h(i, h) for i in P]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.10,  'Failed dvdP_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudP_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdP_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdP_h, state '+str(k+1)+', region 4!')

        ## partial with respect to h
        Pn = [1, 10]
        hn = [numpy.linspace(region4.hf(Pn[0]), region4.hg(Pn[0]), n),
              numpy.linspace(region4.hf(Pn[1]), region4.hg(Pn[1]), n)]
        for k, p in enumerate(Pn):
            h = hn[k]
            dvdhn = [(region4.v_h(p, i + 0.1) - region4.v_h(p, i - 0.1)) / 0.2 for i in h]
            dudhn = [(region4.u_h(p, i + 0.1) - region4.u_h(p, i - 0.1)) / 0.2 for i in h]
            dsdhn = [(region4.s_h(p, i + 0.1) - region4.s_h(p, i - 0.1)) / 0.2 for i in h]
            dgdhn = [(region4.g_h(p, i + 0.1) - region4.g_h(p, i - 0.1)) / 0.2 for i in h]
            dxdhn = [(region4.x_h(p, i + 0.1) - region4.x_h(p, i - 0.1)) / 0.2 for i in h]
        
            dvdh  = [region4.dvdh_h(p, i) for i in h]
            dudh  = [region4.dudh_h(p, i) for i in h]
            dsdh  = [region4.dsdh_h(p, i) for i in h]
            dgdh  = [region4.dgdh_h(p, i) for i in h]
            dxdh  = [region4.dxdh_h(p) for i in h]

            self.assertLessEqual(abs(sum([(abs(dvdh[i] - dvdhn[i]) / dvdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dvdh_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dudh[i] - dudhn[i]) / dudhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dudh_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dsdh[i] - dsdhn[i]) / dsdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dsdh_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dgdh[i] - dgdhn[i]) / dgdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dgdh_h, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dxdh[i] - dxdhn[i]) / dxdhn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dTdh_h, state '+str(k+1)+', region 4!')

    def test_ThermodynamicDerivative_Region4_Ps(self):
        n = 10

        ## partial with respect to p
        sn = [region4.s_s(h2o.satP(623.15), 0.4),
              region4.s_s(h2o.satP(623.15), 0.8)]
        P = numpy.logspace(-3, numpy.log10(h2o.satP(623.15) -  0.1), n)
        for k, s in enumerate(sn):
            dvdpn = [(region4.v_s(i + i/100, s) - region4.v_s(i - i/100, s)) / (2*i/100 * 1e3) for i in P]
            dudpn = [(region4.u_s(i + i/100, s) - region4.u_s(i - i/100, s)) / (2*i/100 * 1e3) for i in P]
            dhdpn = [(region4.h_s(i + i/100, s) - region4.h_s(i - i/100, s)) / (2*i/100 * 1e3) for i in P]
            dgdpn = [(region4.g_s(i + i/100, s) - region4.g_s(i - i/100, s)) / (2*i/100 * 1e3) for i in P]
        
            dvdp  = [region4.dvdP_s(i, s) for i in P]
            dudp  = [region4.dudP_s(i, s) for i in P]
            dhdp  = [region4.dhdP_s(i, s) for i in P]
            dgdp  = [region4.dgdP_s(i, s) for i in P]

            self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.10,  'Failed dvdP_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudP_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdP_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdP_s, state '+str(k+1)+', region 4!')

        ## partial with respect to h
        Pn = [1, 10]
        sn = [numpy.linspace(region4.sf(Pn[0]), region4.sg(Pn[0]), n),
              numpy.linspace(region4.sf(Pn[1]), region4.sg(Pn[1]), n)]
        for k, p in enumerate(Pn):
            s = sn[k]
            dvdsn = [(region4.v_s(p, i + 0.1) - region4.v_s(p, i - 0.1)) / 0.2 for i in s]
            dudsn = [(region4.u_s(p, i + 0.1) - region4.u_s(p, i - 0.1)) / 0.2 for i in s]
            dhdsn = [(region4.h_s(p, i + 0.1) - region4.h_s(p, i - 0.1)) / 0.2 for i in s]
            dgdsn = [(region4.g_s(p, i + 0.1) - region4.g_s(p, i - 0.1)) / 0.2 for i in s]
            dxdsn = [(region4.x_s(p, i + 0.1) - region4.x_s(p, i - 0.1)) / 0.2 for i in s]
        
            dvds  = [region4.dvds_s(p, i) for i in s]
            duds  = [region4.duds_s(p, i) for i in s]
            dhds  = [region4.dhds_s(p, i) for i in s]
            dgds  = [region4.dgds_s(p, i) for i in s]
            dxds  = [region4.dxds_s(p) for i in s]

            self.assertLessEqual(abs(sum([(abs(dvds[i] - dvdsn[i]) / dvdsn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dvds_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(duds[i] - dudsn[i]) / dudsn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed duds_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dhds[i] - dhdsn[i]) / dhdsn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dhds_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dgds[i] - dgdsn[i]) / dgdsn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dgds_s, state '+str(k+1)+', region 4!')
            self.assertLessEqual(abs(sum([(abs(dxds[i] - dxdsn[i]) / dxdsn[i]) for i in range(n)]) * 100), 0.05 * n,  'Failed dTds_s, state '+str(k+1)+', region 4!')

if __name__ == '__main__':
    unittest.main()
