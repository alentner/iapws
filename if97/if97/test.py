import unittest
from if97 import region1, region2, region3, region4, h2o
import numpy

class test_ThermodynamicProperty(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_State1(self):
        self.assertEqual(round(region1.v(3, 300), 11), 0.100215168e-2, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.u(3, 300), 6),  0.112324818e3,  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(3, 300), 9),  0.392294792,    'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.h(3, 300), 6),  0.115331273e3,  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.cp(3, 300), 8), 0.417301218e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(3, 300), 8), 0.412120160e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(3, 300), 5),  0.150773921e4,  'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region1_State2(self):
        self.assertEqual(round(region1.v(80, 300), 12), 0.971180894e-3, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.u(80, 300), 6),  0.106448356e3,  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(80, 300), 9),  0.368563852,    'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.h(80, 300), 6),  0.184142828e3,  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.cp(80, 300), 8), 0.401008987e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(80, 300), 8), 0.391736606e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(80, 300), 5),  0.163469054e4,  'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region1_State3(self):
        self.assertEqual(round(region1.v(3, 500), 11), 0.120241800e-2, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.u(3, 500), 6),  0.971934985e3,  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(3, 500), 8),  0.258041912e1,  'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.h(3, 500), 6),  0.975542239e3,  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.cp(3, 500), 8), 0.465580682e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(3, 500), 8), 0.322139223e1,  'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(3, 500), 5),  0.124071337e4,  'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region2_State1(self):
        self.assertEqual(round(region2.v(0.0035, 300), 7),  0.394913866e2, 'Failed specific volume, state 1, region 2!')
        self.assertEqual(round(region2.u(0.0035, 300), 5),  0.241169160e4, 'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(round(region2.s(0.0035, 300), 8),  0.852238967e1, 'Failed specific entropy, state 1, region 2!')
        self.assertEqual(round(region2.h(0.0035, 300), 5),  0.254991145e4, 'Failed specific enthalpy, state 1, region 2!')
        self.assertEqual(round(region2.cp(0.0035, 300), 8), 0.191300162e1, 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(round(region2.w(0.0035, 300), 6),  0.427920172e3, 'Failed speed of sound, state 1, region 2!')
    def test_ThermodynamicProperty_Region2_State2(self):
        self.assertEqual(round(region2.v(0.0035, 700), 7),  0.923015898e2, 'Failed specific volume, state 2, region 2!')
        self.assertEqual(round(region2.u(0.0035, 700), 5),  0.301262819e4, 'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(round(region2.s(0.0035, 700), 7),  0.101749996e2, 'Failed specific entropy, state 2, region 2!')
        self.assertEqual(round(region2.h(0.0035, 700), 5),  0.333568375e4, 'Failed specific enthalpy, state 2, region 2!')
        self.assertEqual(round(region2.cp(0.0035, 700), 8), 0.208141274e1, 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(round(region2.w(0.0035, 700), 6),  0.644289068e3, 'Failed speed of sound, state 2, region 2!')
    def test_ThermodynamicProperty_Region2_State3(self):
        self.assertEqual(round(region2.v(30, 700), 11), 0.542946619e-2, 'Failed specific volume, state 3, region 2!')
        self.assertEqual(round(region2.u(30, 700), 5),  0.246861076e4,  'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(round(region2.s(30, 700), 8),  0.517540298e1,  'Failed specific entropy, state 3, region 2!')
        self.assertEqual(round(region2.h(30, 700), 5),  0.263149474e4,  'Failed specific enthalpy, state 3, region 2!')
        self.assertEqual(round(region2.cp(30, 700), 7), 0.103505092e2,  'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(round(region2.w(30, 700), 6),  0.480386523e3,  'Failed speed of sound, state 3, region 2!')
    def test_ThermodynamicProperty_Region3_bnd23(self):
        self.assertEqual(round(region3.bnd23T(0.165291643e2), 6), 0.623150000e3, 'Failed boundary equation, region 23 by T!')
        self.assertEqual(round(region3.bnd23P(0.623150000e3), 7), 0.165291643e2, 'Failed boundary equation, region 23 by P!') 
    def test_ThermodynamicProperty_Region3_State1(self):
        self.assertEqual(round(region3.P(1 / 500, 650), 7),  0.255837018e2, 'Failed pressure, state 1, region 3!')
        self.assertEqual(round(region3.u(1 / 500, 650), 5),  0.181226279e4, 'Failed specific internal energy, state 1, region 3!')
        self.assertEqual(round(region3.s(1 / 500, 650), 8),  0.405427273e1, 'Failed specific entropy, state 1, region 3!')
        self.assertEqual(round(region3.h(1 / 500, 650), 5),  0.186343019e4, 'Failed specific enthalpy, state 1, region 3!')
        self.assertEqual(round(region3.cp(1 / 500, 650), 7), 0.138935717e2, 'Failed specific heat capacity, state 1, region 3!')
        self.assertEqual(round(region3.w(1 / 500, 650), 6),  0.502005554e3, 'Failed speed of sound, state 1, region 3!')
    def test_ThermodynamicProperty_Region3_State2(self):
        self.assertEqual(round(region3.P(1 / 200, 650), 7),  0.222930643e2, 'Failed pressure, state 1, region 3!')
        self.assertEqual(round(region3.u(1 / 200, 650), 5),  0.226365868e4, 'Failed specific internal energy, state 1, region 3!')
        self.assertEqual(round(region3.s(1 / 200, 650), 8),  0.485438792e1, 'Failed specific entropy, state 1, region 3!')
        self.assertEqual(round(region3.h(1 / 200, 650), 5),  0.237512401e4, 'Failed specific enthalpy, state 1, region 3!')
        self.assertEqual(round(region3.cp(1 / 200, 650), 7), 0.446579342e2, 'Failed specific heat capacity, state 1, region 3!')
        self.assertEqual(round(region3.w(1 / 200, 650), 6),  0.383444594e3, 'Failed speed of sound, state 1, region 3!')
    def test_ThermodynamicProperty_Region3_State3(self):
        self.assertEqual(round(region3.P(1 / 500, 750), 7), 0.783095639e2, 'Failed pressure, state 1, region 3!')
        self.assertEqual(round(region3.u(1 / 500, 750), 5),  0.210206932e4, 'Failed specific internal energy, state 1, region 3!')
        self.assertEqual(round(region3.s(1 / 500, 750), 8),  0.446971906e1, 'Failed specific entropy, state 1, region 3!')
        self.assertEqual(round(region3.h(1 / 500, 750), 5),  0.225868845e4, 'Failed specific enthalpy, state 1, region 3!')
        self.assertEqual(round(region3.cp(1 / 500, 750), 8), 0.634165359e1, 'Failed specific heat capacity, state 1, region 3!')
        self.assertEqual(round(region3.w(1 / 500, 750), 6),  0.760696041e3, 'Failed speed of sound, state 1, region 3!')
    def test_ThermodynamicProperty_Region4_satP(self):
        self.assertEqual(round(region4.satP(300), 11), 0.353658941e-2, 'Failed satuation pressure, 300K!') 
        self.assertEqual(round(region4.satP(500), 8),  0.263889776e1, 'Failed satuation pressure, 500K!') 
        self.assertEqual(round(region4.satP(600), 7),  0.123443146e2, 'Failed satuation pressure, 600K!') 
    def test_ThermodynamicProperty_Region4_satT(self):
        self.assertEqual(round(region4.satT(0.10), 6), 0.372755919e3, 'Failed satuation pressure, 0.1 MPa!') 
        self.assertEqual(round(region4.satT(1.00), 6), 0.453035632e3, 'Failed satuation pressure, 1.0 MPa!') 
        self.assertEqual(round(region4.satT(10.0), 6), 0.584149488e3, 'Failed satuation pressure, 10 MPa!') 
class test_ThermodynamicDerivative(unittest.TestCase):
    def test_ThermodynamicPartialT_Region1_state1(self):
        n = 100
        p = 3
        T = numpy.linspace(273.15 + 0.5, region4.satT(p) - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state 1, region 1!')
    def test_ThermodynamicPartialT_Region1_state2(self):
        n = 100
        p = 80
        T = numpy.linspace(273.15 + 0.5, 623.15 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state 2, region 1!')
    def test_ThermodynamicPartialP_Region1_state1(self):
        n = 100
        T = 300
        p = numpy.linspace(region4.satP(T) + 0.5, 100 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state 1, region 1!')
    def test_ThermodynamicPartialP_Region1_state2(self):
        n = 100
        T = 500
        p = numpy.linspace(region4.satP(T) + 0.5, 100 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state 2, region 1!')
    def test_ThermodynamicPartialT_Region2_state1(self):
        n = 100
        p = 0.0035
        T = numpy.linspace(region4.satT(p) + 0.5, 1073.15 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state 1, region 2!')
    def test_ThermodynamicPartialT_Region2_state2(self):
        n = 100
        p = 10
        T = numpy.linspace(region4.satT(p) + 0.5, 1073.15 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdt[i] - dvdtn[i]) / dvdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdt, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdt[i] - dgdtn[i]) / dgdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdt, state 2, region 2!')
    def test_ThermodynamicPartialP_Region2_state1(self):
        n = 100
        T = 700.0
        p = numpy.linspace(0.0 + 0.5, region3.bnd23P(T) - 0.5, n)

        dvdpn = [(region2.v(i + 0.001, T) - region2.v(i - 0.001, T)) / (0.002 * 1e3) for i in p]
        dudpn = [(region2.u(i + 0.01, T)  - region2.u(i - 0.01, T))  / (0.02 * 1e3)  for i in p]
        dhdpn = [(region2.h(i + 0.1, T)   - region2.h(i - 0.1, T))   / (0.2 * 1e3)   for i in p]
        dsdpn = [(region2.s(i + 0.01, T)  - region2.s(i - 0.01, T))  / (0.02 * 1e3)  for i in p]
        dgdpn = [(region2.g(i + 0.01, T)  - region2.g(i - 0.01, T))  / (0.02 * 1e3)  for i in p]
        
        dvdp  = [region2.dvdP(i, T) for i in p]
        dudp  = [region2.dudP(i, T) for i in p]
        dhdp  = [region2.dhdP(i, T) for i in p]
        dsdp  = [region2.dsdP(i, T) for i in p]
        dgdp  = [region2.dgdP(i, T) for i in p]

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state 1, region 2!')
    def test_ThermodynamicPartialP_Region2_state2(self):
        n = 50
        T = 554
        p = numpy.linspace(0.0 + 0.5, region4.satP(T) - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dvdp, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudp, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dhdp[i] - dhdpn[i]) / dhdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdp, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdp, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 0.05,  'Failed dgdp, state 2, region 2!')
    def test_ThermodynamicPartialT_Region3_state1(self):
        n = 100
        nu = 0.00144225
        T = numpy.linspace(623.15 + 0.5, 760.688440, n)

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

        self.assertLessEqual(abs(sum([(abs(dpdt[i] - dpdtn[i]) / dpdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dpdt, state 1, region 3!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 1, region 3!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 1, region 3!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 1, region 3!')
        self.assertLessEqual(abs(sum([(abs(dfdt[i] - dfdtn[i]) / dfdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dfdt, state 1, region 3!')
    def test_ThermodynamicPartialT_Region3_state2(self):
        n = 100
        nu = 0.01147
        T = numpy.linspace(623.15 + 0.5, 539.975, n)

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

        self.assertLessEqual(abs(sum([(abs(dpdt[i] - dpdtn[i]) / dpdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dpdt, state 2, region 3!')
        self.assertLessEqual(abs(sum([(abs(dudt[i] - dudtn[i]) / dudtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dudt, state 2, region 3!')
        self.assertLessEqual(abs(sum([(abs(dhdt[i] - dhdtn[i]) / dhdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dhdt, state 2, region 3!')
        self.assertLessEqual(abs(sum([(abs(dsdt[i] - dsdtn[i]) / dsdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dsdt, state 2, region 3!')
        self.assertLessEqual(abs(sum([(abs(dfdt[i] - dfdtn[i]) / dfdtn[i]) for i in range(n)]) * 100), 0.05,  'Failed dfdt, state 2, region 3!')
    def test_ThermodynamicPartialv_Region3_state1(self):
        n = 100
        T = 649.79
        nu = numpy.linspace(0.00156888 + 5e-8, 0.00787900 - 5e-8, n)

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
    def test_ThermodynamicPartialv_Region3_state2(self):
        n = 100
        T = 800.0
        nu = numpy.linspace(0.00459770 + 5e-8, 0.00220000 - 5e-8, n)

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
class test_ThermodynamicPropertyBackwards(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_backwards(self):
        self.assertEqual(round(region1.T_h(3, 500), 6),   0.391798509e3, 'Failed backward temperature, state 1, region 1!') 
        self.assertEqual(round(region1.T_h(80, 500), 6),  0.378108626e3, 'Failed backward temperature, state 2, region 1!')
        self.assertEqual(round(region1.T_h(80, 1500), 6), 0.611041229e3, 'Failed backward temperature, state 3, region 1!')

        self.assertAlmostEqual(region1.v(3, 300)  / region1.v_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed v consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.u(3, 300)  / region1.u_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed u consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.s(3, 300)  / region1.s_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed s consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.g(3, 300)  / region1.g_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed g consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.cp(3, 300) / region1.cp_h(3, region1.h(3, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.cv(3, 300) / region1.cv_h(3, region1.h(3, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.w(3, 300)  / region1.w_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed w consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.a(3, 300)  / region1.a_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed a consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.k(3, 300)  / region1.k_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed k consistancy, state 1, region 1!')

        self.assertAlmostEqual(region1.v(80, 300)  / region1.v_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed v consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.u(80, 300)  / region1.u_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed u consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.s(80, 300)  / region1.s_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed s consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.g(80, 300)  / region1.g_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed g consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.cp(80, 300) / region1.cp_h(80, region1.h(80, 300)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.cv(80, 300) / region1.cv_h(80, region1.h(80, 300)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.w(80, 300)  / region1.w_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed w consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.a(80, 300)  / region1.a_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed a consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.k(80, 300)  / region1.k_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed k consistancy, state 2, region 1!')
        
        self.assertAlmostEqual(region1.v(3, 500)  / region1.v_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed v consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.u(3, 500)  / region1.u_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed u consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.s(3, 500)  / region1.s_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed s consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.g(3, 500)  / region1.g_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed g consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.cp(3, 500) / region1.cp_h(3, region1.h(3, 500)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.cv(3, 500) / region1.cv_h(3, region1.h(3, 500)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.w(3, 500)  / region1.w_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed w consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.a(3, 500)  / region1.a_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed a consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.k(3, 500)  / region1.k_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed k consistancy, state 3, region 1!')
    def test_ThermodynamicProperty_Region2_backwards(self):
        self.assertEqual(round(region2.bnd2b2c(0.100e3), 6), 0.3516004323e4, 'Failed boundary equation, region 2b-2c!')  
        self.assertEqual(round(region2.T_h(0.001, 3000), 6), 0.534433241e3, 'Failed backward temperature, state 1, region 2a!')
        self.assertEqual(round(region2.T_h(3, 3000), 6),     0.575373370e3, 'Failed backward temperature, state 2, region 2a!') 
        self.assertEqual(round(region2.T_h(3, 4000), 5),     0.101077577e4, 'Failed backward temperature, state 3, region 2a!')         
        self.assertEqual(round(region2.T_h(5, 3500), 6),     0.801299102e3, 'Failed backward temperature, state 1, region 2b!')
        self.assertEqual(round(region2.T_h(5, 4000), 5),     0.101531583e4, 'Failed backward temperature, state 2, region 2b!') 
        self.assertEqual(round(region2.T_h(25, 3500), 6),    0.875279054e3, 'Failed backward temperature, state 3, region 2b!')
        self.assertEqual(round(region2.T_h(40, 2700), 6),    0.743056411e3, 'Failed backward temperature, state 1, region 2c!')
        self.assertEqual(round(region2.T_h(60, 2700), 6),    0.791137067e3, 'Failed backward temperature, state 2, region 2c!') 
        self.assertEqual(round(region2.T_h(60, 3200), 6),    0.882756860e3, 'Failed backward temperature, state 3, region 2c!')

        self.assertAlmostEqual(region2.v(0.0035, 300)  / region2.v_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed v consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 300)  / region2.u_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed u consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.s(0.0035, 300)  / region2.s_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed s consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 300)  / region2.g_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=1, msg='Failed g consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 300) / region2.cp_h(0.0035, region2.h(0.0035, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 300) / region2.cv_h(0.0035, region2.h(0.0035, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 300)  / region2.w_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed w consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 300)  / region2.a_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed a consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 300)  / region2.k_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed k consistancy, state 1, region 2!')

        self.assertAlmostEqual(region2.v(0.0035, 700)  / region2.v_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed v consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 700)  / region2.u_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed u consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.s(0.0035, 700)  / region2.s_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed s consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 700)  / region2.g_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed g consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 700) / region2.cp_h(0.0035, region2.h(0.0035, 700)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 700) / region2.cv_h(0.0035, region2.h(0.0035, 700)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 700)  / region2.w_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed w consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 700)  / region2.a_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed a consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 700)  / region2.k_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed k consistancy, state 2, region 2!')
        
        self.assertAlmostEqual(region2.v(30, 700)  / region2.v_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed v consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.u(30, 700)  / region2.u_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed u consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.s(30, 700)  / region2.s_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed s consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.g(30, 700)  / region2.g_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed g consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.cp(30, 700) / region2.cp_h(30, region2.h(30, 700)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.cv(30, 700) / region2.cv_h(30, region2.h(30, 700)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.w(30, 700)  / region2.w_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed w consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.a(30, 700)  / region2.a_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed a consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.k(30, 700)  / region2.k_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed k consistancy, state 3, region 2!')
class test_ThermodynamicDerivativeBackwards(unittest.TestCase):
    def test_ThermodynamicPartialP_Region1_state1_backwards(self):
        n = 100
        h = 1500
        p = numpy.linspace(20 + 0.5, 80 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdp_h, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dudp_h, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dsdp_h, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdp_h, state 1, region 1!')
        self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdp_h, state 1, region 1!')
    def test_ThermodynamicPartialP_Region1_state2_backwards(self):
        n = 100
        h = 50
        p = numpy.linspace(1.0 + 0.5, 50 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dvdp_h, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dudp_h, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dsdp_h, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 1.34 / 2 * n,  'Failed dgdp_h, state 2, region 1!')
        self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 1.34**2  * n,  'Failed dTdp_h, state 2, region 1!')
    def test_ThermodynamicPartialP_Region2_state1_backwards(self):
        n = 100
        h = 2800
        p = numpy.linspace(1.0 + 0.5, 80 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dvdp_h, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dudp_h, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dsdp_h, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dgdp_h, state 1, region 2!')
        self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 2.90**2  * n,  'Failed dTdp_h, state 1, region 2!')
    def test_ThermodynamicPartialP_Region2_state2_backwards(self):
        n = 100
        h = 3500
        p = numpy.linspace(1.0 + 0.5, 80 - 0.5, n)

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

        self.assertLessEqual(abs(sum([(abs(dvdp[i] - dvdpn[i]) / dvdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dvdp_h, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dudp[i] - dudpn[i]) / dudpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dudp_h, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dsdp[i] - dsdpn[i]) / dsdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dsdp_h, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dgdp[i] - dgdpn[i]) / dgdpn[i]) for i in range(n)]) * 100), 2.90 * n,  'Failed dgdp_h, state 2, region 2!')
        self.assertLessEqual(abs(sum([(abs(dTdp[i] - dTdpn[i]) / dTdpn[i]) for i in range(n)]) * 100), 2.90**2  * n,  'Failed dTdp_h, state 2, region 2!')
class test_ThermodynamicWrapper(unittest.TestCase):
    def test_ThermodynamicWrapper_Region1_State1(self):
        self.assertEqual(h2o.v(3, 300),  region1.v(3, 300),  'Failed specific volume, state 1, region 1!')
        self.assertEqual(h2o.u(3, 300),  region1.u(3, 300),  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(h2o.s(3, 300),  region1.s(3, 300),  'Failed specific entropy, state 1, region 1!')
        self.assertEqual(h2o.h(3, 300),  region1.h(3, 300),  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(h2o.cp(3, 300), region1.cp(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.cv(3, 300), region1.cv(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.w(3, 300),  region1.w(3, 300),  'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicWrapper_Region1_State2(self):
        self.assertEqual(h2o.v(80, 300),  region1.v(80, 300),  'Failed specific volume, state 2, region 1!')
        self.assertEqual(h2o.u(80, 300),  region1.u(80, 300),  'Failed specific internal energy, state 2, region 1!')
        self.assertEqual(h2o.s(80, 300),  region1.s(80, 300),  'Failed specific entropy, state 2, region 1!')
        self.assertEqual(h2o.h(80, 300),  region1.h(80, 300),  'Failed specific enthalpy, state 2, region 1!')
        self.assertEqual(h2o.cp(80, 300), region1.cp(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.cv(80, 300), region1.cv(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.w(80, 300),  region1.w(80, 300),  'Failed speed of sound, state 1, region 2!')
    def test_ThermodynamicWrapper_Region1_State3(self):
        self.assertEqual(h2o.v(3, 500),  region1.v(3, 500),  'Failed specific volume, state 3, region 1!')
        self.assertEqual(h2o.u(3, 500),  region1.u(3, 500),  'Failed specific internal energy, state 3, region 1!')
        self.assertEqual(h2o.s(3, 500),  region1.s(3, 500),  'Failed specific entropy, state 3, region 1!')
        self.assertEqual(h2o.h(3, 500),  region1.h(3, 500),  'Failed specific enthalpy, state 3, region 1!')
        self.assertEqual(h2o.cp(3, 500), region1.cp(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.cv(3, 500), region1.cv(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.w(3, 500),  region1.w(3, 500),  'Failed speed of sound, state 1, region 3!')
    def test_ThermodynamicWrapper_Region2_State1(self):
        self.assertEqual(h2o.v(0.0035, 300),  region2.v(0.0035, 300),  'Failed specific volume, state 1, region 2!')
        self.assertEqual(h2o.u(0.0035, 300),  region2.u(0.0035, 300),  'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(h2o.s(0.0035, 300),  region2.s(0.0035, 300),  'Failed specific entropy, state 1, region 2!')
        self.assertEqual(h2o.h(0.0035, 300),  region2.h(0.0035, 300),  'Failed specific enthalpy, state 1, region 2!')
        self.assertEqual(h2o.cp(0.0035, 300), region2.cp(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.cv(0.0035, 300), region2.cv(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.w(0.0035, 300),  region2.w(0.0035, 300),  'Failed speed of sound, state 1, region 2!')
    def test_ThermodynamicWrapper_Region2_State2(self):
        self.assertEqual(h2o.v(0.0035, 700),  region2.v(0.0035, 700),  'Failed specific volume, state 2, region 2!')
        self.assertEqual(h2o.u(0.0035, 700),  region2.u(0.0035, 700),  'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(h2o.s(0.0035, 700),  region2.s(0.0035, 700),  'Failed specific entropy, state 2, region 2!')
        self.assertEqual(h2o.h(0.0035, 700),  region2.h(0.0035, 700),  'Failed specific enthalpy, state 2, region 2!')
        self.assertEqual(h2o.cp(0.0035, 700), region2.cp(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.cv(0.0035, 700), region2.cv(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.w(0.0035, 700),  region2.w(0.0035, 700),  'Failed speed of sound, state 2, region 2!')
    def test_ThermodynamicWrapper_Region2_State3(self):
        self.assertEqual(h2o.v(30, 700),  region2.v(30, 700),  'Failed specific volume, state 3, region 2!')
        self.assertEqual(h2o.u(30, 700),  region2.u(30, 700),  'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(h2o.s(30, 700),  region2.s(30, 700),  'Failed specific entropy, state 3, region 2!')
        self.assertEqual(h2o.h(30, 700),  region2.h(30, 700),  'Failed specific enthalpy, state 3, region 2!')
        self.assertEqual(h2o.cp(30, 700), region2.cp(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.cv(30, 700), region2.cv(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.w(30, 700),  region2.w(30, 700),  'Failed speed of sound, state 3, region 2!')

if __name__ == '__main__':
    unittest.main()