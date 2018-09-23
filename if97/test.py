import unittest
from if97 import region1, region2, region3, region4, h2o, moody
import numpy
import matplotlib
from matplotlib import pyplot, cm, patches, lines

class test_ThermodynamicProperty(unittest.TestCase):
    def test_ThermodynamicProperty_Region1(self):
        self.assertEqual(round(region1.v(3, 300), 11), 0.100215168e-2,  'Failed specific volume,          state 1, region 1!')
        self.assertEqual(round(region1.u(3, 300), 6),  0.112324818e3,   'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(3, 300), 9),  0.392294792,     'Failed specific entropy,         state 1, region 1!')
        self.assertEqual(round(region1.h(3, 300), 6),  0.115331273e3,   'Failed specific enthalpy,        state 1, region 1!')
        self.assertEqual(round(region1.cp(3, 300), 8), 0.417301218e1,   'Failed specific heat capacity,   state 1, region 1!')
        self.assertEqual(round(region1.cv(3, 300), 8), 0.412120160e1,   'Failed specific heat capacity,   state 1, region 1!')
        self.assertEqual(round(region1.w(3, 300), 5),  0.150773921e4,   'Failed speed of sound,           state 1, region 1!')

        self.assertEqual(round(region1.v(80, 300), 12), 0.971180894e-3, 'Failed specific volume,          state 2, region 1!')
        self.assertEqual(round(region1.u(80, 300), 6),  0.106448356e3,  'Failed specific internal energy, state 2, region 1!')
        self.assertEqual(round(region1.s(80, 300), 9),  0.368563852,    'Failed specific entropy,         state 2, region 1!')
        self.assertEqual(round(region1.h(80, 300), 6),  0.184142828e3,  'Failed specific enthalpy,        state 2, region 1!')
        self.assertEqual(round(region1.cp(80, 300), 8), 0.401008987e1,  'Failed specific heat capacity,   state 2, region 1!')
        self.assertEqual(round(region1.cv(80, 300), 8), 0.391736606e1,  'Failed specific heat capacity,   state 2, region 1!')
        self.assertEqual(round(region1.w(80, 300), 5),  0.163469054e4,  'Failed speed of sound,           state 2, region 1!')

        self.assertEqual(round(region1.v(3, 500), 11), 0.120241800e-2,  'Failed specific volume,          state 3, region 1!')
        self.assertEqual(round(region1.u(3, 500), 6),  0.971934985e3,   'Failed specific internal energy, state 3, region 1!')
        self.assertEqual(round(region1.s(3, 500), 8),  0.258041912e1,   'Failed specific entropy,         state 3, region 1!')
        self.assertEqual(round(region1.h(3, 500), 6),  0.975542239e3,   'Failed specific enthalpy,        state 3, region 1!')
        self.assertEqual(round(region1.cp(3, 500), 8), 0.465580682e1,   'Failed specific heat capacity,   state 3, region 1!')
        self.assertEqual(round(region1.cv(3, 500), 8), 0.322139223e1,   'Failed specific heat capacity,   state 3, region 1!')
        self.assertEqual(round(region1.w(3, 500), 5),  0.124071337e4,   'Failed speed of sound,           state 3, region 1!')
    def test_ThermodynamicProperty_Region2(self):
        self.assertEqual(round(region2.v(0.0035, 300), 7),  0.394913866e2, 'Failed specific volume,          state 1, region 2!')
        self.assertEqual(round(region2.u(0.0035, 300), 5),  0.241169160e4, 'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(round(region2.s(0.0035, 300), 8),  0.852238967e1, 'Failed specific entropy,         state 1, region 2!')
        self.assertEqual(round(region2.h(0.0035, 300), 5),  0.254991145e4, 'Failed specific enthalpy,        state 1, region 2!')
        self.assertEqual(round(region2.cp(0.0035, 300), 8), 0.191300162e1, 'Failed specific heat capacity,   state 1, region 2!')
        self.assertEqual(round(region2.w(0.0035, 300), 6),  0.427920172e3, 'Failed speed of sound,           state 1, region 2!')
        
        self.assertEqual(round(region2.v(0.0035, 700), 7),  0.923015898e2, 'Failed specific volume,          state 2, region 2!')
        self.assertEqual(round(region2.u(0.0035, 700), 5),  0.301262819e4, 'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(round(region2.s(0.0035, 700), 7),  0.101749996e2, 'Failed specific entropy,         state 2, region 2!')
        self.assertEqual(round(region2.h(0.0035, 700), 5),  0.333568375e4, 'Failed specific enthalpy,        state 2, region 2!')
        self.assertEqual(round(region2.cp(0.0035, 700), 8), 0.208141274e1, 'Failed specific heat capacity,   state 2, region 2!')
        self.assertEqual(round(region2.w(0.0035, 700), 6),  0.644289068e3, 'Failed speed of sound,           state 2, region 2!')

        self.assertEqual(round(region2.v(30, 700), 11), 0.542946619e-2,    'Failed specific volume,          state 3, region 2!')
        self.assertEqual(round(region2.u(30, 700), 5),  0.246861076e4,     'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(round(region2.s(30, 700), 8),  0.517540298e1,     'Failed specific entropy,         state 3, region 2!')
        self.assertEqual(round(region2.h(30, 700), 5),  0.263149474e4,     'Failed specific enthalpy,        state 3, region 2!')
        self.assertEqual(round(region2.cp(30, 700), 7), 0.103505092e2,     'Failed specific heat capacity,   state 3, region 2!')
        self.assertEqual(round(region2.w(30, 700), 6),  0.480386523e3,     'Failed speed of sound,           state 3, region 2!')
    def test_ThermodynamicProperty_Region3_bnd23(self):
        self.assertEqual(round(region3.bnd23T(0.165291643e2), 6), 0.623150000e3, 'Failed boundary equation, region 23 by T!')
        self.assertEqual(round(region3.bnd23P(0.623150000e3), 7), 0.165291643e2, 'Failed boundary equation, region 23 by P!') 
    def test_ThermodynamicProperty_Region3(self):
        self.assertEqual(round(region3.P(1 / 500, 650), 7),  0.255837018e2, 'Failed pressure,                 state 1, region 3!')
        self.assertEqual(round(region3.u(1 / 500, 650), 5),  0.181226279e4, 'Failed specific internal energy, state 1, region 3!')
        self.assertEqual(round(region3.s(1 / 500, 650), 8),  0.405427273e1, 'Failed specific entropy,         state 1, region 3!')
        self.assertEqual(round(region3.h(1 / 500, 650), 5),  0.186343019e4, 'Failed specific enthalpy,        state 1, region 3!')
        self.assertEqual(round(region3.cp(1 / 500, 650), 7), 0.138935717e2, 'Failed specific heat capacity,   state 1, region 3!')
        self.assertEqual(round(region3.w(1 / 500, 650), 6),  0.502005554e3, 'Failed speed of sound,           state 1, region 3!')

        self.assertEqual(round(region3.P(1 / 200, 650), 7),  0.222930643e2, 'Failed pressure,                 state 2, region 3!')
        self.assertEqual(round(region3.u(1 / 200, 650), 5),  0.226365868e4, 'Failed specific internal energy, state 2, region 3!')
        self.assertEqual(round(region3.s(1 / 200, 650), 8),  0.485438792e1, 'Failed specific entropy,         state 2, region 3!')
        self.assertEqual(round(region3.h(1 / 200, 650), 5),  0.237512401e4, 'Failed specific enthalpy,        state 2, region 3!')
        self.assertEqual(round(region3.cp(1 / 200, 650), 7), 0.446579342e2, 'Failed specific heat capacity,   state 2, region 3!')
        self.assertEqual(round(region3.w(1 / 200, 650), 6),  0.383444594e3, 'Failed speed of sound,           state 2, region 3!')

        self.assertEqual(round(region3.P(1 / 500, 750), 7), 0.783095639e2,  'Failed pressure,                 state 3, region 3!')
        self.assertEqual(round(region3.u(1 / 500, 750), 5),  0.210206932e4, 'Failed specific internal energy, state 3, region 3!')
        self.assertEqual(round(region3.s(1 / 500, 750), 8),  0.446971906e1, 'Failed specific entropy,         state 3, region 3!')
        self.assertEqual(round(region3.h(1 / 500, 750), 5),  0.225868845e4, 'Failed specific enthalpy,        state 3, region 3!')
        self.assertEqual(round(region3.cp(1 / 500, 750), 8), 0.634165359e1, 'Failed specific heat capacity,   state 3, region 3!')
        self.assertEqual(round(region3.w(1 / 500, 750), 6),  0.760696041e3, 'Failed speed of sound,           state 3, region 3!')
    def test_ThermodynamicProperty_Region4_satP(self):
        self.assertEqual(round(region4.satP(300), 11), 0.353658941e-2, 'Failed satuation pressure, 300K!') 
        self.assertEqual(round(region4.satP(500), 8),  0.263889776e1,  'Failed satuation pressure, 500K!') 
        self.assertEqual(round(region4.satP(600), 7),  0.123443146e2,  'Failed satuation pressure, 600K!') 
    def test_ThermodynamicProperty_Region4_satT(self):
        self.assertEqual(round(region4.satT(0.10), 6), 0.372755919e3,  'Failed satuation pressure, 0.1 MPa!') 
        self.assertEqual(round(region4.satT(1.00), 6), 0.453035632e3,  'Failed satuation pressure, 1.0 MPa!') 
        self.assertEqual(round(region4.satT(10.0), 6), 0.584149488e3,  'Failed satuation pressure, 10  MPa!') 

class test_ThermodynamicPropertyBackwards(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_Ph(self):
        self.assertEqual(round(region1.T_h(3, 500), 6),   0.391798509e3, 'Failed backward temperature, state 1, region 1!') 
        self.assertEqual(round(region1.T_h(80, 500), 6),  0.378108626e3, 'Failed backward temperature, state 2, region 1!')
        self.assertEqual(round(region1.T_h(80, 1500), 6), 0.611041229e3, 'Failed backward temperature, state 3, region 1!')

        self.assertAlmostEqual(region1.v(3, 300)  / region1.v_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed v consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.u(3, 300)  / region1.u_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed u consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.s(3, 300)  / region1.s_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed s consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.g(3, 300)  / region1.g_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed g consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.cp(3, 300) / region1.cp_h(3, region1.h(3, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.cv(3, 300) / region1.cv_h(3, region1.h(3, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.w(3, 300)  / region1.w_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed w consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.a(3, 300)  / region1.a_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed a consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.k(3, 300)  / region1.k_h(3,  region1.h(3, 300)), 1.000, places=2, msg='Failed k consistancy,  state 1, region 1!')

        self.assertAlmostEqual(region1.v(80, 300)  / region1.v_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed v consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.u(80, 300)  / region1.u_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed u consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.s(80, 300)  / region1.s_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed s consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.g(80, 300)  / region1.g_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed g consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.cp(80, 300) / region1.cp_h(80, region1.h(80, 300)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.cv(80, 300) / region1.cv_h(80, region1.h(80, 300)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.w(80, 300)  / region1.w_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed w consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.a(80, 300)  / region1.a_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed a consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.k(80, 300)  / region1.k_h(80,  region1.h(80, 300)), 1.000, places=2, msg='Failed k consistancy,  state 2, region 1!')
        
        self.assertAlmostEqual(region1.v(3, 500)  / region1.v_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed v consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.u(3, 500)  / region1.u_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed u consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.s(3, 500)  / region1.s_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed s consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.g(3, 500)  / region1.g_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed g consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.cp(3, 500) / region1.cp_h(3, region1.h(3, 500)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.cv(3, 500) / region1.cv_h(3, region1.h(3, 500)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.w(3, 500)  / region1.w_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed w consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.a(3, 500)  / region1.a_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed a consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.k(3, 500)  / region1.k_h(3,  region1.h(3, 500)), 1.000, places=2, msg='Failed k consistancy,  state 3, region 1!')
    def test_ThermodynamicProperty_Region1_Ps(self):
        self.assertEqual(round(region1.T_s(3,  0.5), 6), 0.307842258e3, 'Failed backward temperature, state 1, region 1!') 
        self.assertEqual(round(region1.T_s(80, 0.5), 6), 0.309979785e3, 'Failed backward temperature, state 2, region 1!')
        self.assertEqual(round(region1.T_s(80, 3.0), 6), 0.565899909e3, 'Failed backward temperature, state 3, region 1!')

        self.assertAlmostEqual(region1.v(3, 300)  / region1.v_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed v consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.u(3, 300)  / region1.u_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed u consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.h(3, 300)  / region1.h_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed s consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.g(3, 300)  / region1.g_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed g consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.cp(3, 300) / region1.cp_s(3, region1.s(3, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.cv(3, 300) / region1.cv_s(3, region1.s(3, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 1!')
        self.assertAlmostEqual(region1.w(3, 300)  / region1.w_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed w consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.a(3, 300)  / region1.a_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed a consistancy,  state 1, region 1!')
        self.assertAlmostEqual(region1.k(3, 300)  / region1.k_s(3,  region1.s(3, 300)), 1.000, places=2, msg='Failed k consistancy,  state 1, region 1!')

        self.assertAlmostEqual(region1.v(80, 300)  / region1.v_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed v consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.u(80, 300)  / region1.u_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed u consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.h(80, 300)  / region1.h_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed s consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.g(80, 300)  / region1.g_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed g consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.cp(80, 300) / region1.cp_s(80, region1.s(80, 300)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.cv(80, 300) / region1.cv_s(80, region1.s(80, 300)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 1!')
        self.assertAlmostEqual(region1.w(80, 300)  / region1.w_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed w consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.a(80, 300)  / region1.a_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed a consistancy,  state 2, region 1!')
        self.assertAlmostEqual(region1.k(80, 300)  / region1.k_s(80,  region1.s(80, 300)), 1.000, places=2, msg='Failed k consistancy,  state 2, region 1!')
        
        self.assertAlmostEqual(region1.v(3, 500)  / region1.v_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed v consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.u(3, 500)  / region1.u_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed u consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.h(3, 500)  / region1.h_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed s consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.g(3, 500)  / region1.g_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed g consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.cp(3, 500) / region1.cp_s(3, region1.s(3, 500)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.cv(3, 500) / region1.cv_s(3, region1.s(3, 500)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 1!')
        self.assertAlmostEqual(region1.w(3, 500)  / region1.w_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed w consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.a(3, 500)  / region1.a_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed a consistancy,  state 3, region 1!')
        self.assertAlmostEqual(region1.k(3, 500)  / region1.k_s(3,  region1.s(3, 500)), 1.000, places=2, msg='Failed k consistancy,  state 3, region 1!')
    def test_ThermodynamicProperty_Region2_Ph(self):
        self.assertEqual(round(region2.bnd2b2c(0.100e3), 6), 0.3516004323e4, 'Failed boundary equation, region 2b-2c!')  
        self.assertEqual(round(region2.T_h(0.001, 3000), 6), 0.534433241e3,  'Failed backward temperature, state 1, region 2a!')
        self.assertEqual(round(region2.T_h(3, 3000), 6),     0.575373370e3,  'Failed backward temperature, state 2, region 2a!') 
        self.assertEqual(round(region2.T_h(3, 4000), 5),     0.101077577e4,  'Failed backward temperature, state 3, region 2a!')         
        self.assertEqual(round(region2.T_h(5, 3500), 6),     0.801299102e3,  'Failed backward temperature, state 1, region 2b!')
        self.assertEqual(round(region2.T_h(5, 4000), 5),     0.101531583e4,  'Failed backward temperature, state 2, region 2b!') 
        self.assertEqual(round(region2.T_h(25, 3500), 6),    0.875279054e3,  'Failed backward temperature, state 3, region 2b!')
        self.assertEqual(round(region2.T_h(40, 2700), 6),    0.743056411e3,  'Failed backward temperature, state 1, region 2c!')
        self.assertEqual(round(region2.T_h(60, 2700), 6),    0.791137067e3,  'Failed backward temperature, state 2, region 2c!') 
        self.assertEqual(round(region2.T_h(60, 3200), 6),    0.882756860e3,  'Failed backward temperature, state 3, region 2c!')

        self.assertAlmostEqual(region2.v(0.0035, 300)  / region2.v_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed v consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 300)  / region2.u_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed u consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.s(0.0035, 300)  / region2.s_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed s consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 300)  / region2.g_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=1, msg='Failed g consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 300) / region2.cp_h(0.0035, region2.h(0.0035, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 300) / region2.cv_h(0.0035, region2.h(0.0035, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 300)  / region2.w_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed w consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 300)  / region2.a_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed a consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 300)  / region2.k_h(0.0035,  region2.h(0.0035, 300)), 1.000, places=2, msg='Failed k consistancy,  state 1, region 2!')

        self.assertAlmostEqual(region2.v(0.0035, 700)  / region2.v_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed v consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 700)  / region2.u_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed u consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.s(0.0035, 700)  / region2.s_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed s consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 700)  / region2.g_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed g consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 700) / region2.cp_h(0.0035, region2.h(0.0035, 700)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 700) / region2.cv_h(0.0035, region2.h(0.0035, 700)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 700)  / region2.w_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed w consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 700)  / region2.a_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed a consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 700)  / region2.k_h(0.0035,  region2.h(0.0035, 700)), 1.000, places=2, msg='Failed k consistancy,  state 2, region 2!')
        
        self.assertAlmostEqual(region2.v(30, 700)  / region2.v_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed v consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.u(30, 700)  / region2.u_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed u consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.s(30, 700)  / region2.s_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed s consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.g(30, 700)  / region2.g_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed g consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.cp(30, 700) / region2.cp_h(30, region2.h(30, 700)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.cv(30, 700) / region2.cv_h(30, region2.h(30, 700)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.w(30, 700)  / region2.w_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed w consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.a(30, 700)  / region2.a_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed a consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.k(30, 700)  / region2.k_h(30,  region2.h(30, 700)), 1.000, places=2, msg='Failed k consistancy,  state 3, region 2!')
    def test_ThermodynamicProperty_Region2_Ps(self):
        self.assertEqual(round(region2.T_s(0.1, 7.5), 6), 0.399517097e3,  'Failed backward temperature, state 1, region 2a!')
        self.assertEqual(round(region2.T_s(0.1, 8.0), 6), 0.514127081e3,  'Failed backward temperature, state 2, region 2a!') 
        self.assertEqual(round(region2.T_s(2.5, 8.0), 5), 0.103984917e4,  'Failed backward temperature, state 3, region 2a!')         
        self.assertEqual(round(region2.T_s(8.0, 6.0), 6), 0.600484040e3,  'Failed backward temperature, state 1, region 2b!')
        self.assertEqual(round(region2.T_s(8.0, 7.5), 5), 0.106495556e4,  'Failed backward temperature, state 2, region 2b!') 
        self.assertEqual(round(region2.T_s(90, 6.0), 5),  0.103801126e4,  'Failed backward temperature, state 3, region 2b!')
        self.assertEqual(round(region2.T_s(20, 5.75), 6), 0.697992849e3,  'Failed backward temperature, state 1, region 2c!')
        self.assertEqual(round(region2.T_s(80, 5.25), 6), 0.854011484e3,  'Failed backward temperature, state 2, region 2c!') 
        self.assertEqual(round(region2.T_s(80, 5.75), 6), 0.949017998e3,  'Failed backward temperature, state 3, region 2c!')

        self.assertAlmostEqual(region2.v(0.0035, 300)  / region2.v_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed v consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 300)  / region2.u_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed u consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.h(0.0035, 300)  / region2.h_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed s consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 300)  / region2.g_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=1, msg='Failed g consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 300) / region2.cp_s(0.0035, region2.s(0.0035, 300)), 1.000, places=2, msg='Failed cp consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 300) / region2.cv_s(0.0035, region2.s(0.0035, 300)), 1.000, places=2, msg='Failed cv consistancy, state 1, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 300)  / region2.w_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed w consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 300)  / region2.a_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed a consistancy,  state 1, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 300)  / region2.k_s(0.0035,  region2.s(0.0035, 300)), 1.000, places=2, msg='Failed k consistancy,  state 1, region 2!')

        self.assertAlmostEqual(region2.v(0.0035, 700)  / region2.v_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed v consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.u(0.0035, 700)  / region2.u_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed u consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.h(0.0035, 700)  / region2.h_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed s consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.g(0.0035, 700)  / region2.g_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed g consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.cp(0.0035, 700) / region2.cp_s(0.0035, region2.s(0.0035, 700)), 1.000, places=2, msg='Failed cp consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.cv(0.0035, 700) / region2.cv_s(0.0035, region2.s(0.0035, 700)), 1.000, places=2, msg='Failed cv consistancy, state 2, region 2!')
        self.assertAlmostEqual(region2.w(0.0035, 700)  / region2.w_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed w consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.a(0.0035, 700)  / region2.a_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed a consistancy,  state 2, region 2!')
        self.assertAlmostEqual(region2.k(0.0035, 700)  / region2.k_s(0.0035,  region2.s(0.0035, 700)), 1.000, places=2, msg='Failed k consistancy,  state 2, region 2!')
        
        self.assertAlmostEqual(region2.v(30, 700)  / region2.v_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed v consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.u(30, 700)  / region2.u_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed u consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.h(30, 700)  / region2.h_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed s consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.g(30, 700)  / region2.g_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed g consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.cp(30, 700) / region2.cp_s(30, region2.s(30, 700)), 1.000, places=2, msg='Failed cp consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.cv(30, 700) / region2.cv_s(30, region2.s(30, 700)), 1.000, places=2, msg='Failed cv consistancy, state 3, region 2!')
        self.assertAlmostEqual(region2.w(30, 700)  / region2.w_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed w consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.a(30, 700)  / region2.a_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed a consistancy,  state 3, region 2!')
        self.assertAlmostEqual(region2.k(30, 700)  / region2.k_s(30,  region2.s(30, 700)), 1.000, places=2, msg='Failed k consistancy,  state 3, region 2!')

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

class test_ThermodynamicWrapper(unittest.TestCase):
    def test_ThermodynamicWrapper_Region1_State1(self):
        self.assertEqual(h2o.v(3, 300),  region1.v(3, 300),  'Failed specific volume, state 1, region 1!')
        self.assertEqual(h2o.u(3, 300),  region1.u(3, 300),  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(h2o.s(3, 300),  region1.s(3, 300),  'Failed specific entropy, state 1, region 1!')
        self.assertEqual(h2o.h(3, 300),  region1.h(3, 300),  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(h2o.cp(3, 300), region1.cp(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.cv(3, 300), region1.cv(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.w(3, 300),  region1.w(3, 300),  'Failed speed of sound, state 1, region 1!')
        self.assertEqual(h2o.a(3, 300),  region1.a(3, 300),  'Failed expansion coefficient, state 1, region 1!')
        self.assertEqual(h2o.k(3, 300),  region1.k(3, 300),  'Failed compressability, state 1, region 1!')
    def test_ThermodynamicWrapper_Region1_State2(self):
        self.assertEqual(h2o.v(80, 300),  region1.v(80, 300),  'Failed specific volume, state 2, region 1!')
        self.assertEqual(h2o.u(80, 300),  region1.u(80, 300),  'Failed specific internal energy, state 2, region 1!')
        self.assertEqual(h2o.s(80, 300),  region1.s(80, 300),  'Failed specific entropy, state 2, region 1!')
        self.assertEqual(h2o.h(80, 300),  region1.h(80, 300),  'Failed specific enthalpy, state 2, region 1!')
        self.assertEqual(h2o.cp(80, 300), region1.cp(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.cv(80, 300), region1.cv(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.w(80, 300),  region1.w(80, 300),  'Failed speed of sound, state 2, region 1!')
        self.assertEqual(h2o.a(80, 300),  region1.a(80, 300),  'Failed expansion coefficient, state 2, region 1!')
        self.assertEqual(h2o.k(80, 300),  region1.k(80, 300),  'Failed compressability, state 2, region 1!')
    def test_ThermodynamicWrapper_Region1_State3(self):
        self.assertEqual(h2o.v(3, 500),  region1.v(3, 500),  'Failed specific volume, state 3, region 1!')
        self.assertEqual(h2o.u(3, 500),  region1.u(3, 500),  'Failed specific internal energy, state 3, region 1!')
        self.assertEqual(h2o.s(3, 500),  region1.s(3, 500),  'Failed specific entropy, state 3, region 1!')
        self.assertEqual(h2o.h(3, 500),  region1.h(3, 500),  'Failed specific enthalpy, state 3, region 1!')
        self.assertEqual(h2o.cp(3, 500), region1.cp(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.cv(3, 500), region1.cv(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.w(3, 500),  region1.w(3, 500),  'Failed speed of sound, state 3, region 1!')
        self.assertEqual(h2o.a(3, 500),  region1.a(3, 500),  'Failed expansion coefficient, state 3, region 1!')
        self.assertEqual(h2o.k(3, 500),  region1.k(3, 500),  'Failed compressability, state 3, region 1!')
    def test_ThermodynamicWrapper_Region2_State1(self):
        self.assertEqual(h2o.v(0.0035, 300),  region2.v(0.0035, 300),  'Failed specific volume, state 1, region 2!')
        self.assertEqual(h2o.u(0.0035, 300),  region2.u(0.0035, 300),  'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(h2o.s(0.0035, 300),  region2.s(0.0035, 300),  'Failed specific entropy, state 1, region 2!')
        self.assertEqual(h2o.h(0.0035, 300),  region2.h(0.0035, 300),  'Failed specific enthalpy, state 1, region 2!')
        self.assertEqual(h2o.cp(0.0035, 300), region2.cp(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.cv(0.0035, 300), region2.cv(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.w(0.0035, 300),  region2.w(0.0035, 300),  'Failed speed of sound, state 1, region 2!')
        self.assertEqual(h2o.a(0.0035, 300),  region2.a(0.0035, 300),  'Failed expansion coefficient, state 1, region 2!')
        self.assertEqual(h2o.k(0.0035, 300),  region2.k(0.0035, 300),  'Failed compressability, state 1, region 2!')
    def test_ThermodynamicWrapper_Region2_State2(self):
        self.assertEqual(h2o.v(0.0035, 700),  region2.v(0.0035, 700),  'Failed specific volume, state 2, region 2!')
        self.assertEqual(h2o.u(0.0035, 700),  region2.u(0.0035, 700),  'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(h2o.s(0.0035, 700),  region2.s(0.0035, 700),  'Failed specific entropy, state 2, region 2!')
        self.assertEqual(h2o.h(0.0035, 700),  region2.h(0.0035, 700),  'Failed specific enthalpy, state 2, region 2!')
        self.assertEqual(h2o.cp(0.0035, 700), region2.cp(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.cv(0.0035, 700), region2.cv(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.w(0.0035, 700),  region2.w(0.0035, 700),  'Failed speed of sound, state 2, region 2!')
        self.assertEqual(h2o.a(0.0035, 700),  region2.a(0.0035, 700),  'Failed expansion coefficient, state 2, region 2!')
        self.assertEqual(h2o.k(0.0035, 700),  region2.k(0.0035, 700),  'Failed compressability, state 2, region 2!')
    def test_ThermodynamicWrapper_Region2_State3(self):
        self.assertEqual(h2o.v(30, 700),  region2.v(30, 700),  'Failed specific volume, state 3, region 2!')
        self.assertEqual(h2o.u(30, 700),  region2.u(30, 700),  'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(h2o.s(30, 700),  region2.s(30, 700),  'Failed specific entropy, state 3, region 2!')
        self.assertEqual(h2o.h(30, 700),  region2.h(30, 700),  'Failed specific enthalpy, state 3, region 2!')
        self.assertEqual(h2o.cp(30, 700), region2.cp(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.cv(30, 700), region2.cv(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.w(30, 700),  region2.w(30, 700),  'Failed speed of sound, state 3, region 2!')
        self.assertEqual(h2o.a(30, 700),  region2.a(30, 700),  'Failed expansion coefficient, state 3, region 2!')
        self.assertEqual(h2o.k(30, 700),  region2.k(30, 700),  'Failed compressability, state 3, region 2!')

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
        h1 = numpy.array([[h2o.h(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])
        v1 = numpy.array([[h2o.v(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])
        u1 = numpy.array([[h2o.u(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])
        dhdp1 = numpy.array([[h2o.dhdP(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])
        dvdp1 = numpy.array([[h2o.dvdP(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])
        dudp1 = numpy.array([[h2o.dudP(P1[j, i], Ts[j], 1) for i in range(len(P1))] for j in range(len(Ts))])

        P2 = numpy.array([numpy.linspace(10**-3, (h2o.satP(n) if n <= 623.15 else region3.bnd23P(n)), len(Ts)) for n in Tl])
        T2 = numpy.array([numpy.ones(len(Ts)) * n for n in Tl])
        h2 = numpy.array([[h2o.h(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        v2 = numpy.array([[h2o.v(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        u2 = numpy.array([[h2o.u(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dhdp2 = numpy.array([[h2o.dhdP(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dvdp2 = numpy.array([[h2o.dvdP(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])
        dudp2 = numpy.array([[h2o.dudP(P2[j, i], Tl[j], 2) for i in range(len(P2[0]))] for j in range(len(Tl))])

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
        clrmp = cm.viridis
        clrmp.set_under(cm.viridis.colors[0])
        clrmp.set_over(cm.viridis.colors[-1])
        clrmp.set_bad(color='black')
        fig = pyplot.figure(figsize=(13,11))

        pyplot.subplot(321, axisbg='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dhdp1, levels=lvlh/15, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dhdp2[:], levels=lvlh, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific enthalpy w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/15th', style='italic')

        pyplot.subplot(322, axisbg='darkgrey')
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

        pyplot.subplot(323, axisbg='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dvdp1, levels=lvlv/15e9, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dvdp2[:], levels=lvlv, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific volume w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/15e9', style='italic')

        pyplot.subplot(324, axisbg='darkgrey')
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

        pyplot.subplot(325, axisbg='darkgrey')
        pyplot.semilogy(Tsat, Pfg, 'k')
        pyplot.contourf(T1, P1, dudp1, levels=lvlu/20, cmap=clrmp, extend="both")
        pyplot.contourf(T2[:], P2[:], dudp2[:], levels=lvlu, cmap=clrmp, extend="both")
        pyplot.title('Partial derivative of specific internal energy w.r.t pressure', fontsize=8)
        pyplot.xlabel('Temperature [K]')
        pyplot.ylabel('Pressure [MPa]')
        pyplot.ylim(10**-3, 50)
        pyplot.text(300, 10, 'scale 1/20th', style='italic')

        pyplot.subplot(326, axisbg='darkgrey')
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
        #pyplot.show()
        pyplot.savefig("if97/__testout__/test_PartialP", dpi=300)
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
        #pyplot.show()
        pyplot.savefig("if97/__testout__/test_Properties", dpi=300)
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
        pyplot.savefig("if97/__testout__/test_moody", dpi=300)
        #pyplot.show()
    def test_CriticalFlowModDerHEM_Plot(self):
        P0 = [0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 
              10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
        h0 = numpy.linspace(250, 2750, 50)
        y = [1, 1e9]
        n = 2

        G = [[moody.Gcritical(j, i) for i in h0] for j in P0[:n]]
        dGdP = [[moody.dGdPcritical(j, i) for i in h0] for j in P0[:n]]

        #pyplot.xkcd()
        fig = pyplot.figure(figsize=(11, 7))
        ax1 = fig.add_subplot(111)
        [pyplot.semilogy(h0, G[j], markersize=5, marker='.', label=str(P0[j])) for j in range(len(P0[:n]))]
        ax1.fill_betweenx(y, [h2o.hf(0.2), h2o.hf(0.2)], [h2o.hg(0.2), h2o.hg(0.2)], facecolor='blue')
        pyplot.ylim(10, 10**6)
        pyplot.ylabel('Critical Mass Flux [kg / s m**2]')
        pyplot.xlim(0, 3000)
        pyplot.xlabel('Stagnation Enthalpy [kJ / kg]')
        pyplot.legend(loc='lower left', ncol=5)
        pyplot.text(750, 8e1, 'Stagnation pressures [Mpa]', backgroundcolor="white", ha='center', va='top', color='black')
        pyplot.grid(linestyle="--", linewidth=0.5, color='.25', zorder=-10)
        pyplot.title('Critical Flow')

        ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position("right")
        pyplot.ylim(1e-4, 1e0)
        pyplot.ylabel('Derivative of Crit. Mass Flux [kg m**3 / s m**2 kJ]')
        [pyplot.semilogy(h0, dGdP[j], '--', markersize=5, marker='.', label=str(P0[j])) for j in range(len(P0[:n]))]
        
        #pyplot.savefig("if97/__testout__/test_moody", dpi=300)
        pyplot.show()

if __name__ == '__main__':
    unittest.main()