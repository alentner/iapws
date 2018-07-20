import unittest
from if97 import region1, region2, region4
from if97 import region1_ext, region2_ext, region4_ext
from if97 import h2o

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
    def test_ThermodynamicProperty_Region4_satP(self):
        self.assertEqual(round(region4.satP(300), 11), 0.353658941e-2, 'Failed satuation pressure, 300K!') 
        self.assertEqual(round(region4.satP(500), 8),  0.263889776e1, 'Failed satuation pressure, 500K!') 
        self.assertEqual(round(region4.satP(600), 7),  0.123443146e2, 'Failed satuation pressure, 600K!') 
    def test_ThermodynamicProperty_Region4_satT(self):
        self.assertEqual(round(region4.satT(0.10), 6), 0.372755919e3, 'Failed satuation pressure, 0.1 MPa!') 
        self.assertEqual(round(region4.satT(1.00), 6), 0.453035632e3, 'Failed satuation pressure, 1.0 MPa!') 
        self.assertEqual(round(region4.satT(10.0), 6), 0.584149488e3, 'Failed satuation pressure, 10 MPa!') 

class test_ThermodynamicProperty_wrapper(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_State1(self):
        self.assertEqual(h2o.v(3, 300),  region1.v(3, 300),  'Failed specific volume, state 1, region 1!')
        self.assertEqual(h2o.u(3, 300),  region1.u(3, 300),  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(h2o.s(3, 300),  region1.s(3, 300),  'Failed specific entropy, state 1, region 1!')
        self.assertEqual(h2o.h(3, 300),  region1.h(3, 300),  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(h2o.cp(3, 300), region1.cp(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.cv(3, 300), region1.cv(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.w(3, 300),  region1.w(3, 300),  'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region1_State2(self):
        self.assertEqual(h2o.v(80, 300),  region1.v(80, 300),  'Failed specific volume, state 2, region 1!')
        self.assertEqual(h2o.u(80, 300),  region1.u(80, 300),  'Failed specific internal energy, state 2, region 1!')
        self.assertEqual(h2o.s(80, 300),  region1.s(80, 300),  'Failed specific entropy, state 2, region 1!')
        self.assertEqual(h2o.h(80, 300),  region1.h(80, 300),  'Failed specific enthalpy, state 2, region 1!')
        self.assertEqual(h2o.cp(80, 300), region1.cp(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.cv(80, 300), region1.cv(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.w(80, 300),  region1.w(80, 300),  'Failed speed of sound, state 1, region 2!')
    def test_ThermodynamicProperty_Region1_State3(self):
        self.assertEqual(h2o.v(3, 500),  region1.v(3, 500),  'Failed specific volume, state 3, region 1!')
        self.assertEqual(h2o.u(3, 500),  region1.u(3, 500),  'Failed specific internal energy, state 3, region 1!')
        self.assertEqual(h2o.s(3, 500),  region1.s(3, 500),  'Failed specific entropy, state 3, region 1!')
        self.assertEqual(h2o.h(3, 500),  region1.h(3, 500),  'Failed specific enthalpy, state 3, region 1!')
        self.assertEqual(h2o.cp(3, 500), region1.cp(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.cv(3, 500), region1.cv(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.w(3, 500),  region1.w(3, 500),  'Failed speed of sound, state 1, region 3!')
    def test_ThermodynamicProperty_Region2_State1(self):
        self.assertEqual(h2o.v(0.0035, 300),  region2.v(0.0035, 300),  'Failed specific volume, state 1, region 2!')
        self.assertEqual(h2o.u(0.0035, 300),  region2.u(0.0035, 300),  'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(h2o.s(0.0035, 300),  region2.s(0.0035, 300),  'Failed specific entropy, state 1, region 2!')
        self.assertEqual(h2o.h(0.0035, 300),  region2.h(0.0035, 300),  'Failed specific enthalpy, state 1, region 2!')
        self.assertEqual(h2o.cp(0.0035, 300), region2.cp(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.cv(0.0035, 300), region2.cv(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.w(0.0035, 300),  region2.w(0.0035, 300),  'Failed speed of sound, state 1, region 2!')
    def test_ThermodynamicProperty_Region2_State2(self):
        self.assertEqual(h2o.v(0.0035, 700),  region2.v(0.0035, 700),  'Failed specific volume, state 2, region 2!')
        self.assertEqual(h2o.u(0.0035, 700),  region2.u(0.0035, 700),  'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(h2o.s(0.0035, 700),  region2.s(0.0035, 700),  'Failed specific entropy, state 2, region 2!')
        self.assertEqual(h2o.h(0.0035, 700),  region2.h(0.0035, 700),  'Failed specific enthalpy, state 2, region 2!')
        self.assertEqual(h2o.cp(0.0035, 700), region2.cp(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.cv(0.0035, 700), region2.cv(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.w(0.0035, 700),  region2.w(0.0035, 700),  'Failed speed of sound, state 2, region 2!')
    def test_ThermodynamicProperty_Region2_State3(self):
        self.assertEqual(h2o.v(30, 700),  region2.v(30, 700),  'Failed specific volume, state 3, region 2!')
        self.assertEqual(h2o.u(30, 700),  region2.u(30, 700),  'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(h2o.s(30, 700),  region2.s(30, 700),  'Failed specific entropy, state 3, region 2!')
        self.assertEqual(h2o.h(30, 700),  region2.h(30, 700),  'Failed specific enthalpy, state 3, region 2!')
        self.assertEqual(h2o.cp(30, 700), region2.cp(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.cv(30, 700), region2.cv(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.w(30, 700),  region2.w(30, 700),  'Failed speed of sound, state 3, region 2!')

class test_ThermodynamicProperty_ext(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_T(self):
        self.assertEqual(round(region1_ext.T(3, 500), 6),   0.391798509e3, 'Failed backward temperature, state 1, region 1!') 
        self.assertEqual(round(region1_ext.T(80, 500), 6),  0.378108626e3, 'Failed backward temperature, state 2, region 1!')
        self.assertEqual(round(region1_ext.T(80, 1500), 6), 0.611041229e3, 'Failed backward temperature, state 3, region 1!')
    def test_ThermodynamicProperty_Region1_v(self):
        self.assertEqual(round(region1_ext.v(3, region1.h(3, 300)), 9),   round(region1.v(3, region1_ext.T(3, region1.h(3, 300))), 9),    'Failed specific volume consistancy, state 1, region 1!')
        self.assertEqual(round(region1_ext.v(80, region1.h(80, 300)), 9), round(region1.v(80, region1_ext.T(80, region1.h(80, 300))), 9), 'Failed specific volume consistancy, state 2, region 1!')
        self.assertEqual(round(region1_ext.v(3, region1.h(3, 500)), 9),   round(region1.v(3, region1_ext.T(3, region1.h(3, 500))), 9),    'Failed specific volume consistancy, state 3, region 1!')
    def test_ThermodynamicProperty_Region2_Boundary2b2c(self):
        self.assertEqual(round(region2_ext.bnd2b2c(0.100e3), 6), 0.3516004323e4, 'Failed boundary equation, region 2b-2c!')    
    def test_ThermodynamicProperty_Region2_T(self):
        self.assertEqual(round(region2_ext.T(0.001, 3000), 6), 0.534433241e3, 'Failed backward temperature, state 1, region 2a!')
        self.assertEqual(round(region2_ext.T(3, 3000), 6),     0.575373370e3, 'Failed backward temperature, state 2, region 2a!') 
        self.assertEqual(round(region2_ext.T(3, 4000), 5),     0.101077577e4, 'Failed backward temperature, state 3, region 2a!')         
        self.assertEqual(round(region2_ext.T(5, 3500), 6),     0.801299102e3, 'Failed backward temperature, state 1, region 2b!')
        self.assertEqual(round(region2_ext.T(5, 4000), 5),     0.101531583e4, 'Failed backward temperature, state 2, region 2b!') 
        self.assertEqual(round(region2_ext.T(25, 3500), 6),    0.875279054e3, 'Failed backward temperature, state 3, region 2b!')
        self.assertEqual(round(region2_ext.T(40, 2700), 6),    0.743056411e3, 'Failed backward temperature, state 1, region 2c!')
        self.assertEqual(round(region2_ext.T(60, 2700), 6),    0.791137067e3, 'Failed backward temperature, state 2, region 2c!') 
        self.assertEqual(round(region2_ext.T(60, 3200), 6),    0.882756860e3, 'Failed backward temperature, state 3, region 2c!')
    def test_ThermodynamicProperty_Region2_v(self):
        self.assertEqual(round(region2_ext.v(0.0035, region2.h(0.0035, 300)), 9), round(region2.v(0.0035, region2_ext.T(0.0035, region2.h(0.0035, 300))), 9), 'Failed specific volume consistancy, state 1, region 2!')
        self.assertEqual(round(region2_ext.v(0.0035, region2.h(0.0035, 700)), 9), round(region2.v(0.0035, region2_ext.T(0.0035, region2.h(0.0035, 700))), 9), 'Failed specific volume consistancy, state 2, region 2!')
        self.assertEqual(round(region2_ext.v(30, region2.h(30, 700)), 9),         round(region2.v(30, region2_ext.T(30, region2.h(30, 700))), 9),             'Failed specific volume consistancy, state 3, region 2!')

if __name__ == '__main__':
    unittest.main()
