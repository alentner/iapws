import unittest
from if97 import region1, region4

class test_ThermodynamicProperty(unittest.TestCase):
    def test_ThermodynamicProperty_Region1_state1(self):
        self.assertEqual(round(region1.v(3, 300), 11), 0.100215168E-2, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.h(3, 300), 6), 0.115331273E3, 'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.u(3, 300), 6), 0.112324818E3, 'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(3, 300), 9), 0.392294792, 'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.cp(3, 300), 8), 0.417301218E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(3, 300), 8), 0.412120160E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(3, 300), 5), 0.150773921E4, 'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region1_state2(self):
        self.assertEqual(round(region1.v(80, 300), 12), 0.971180894E-3, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.h(80, 300), 6), 0.184142828E3, 'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.u(80, 300), 6), 0.106448356E3, 'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(80, 300), 9), 0.368563852, 'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.cp(80, 300), 8), 0.401008987E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(80, 300), 8), 0.391736606E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(80, 300), 5), 0.163469054E4, 'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region1_state3(self):
        self.assertEqual(round(region1.v(3, 500), 11), 0.120241800E-2, 'Failed specific volume, state 1, region 1!')
        self.assertEqual(round(region1.h(3, 500), 6), 0.975542239E3, 'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(round(region1.u(3, 500), 6), 0.971934985E3, 'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(round(region1.s(3, 500), 8), 0.258041912E1, 'Failed specific entropy, state 1, region 1!')
        self.assertEqual(round(region1.cp(3, 500), 8), 0.465580682E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.cv(3, 500), 8), 0.322139223E1, 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(round(region1.w(3, 500), 5), 0.124071337E4, 'Failed speed of sound, state 1, region 1!')
    def test_ThermodynamicProperty_Region4_satP(self):
        self.assertEqual(round(region4.satP(300), 11), 0.353658941E-2, 'Failed satuation pressure, 300K!') 
        self.assertEqual(round(region4.satP(500), 8), 0.263889776E1, 'Failed satuation pressure, 500K!') 
        self.assertEqual(round(region4.satP(600), 7), 0.123443146E2, 'Failed satuation pressure, 600K!') 
    def test_ThermodynamicProperty_Region4_satT(self):
        self.assertEqual(round(region4.satT(0.10), 6), 0.372755919E3, 'Failed satuation pressure, 0.1 MPa!') 
        self.assertEqual(round(region4.satT(1.00), 6), 0.453035632E3, 'Failed satuation pressure, 1.0 MPa!') 
        self.assertEqual(round(region4.satT(10.0), 6), 0.584149488E3, 'Failed satuation pressure, 10 MPa!') 

if __name__ == '__main__':
    unittest.main()
