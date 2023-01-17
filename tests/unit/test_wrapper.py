import unittest
from iapws.if97 import region1, region2, region3, region4, h2o

class test_ThermodynamicWrapper(unittest.TestCase):

    def test_ThermodynamicWrapper_Region1_State1(self):
        self.assertEqual(h2o.v(3, 300),  region1.v(3, 300),  'Failed specific volume, state 1, region 1!')
        self.assertEqual(h2o.u(3, 300),  region1.u(3, 300),  'Failed specific internal energy, state 1, region 1!')
        self.assertEqual(h2o.s(3, 300),  region1.s(3, 300),  'Failed specific entropy, state 1, region 1!')
        self.assertEqual(h2o.h(3, 300),  region1.h(3, 300),  'Failed specific enthalpy, state 1, region 1!')
        self.assertEqual(h2o.cp(3, 300), region1.cp(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.cv(3, 300), region1.cv(3, 300), 'Failed specific heat capacity, state 1, region 1!')
        self.assertEqual(h2o.w(3, 300),  region1.w(3, 300),  'Failed speed of sound, state 1, region 1!')
        self.assertEqual(h2o.av(3, 300), region1.av(3, 300), 'Failed expansion coefficient, state 1, region 1!')
        self.assertEqual(h2o.kT(3, 300), region1.kT(3, 300), 'Failed compressability, state 1, region 1!')

    def test_ThermodynamicWrapper_Region1_State2(self):
        self.assertEqual(h2o.v(80, 300),  region1.v(80, 300),  'Failed specific volume, state 2, region 1!')
        self.assertEqual(h2o.u(80, 300),  region1.u(80, 300),  'Failed specific internal energy, state 2, region 1!')
        self.assertEqual(h2o.s(80, 300),  region1.s(80, 300),  'Failed specific entropy, state 2, region 1!')
        self.assertEqual(h2o.h(80, 300),  region1.h(80, 300),  'Failed specific enthalpy, state 2, region 1!')
        self.assertEqual(h2o.cp(80, 300), region1.cp(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.cv(80, 300), region1.cv(80, 300), 'Failed specific heat capacity, state 2, region 1!')
        self.assertEqual(h2o.w(80, 300),  region1.w(80, 300),  'Failed speed of sound, state 2, region 1!')
        self.assertEqual(h2o.av(80, 300), region1.av(80, 300), 'Failed expansion coefficient, state 2, region 1!')
        self.assertEqual(h2o.kT(80, 300), region1.kT(80, 300), 'Failed compressability, state 2, region 1!')

    def test_ThermodynamicWrapper_Region1_State3(self):
        self.assertEqual(h2o.v(3, 500),  region1.v(3, 500),  'Failed specific volume, state 3, region 1!')
        self.assertEqual(h2o.u(3, 500),  region1.u(3, 500),  'Failed specific internal energy, state 3, region 1!')
        self.assertEqual(h2o.s(3, 500),  region1.s(3, 500),  'Failed specific entropy, state 3, region 1!')
        self.assertEqual(h2o.h(3, 500),  region1.h(3, 500),  'Failed specific enthalpy, state 3, region 1!')
        self.assertEqual(h2o.cp(3, 500), region1.cp(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.cv(3, 500), region1.cv(3, 500), 'Failed specific heat capacity, state 3, region 1!')
        self.assertEqual(h2o.w(3, 500),  region1.w(3, 500),  'Failed speed of sound, state 3, region 1!')
        self.assertEqual(h2o.av(3, 500), region1.av(3, 500), 'Failed expansion coefficient, state 3, region 1!')
        self.assertEqual(h2o.kT(3, 500), region1.kT(3, 500), 'Failed compressability, state 3, region 1!')

    def test_ThermodynamicWrapper_Region2_State1(self):
        self.assertEqual(h2o.v(0.0035, 300),  region2.v(0.0035, 300),  'Failed specific volume, state 1, region 2!')
        self.assertEqual(h2o.u(0.0035, 300),  region2.u(0.0035, 300),  'Failed specific internal energy, state 1, region 2!')
        self.assertEqual(h2o.s(0.0035, 300),  region2.s(0.0035, 300),  'Failed specific entropy, state 1, region 2!')
        self.assertEqual(h2o.h(0.0035, 300),  region2.h(0.0035, 300),  'Failed specific enthalpy, state 1, region 2!')
        self.assertEqual(h2o.cp(0.0035, 300), region2.cp(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.cv(0.0035, 300), region2.cv(0.0035, 300), 'Failed specific heat capacity, state 1, region 2!')
        self.assertEqual(h2o.w(0.0035, 300),  region2.w(0.0035, 300),  'Failed speed of sound, state 1, region 2!')
        self.assertEqual(h2o.av(0.0035, 300), region2.av(0.0035, 300), 'Failed expansion coefficient, state 1, region 2!')
        self.assertEqual(h2o.kT(0.0035, 300), region2.kT(0.0035, 300), 'Failed compressability, state 1, region 2!')

    def test_ThermodynamicWrapper_Region2_State2(self):
        self.assertEqual(h2o.v(0.0035, 700),  region2.v(0.0035, 700),  'Failed specific volume, state 2, region 2!')
        self.assertEqual(h2o.u(0.0035, 700),  region2.u(0.0035, 700),  'Failed specific internal energy, state 2, region 2!')
        self.assertEqual(h2o.s(0.0035, 700),  region2.s(0.0035, 700),  'Failed specific entropy, state 2, region 2!')
        self.assertEqual(h2o.h(0.0035, 700),  region2.h(0.0035, 700),  'Failed specific enthalpy, state 2, region 2!')
        self.assertEqual(h2o.cp(0.0035, 700), region2.cp(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.cv(0.0035, 700), region2.cv(0.0035, 700), 'Failed specific heat capacity, state 2, region 2!')
        self.assertEqual(h2o.w(0.0035, 700),  region2.w(0.0035, 700),  'Failed speed of sound, state 2, region 2!')
        self.assertEqual(h2o.av(0.0035, 700), region2.av(0.0035, 700), 'Failed expansion coefficient, state 2, region 2!')
        self.assertEqual(h2o.kT(0.0035, 700), region2.kT(0.0035, 700), 'Failed compressability, state 2, region 2!')

    def test_ThermodynamicWrapper_Region2_State3(self):
        self.assertEqual(h2o.v(30, 700),  region2.v(30, 700),  'Failed specific volume, state 3, region 2!')
        self.assertEqual(h2o.u(30, 700),  region2.u(30, 700),  'Failed specific internal energy, state 3, region 2!')
        self.assertEqual(h2o.s(30, 700),  region2.s(30, 700),  'Failed specific entropy, state 3, region 2!')
        self.assertEqual(h2o.h(30, 700),  region2.h(30, 700),  'Failed specific enthalpy, state 3, region 2!')
        self.assertEqual(h2o.cp(30, 700), region2.cp(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.cv(30, 700), region2.cv(30, 700), 'Failed specific heat capacity, state 3, region 2!')
        self.assertEqual(h2o.w(30, 700),  region2.w(30, 700),  'Failed speed of sound, state 3, region 2!')
        self.assertEqual(h2o.av(30, 700), region2.av(30, 700), 'Failed expansion coefficient, state 3, region 2!')
        self.assertEqual(h2o.kT(30, 700), region2.kT(30, 700), 'Failed compressability, state 3, region 2!')

if __name__ == '__main__':
    unittest.main()
