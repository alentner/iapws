import unittest
from if97 import region1, region4, conversion

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
class test_UnitConversion(unittest.TestCase):
    def test_Conversion_Temperature_0K(self):
        self.assertEqual(round(conversion.temp_K_C(0), 2), -273.15, 'Failed Kelvin to Celsius, Absolute zero!') 
        self.assertEqual(round(conversion.temp_C_F(conversion.temp_K_C(0)), 2), -459.67, 'Failed Kelvin to Fahrenheit, Absolute zero!')
        self.assertEqual(round(conversion.temp_F_R(conversion.temp_C_F(conversion.temp_K_C(0))), 2), 0.00, 'Failed Kelvin to Rankine, Absolute zero!')
    def test_Conversion_Temperature_0F(self):
        self.assertEqual(round(conversion.temp_C_K(conversion.temp_F_C(0)), 2), 255.37, 'Failed Fahrenheit to Kelvin , Freezing point of brine!')
        self.assertEqual(round(conversion.temp_F_C(0), 2), -17.78, 'Failed Fahrenheit to Celsius, Freezing point of brine!') 
        self.assertEqual(round(conversion.temp_F_R(0), 2), 459.67, 'Failed Fahrenheit to Rankine, Freezing point of brine!')
    def test_Conversion_Temperature_0C(self):
        self.assertEqual(round(conversion.temp_C_K(0), 2), 273.15, 'Failed Celsius to Kelvin , Freezing point of water!')
        self.assertEqual(round(conversion.temp_C_F(0), 2), 32.00, 'Failed Celsius to Fahrenheit, Freezing point of water!') 
        self.assertEqual(round(conversion.temp_F_R(conversion.temp_C_F(0)), 2), 491.67, 'Failed Celsius to Rankine, Freezing point of water!')
    def test_Conversion_Temperature_100C(self):
        self.assertEqual(round(conversion.temp_C_K(conversion.temp_F_C(conversion.temp_R_F(671.64102))), 2), 373.13, 'Failed Rankine to Kelvin , Boiling point of water!')
        self.assertEqual(round(conversion.temp_R_F(671.64102), 2), 211.97, 'Failed Rankine to Fahrenheit, Boiling point of water!')
    def test_Conversion_Pressure(self):
        self.assertEqual(round(conversion.press_psia_MPa(1300), 2), 8.96, 'Failed psia to Mpa, Atmosphere of Venus!')
        self.assertEqual(round(conversion.press_MPa_psia(0.08), 2), 11.60, 'Failed MPa to psia, Vacuum cleaner!')

if __name__ == '__main__':
    unittest.main()
