#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 14:30:01 2017

@author: andrew

Vapour pressure calculator 
"""
import math
from thermodynamicGas import UnitConvertor

convertor = UnitConvertor()

class SaturationVapourPressureCalculator(object):
    
    def calculate_water_vapour_pressure(self, T=None, units='atm'): # temp in Kelvin
        """
        Uses the Antoine equation to calculate the vapour pressure at a given 
        temperature, originally in bar so conversion is performed.
        
        """
        A,B,C = self.get_ABC(T=T)
        
        if A is not None and B is not None and C is not None:
            # bar 
            p_vap_bar = math.pow(10, (A-B/(C+T)))
            if units=='bar':
                return p_vap_bar
            
            # atm
            elif units=='atm':    
                p_vap_atm = convertor.convert(
                        p_vap_bar, 
                        currentUnits='bar', 
                        newUnits='atm')
                return p_vap_atm
            
            else:
                return None
        else:
            return None
    
    def get_ABC(self, T=None):
        """
        From: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4#Thermo-Phase
        T(K)        ||    A    ||  B    ||   C    || Ref                   ||Comment
        379. - 573.	3.55959	643.748	-198.043	Liu and Lindsay, 1970	Coefficents calculated by NIST from author's data.
        273. - 303.	5.40221	1838.675	-31.737 Bridgeman and Aldrich, 1964	Coefficents calculated by NIST from author's data.
        304. - 333.	5.20389	1733.926	-39.485 Bridgeman and Aldrich, 1964	Coefficents calculated by NIST from author's data.
        334. - 363.	5.0768	1659.793	-45.854	Bridgeman and Aldrich, 1964	Coefficents calculated by NIST from author's data.
        344. - 373.	5.08354	1663.125	-45.622	Bridgeman and Aldrich, 1964	Coefficents calculated by NIST from author's data.
        293. - 343.	6.20963	2354.731	7.559	Gubkov, Fermor, et al., 1964	Coefficents calculated by NIST from author's data.
        255.9 - 373.	4.6543	1435.264	-64.848	Stull, 1947	Coefficents calculated by NIST from author's data.
        """
        A = None
        B = None
        C = None
        # quickest way to get values for most temperatures.
        
        if (T>=255.0 and T<=373.0):
            """ Stull, 1947 """
            A = 4.6543	
            B = 1435.264
            C = -64.848
            
        elif (T>=379.0 and T<=573.0):
            """ Liu and Lindasay, 1970 """
            A = 3.55959
            B = 643.748
            C = -198.043
        
        return A,B,C
    
    def convert_pressure_to_rh(self, T,P): #assume p in atm
        return P/self.calculate_water_vapour_pressure(T)