#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 17:38:56 2017

@author: andrew

thermodynamicGasCalculation
"""

import math


# key water data taken from CODATA tables
# http://www.science.uwaterloo.ca/~cchieh/cact/tools/thermodata.html
key_water_data_CODATA = {'S_stp(J/K/mol)': 188.835,
                  'H_stp - H_sp_0K(kJ/mol)': 9.905}



class UnitConvertor(object):
    """
    Class for conversion between frequently used chemical units.
    """
    
    def convert(self, value, currentUnits, newUnits):
        
        newValue = None
        
        # 1.0365327133842577e-05 = conversion for J/mol to eV/fu
        
        if currentUnits=='J/mol' and newUnits=='eV/fu':
            newValue = value*1.0365327133842577e-05
        
        if currentUnits=='Hartree/fu' and newUnits=='eV/fu':
            newValue = value*27.2114
        
        if currentUnits=='bar' and newUnits=='atm':
            newValue = value*0.986923
            
        return newValue
            
            
   
conv=27.2114 # 1 hartree = 27.2114 eV

U_dft = {'water': -0.5*1.5221135574129E+02*conv,
         'mono': -0.5*7.9836512167779E+02*conv,
         'mapi': -3.2304127634186E+02*conv}
     
class thermodynamicGasCalculator(object):
    """
    Class for calculation of the thermodynamic functions of gas (water vapour)
    on the scale of performed DFT calculations.
    """
    
    def __init__(self, H_relative_to_stp, S_relative_to_stp, U_DFT=None):
        if U_DFT is None:
            self.U_DFT = -471 #eV/fu - I guess this was taken from a crystal calc?
        else:
            self.U_DFT = U_DFT
            
        self.kB = 8.617E-5 #eV/K
        
        # instantiate a unit convertor for possible later use
        self.unitConvertor = UnitConvertor()
        
        # convert the zero temp enthalpy offset from CODATA tables to eV/fu 
        self.H_zero_kelvin_offset = self.unitConvertor.convert(
                value=key_water_data_CODATA['H_stp - H_sp_0K(kJ/mol)']*1000, #kJ/mol->J/mol
                currentUnits = 'J/mol',
                newUnits = 'eV/fu'
                )
        # enthalpy of water vapour given as a function of temperature with stp
        # as the thermodynamic reference point.
        self.H_relative_to_stp = H_relative_to_stp
        
        # entropy of water vapour given as a function of temperature with stp
        # as the thermodynamic reference point.
        self.S_relative_to_stp = S_relative_to_stp
        
    def dMu(self, T):
        """
        Method to calculate the temperature dependent offset of the chemical potential
        for water vapour on the DFT scale.
        """
        # assume H and S are pandas dataframes, taken from Nist
        
        # get the enthalpy from kJ/mol to J/mol
        dMu = (self.H_relative_to_stp[self.H_relative_to_stp['T(K)']==T]['H-H(Tr)']*1000).values[0] # J/mol
        
        # subtract away the entropy term at the given temp
        dMu = dMu - T*self.S_relative_to_stp[self.S_relative_to_stp['T(K)']==T]['S'].values[0] # J/mol
        
        # everything was in J/mol up until now so convert dMu to eV/fu
        dMu = self.unitConvertor.convert(value=dMu, currentUnits='J/mol', newUnits='eV/fu')
        
        # add the previously calculated (see constructor of this class) zero 
        # kelvin enthalpy offset.
        dMu = dMu + self.H_zero_kelvin_offset
        
        return dMu
        
    def mu_gas(self, T, P):
        """
        Method to calculate the final temperature and pressure dependent
        chemical potential of water vapour on the DFT scale (assuming ideal gas
        water vapour as always in this analysis).
        """
        return self.U_DFT + self.dMu(T) + self.kB*T*math.log(P) #assume p0=1 atm.
   
    
    def critical_P_calculator_mono(self, T):
        """
        Critical P = phase boundary value of P for a given temp
        """
        enthalpy = U_dft['mono'] - U_dft['mapi'] - U_dft['water']
        #free_energy_vib =  
        
        return math.exp( (enthalpy-self.dMu(T))/(self.kB*T))
    
    def critical_P_calculator_mono_with_vib(self, T, F_mono, F_mapi):
        """
        Critical P = phase boundary value of P for a given temp
        """
        free_energy_vib = F_mono[F_mono['T(K)']==T]['G_solid(eV/cell)']-\
        F_mapi[F_mapi['T(K)']==T]['G_solid(eV/cell)']- U_dft['water']
        print free_energy_vib.values[0]
        
        return math.exp( (free_energy_vib.values[0]-self.dMu(T))/(self.kB*T))

    
        