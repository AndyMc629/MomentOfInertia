#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 15:14:25 2017

@author: andrew


thermodynamicAnalyzer class for analysing crystal results
"""

import re
import pandas as pd
import numpy as np
from thermodynamicGas import UnitConvertor

class ThermodynamicAnalyzer(object):
    
    def __init__(self, U_dft = None):
    
        if U_dft is not None:
            self.U_solid = U_dft['solid']
            self.U_gas = U_dft['gas']
        
        else:
            self.U_solid = None
            self.U_gas = None
            
        self.df_thermo_solid = None
        self.df_thermo_gas = None
        
        # The information about the linear extrapolation of the free energy for the gas.
        self.gas_gibbs_energy_regression = None
        
        # for conversion of units
        self.unitConvertor = UnitConvertor()
        
    def get_solid_thermodynamic_results_crystal(self, path=None):
        
        if path==None or not isinstance(path, str):
            raise ValueError('require path to a .out file to read')
        else:
            pass
        
        crystal_results_file = open(path, "r")
    
        df_G = []
        df_ET = []
        df_TP = []
        
        for line in crystal_results_file:
            if re.search("T =  ", line):
                 data = re.sub("[^0-9\+.E]", " ", line).split()
                 dictionary_TP = {'T(K)': data[0], 'P(MPA)': data[1]}
                 df_TP.append(dictionary_TP)

            # total Gibbs
            if re.search("E0\+", line) and not re.search("ZERO-POINT", line):#:#+E0+ET+PV-TS", line):
                data = line.split()
                print data
                if len(data) == 5:
                    pass
                else:
                    dictionary_G = {'G_solid(Hartree/cell)': data[1],
                                    'G_solid(eV/cell)': data[2],
                                    'G_solid(kJ/mol)': data[3]}
                    df_G.append(dictionary_G)

            # Thermal contribution
            if re.search("ET  ", line) and not (re.search("THERMAL", line) or re.search("NET", line)):
                data = line.split()
                dictionary_ET = {'ET(Hartree/cell)': data[2],
                                    'ET(eV/cell)': data[3],
                                    'ET(kJ/mol)': data[4]}
                df_ET.append(dictionary_ET)


        crystal_results_file.close()
        
        dftp = pd.DataFrame(df_TP)
        dfg = pd.DataFrame(df_G)
        dfet = pd.DataFrame(df_ET)

        df_thermo_solid = pd.concat([dftp, dfg, dfet], axis=1, join='inner')
        #df_thermo_solid = pd.concat([dftp, dfg], axis=1, join='inner')
        self.df_thermo_solid = df_thermo_solid.astype(float)
        #print df_thermo_solid
        return df_thermo_solid.astype(float)

    def get_gas_thermodynamic_results_crystal(self, path=None):
        if path == None or not isinstance(path, str):
            raise ValueError('require path to a .out file to read')
        else:
            pass

        crystal_results_file = open(path, "r")

        df_G = []
        df_TP = []

        for line in crystal_results_file:
            if re.search("T =  ", line):
                data = re.sub("[^0-9\+.E]", " ", line).split()
                dictionary_TP = {'T(K)': data[0], 'P(MPA)': data[1]}
                df_TP.append(dictionary_TP)

            # total Gibbs
            if re.search("E0\+", line) and not re.search("ZERO-POINT", line):  #:#+E0+ET+PV-TS", line):
                data = line.split()
                print data
                if len(data) == 5:
                    pass
                else:
                    dictionary_G = {'G_solid(Hartree/cell)': data[1],
                                    'G_solid(eV/cell)': data[2],
                                    'G_solid(kJ/mol)': data[3]}
                    df_G.append(dictionary_G)
        crystal_results_file.close()

        dftp = pd.DataFrame(df_TP)
        dfg = pd.DataFrame(df_G)


        df_thermo_gas = pd.concat([dftp, dfg], axis=1, join='inner')
        # self.df_thermo_solid = df_thermo_gas.astype(float)
        # print df_thermo_solid
        return df_thermo_gas.astype(float)

    def get_gas_thermodynamic_results_nist(self, path=None, U_H2O=None, units=None):
        
        if U_H2O is None and self.U_gas is None:
            raise ValueError("Must supply a T=0K energy value (eV/fu) to this method.")
        
        elif U_H2O is None and self.U_gas is not None:
            U_H2O = self.U_gas
        
        else: #i.e E0 is not none
            pass
        
        if path is None:
            raise ValueError('require a path to downloaded NIST-JANAF gas data')
            
        if units is None:
            raise ValueError("""must supply units for U_H2O, please choose one of:
                'eV/fu',
                'kJ/mol',
                'Hartree/fu'""")
        
        # I am assuming a very specific format here from the data I found.
        data = pd.read_csv(path, 
                       sep='\t', 
                       skiprows=1, 
                       dtype=np.float64,
                       na_values='INFINITE')
    
        # need the values where the free energy is linear, 5-11 obv empirical
        # for the nist water data.
        regression=pd.ols(x=data['T(K)'][5:11], y=data['-[G-H(Tr)]/T'][5:11])
        
        # keep regression result as class attribute
        self.gas_gibbs_energy_regression = regression
        
        # make the extrapolation
        trend = regression.predict(beta=regression.beta, x=data['T(K)'][0:100])
        
        #save the extrapolated data, think we need a good old minus here?
        data['-G_regression(J/K/mol)'] = trend
        
        #plot to check that it looks correct.
        #data.plot(x='T(K)', y=['-[G-H(Tr)]/T', '-G_regression(J/K/mol)'], xlim=[0,1000])
        
        # convert to more useful units (eV/fu)
        # CLOSELY INSPECT NEXT LINE IF THINGS START LOOKING WRONG
        data['G_regression(eV/fu)'] = data.apply(
                lambda row: self.unitConvertor.convert(
                        -1*row['-G_regression(J/K/mol)']*row['T(K)'], 
                        currentUnits='J/mol',
                        newUnits='eV/fu'),
                        axis=1)
        
        self.df_thermo_gas = data
        
        return data.astype(float)
        
        
        
        
        
        
        
            
    
