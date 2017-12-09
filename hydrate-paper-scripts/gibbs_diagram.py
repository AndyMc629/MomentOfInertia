#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 13:55:21 2017

@author: andrew


Phase diag test

"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from thermodynamicanalyzer import ThermodynamicAnalyzer
from thermodynamicGas import thermodynamicGasCalculator
from vapourPressureCalculator import SaturationVapourPressureCalculator

conv=27.2114 # 1 hartree = 27.2114 eV

U_dft = {'water': -0.5*1.5221135574129E+02*conv,
         'mono': -0.5*7.9836512167779E+02*conv,
         'mapi': -3.2304127634186E+02*conv} 


if __name__ == "__main__":
    # instantiate a thermodynamic analyzer instance
    analyzer = ThermodynamicAnalyzer()
    
    # get mono data
    df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan_1x3x1.out')
    df_thermo_mono = df_thermo_mono[['T(K)', 'G_solid(eV/cell)']]
    # supercell
    df_thermo_mono['G_solid(eV/cell)'] = df_thermo_mono['G_solid(eV/cell)']/3
    
    # 2 f.u's in the crystal calc
    df_thermo_mono['G_solid(eV/cell)'] = df_thermo_mono['G_solid(eV/cell)'].apply(lambda x: 0.5*x)
    
    #get mapi data
    df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')
    df_thermo_mapi = df_thermo_mapi[['T(K)', 'G_solid(eV/cell)']]
    # used a supercell but don't divide it seems, need to check this.

    
    # get thermodynamic water data from preprocessed file.
    df_thermo_gas = analyzer.get_gas_thermodynamic_results_nist(path='data/NIST_JANAF_WATER.dat', 
                                                U_H2O=U_dft['water'],
                                                units='eV/fu')

    
    gas_calc = thermodynamicGasCalculator(
            H_relative_to_stp=df_thermo_gas[['T(K)','H-H(Tr)']],
            S_relative_to_stp=df_thermo_gas[['T(K)','S']],
            U_DFT=U_dft['water'] #eV/fu apparently.
            )
    
    
    vapcalc = SaturationVapourPressureCalculator()
    
    ts = range(100,750,100)
    df_thermo_mono = df_thermo_mono[df_thermo_mono['T(K)'].isin(ts)]
    df_thermo_mapi = df_thermo_mono[df_thermo_mono['T(K)'].isin(ts)]
    df_thermo_gas_final=df_thermo_mono.copy()
    
    df_thermo_gas_final['G_solid(eV/cell)'] = df_thermo_gas_final['T(K)'].apply(lambda x: gas_calc.mu_gas(T=x, P=0.03507718399736391))
    
    
    
    # =============================================================================
    #  G = E_mono - E_mapi - mu    
    # =============================================================================
    x = np.arange(-1,0.1,0.1)
    
    plt.figure()
    plt.plot(x, U_dft['mono']-U_dft['mapi']-U_dft['water']-x)
    plt.xlabel(r'$\Delta \mu(T, p)$')
    plt.ylabel(r'$\Delta G (eV)$')
    plt.show()
    
    
    # =============================================================================
    #     
    # =============================================================================
    
    
    
    
    """
    import numpy as np 
    from pandas import DataFrame
    import seaborn as sns
    %matplotlib inline
    
    Index= ['aaa', 'bbb', 'ccc', 'ddd', 'eee']
    Cols = ['A', 'B', 'C', 'D']
    df = DataFrame(abs(np.random.randn(5, 4)), index=Index, columns=Cols)
    
    sns.heatmap(df, annot=True)
    
    
    #OR
    
    flights = sns.load_dataset("flights")
    flights = flights.pivot("month", "year", "passengers")
    ax = sns.heatmap(flights, annot=True, fmt="d")
    """
    
    
    """
    # for this to work I have to re-run crystal calc's at temp intervals of 100K
    # as the data does not currently match with the gas data - INTERPOLATE????.
    df_thermo_gas['G_H2O(eV/fu)'] = df_thermo_gas['T(K)'].apply(lambda x: gas_calc.mu_gas(T=x,P=1))
    
    df_test=df_thermo_gas[df_thermo_gas['T(K)'].isin([300, 400, 500])][['T(K)', 'G_H2O(eV/fu)']]  
    
    vapcalc = SaturationVapourPressureCalculator()
    
    df_test['G_H2O(eV/fu)_sat_pressure'] = df_test['T(K)'].apply(
            lambda x: gas_calc.mu_gas(T=x,P=vapcalc.calculate_water_vapour_pressure(x))
            )
    
    df_test['dG_sat_pressure'] = df_test['G_H2O(eV/fu)_sat_pressure'].apply(
            lambda x: U_dft['mono'] - U_dft['mapi'] - x)
    
    def convert_pressure_to_rh(T,P): #assume p in atm
        return P/vapcalc.calculate_water_vapour_pressure(T)
    
    df_test['RH_crit']=df_test['T(K)'].apply(lambda x: convert_pressure_to_rh(T=x, P=gas_calc.critical_P_calculator_mono(T=x)))
    """

    
    
    
    