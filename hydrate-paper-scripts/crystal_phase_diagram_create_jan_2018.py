#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 20/1/2018

@author: andrew

Script to calculate and plot monohydrate phase diagram using crystal results.
"""

import pandas as pd

from thermodynamicanalyzer import ThermodynamicAnalyzer
from thermodynamicGas import thermodynamicGasCalculator
from vapourPressureCalculator import SaturationVapourPressureCalculator
import matplotlib.pyplot as plt
import numpy as np

conv=27.2114 # 1 hartree = 27.2114 eV

# the U_DFT for water below was found by deleting all the atoms inside a monohydrate unit cell and leaving the waters
# and running the optimisation with all of the same parameters.
"""
U_dft = {#'water': -0.5*1.5221135574129E+02*conv, # mono water not converged properly!
         #'water': -2040.240234654193, # molecule water not converged properly!
         'water': -7.6365717564858E+01*conv,
         'mono': -0.5*2.3951108085020E+03*conv,
         'mapi': -0.5*7.983651216778E+02*conv,
         'di': -86466.55,
         'pbi2': -2.154603426126E+02*conv}
"""
U_dft={#'water': -7.636646941130E+01*conv, # molecule tightly optimised
       'water': -7.610672734852E+01*conv, #molecule optimised to same standard as other crystals
       #'water': -0.5*1.522122897990E+02*conv, # in the monohydrate crystal
      'mono': -0.5*7.983651378652E+02*conv,
      'mapi': -3.230413016696E+02*conv}

# using crystal's internal conversion factor, not much diff at all.
"""
U_dft2 = {'water':-2070.969337199977,
         'mono': -0.5*21724.619859643008,
          'mapi': -8790.400713767818}
"""
# Uses Antoine equation values from
# NIST Chemweb-book.
# NB: Antoine equation returns
def SaturationVapourPressure(T):
    if T > 255.9 and T <= 373:  # not sure about equals?
        A = 4.6543
        B = 1435.264
        C = -64.848
    elif T >= 379 and T < 573:  # not sure about equals?
        A = 3.55959
        B = 643.748
        C = -198.043
    elif T>373 and T<379:  # need something better than setting to zero = wrong ...
        A = (4.6543+3.55959)/2  # 0
        B = (1435.264+643.748)/2  # 0
        C = (-64.848-198.043)/2  # 0
    elif T<=255.9:
        A = 4.6543
        B = 1435.264
        C = -64.848
    elif T>=573:
        A = 3.55959
        B = 643.748
        C = -198.043

    log10SatVapPress = A-(B/(T+C))
    # print "Psat = %f, T= %f" % (10**(log10SatVapPress), T)
    return 10**(log10SatVapPress)
    #return log10SatVapPress

import math
kB = 8.617E-5 #eV/K

def critical_p_line(T, dFT, dMuT):
    """

    :param T: numpy array
    :param dFT: numpy array
    :param dMuT: numpy array
    :return: numpy array of p values corresponding to the T's passed in.
    """
    T = kB*T
    dU = U_dft['mono']-U_dft['water']-U_dft['mapi']

    return np.exp( 1.0*(dU + dFT - dMuT)/T )


if __name__ == "__main__":

    analyzer = ThermodynamicAnalyzer()
    df_thermo_mono_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan_1x3x1.out')
    df_thermo_mapi_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')
    #df_thermo_water_crystal = analyzer.get_gas_thermodynamic_results_crystal(path='data/water_freq_tscan.out')
    df_thermo_water_crysta = analyzer.get_gas_thermodynamic_results_crystal(path='data/water_opt_molecule.out')
    # supercell
    df_thermo_mono_disp[['G_solid(eV/cell)', 'ET(eV/cell)']] = df_thermo_mono_disp[['G_solid(eV/cell)','ET(eV/cell)']] / 3

    # 2 f.u's in the crystal calc
    df_thermo_mono_disp[['G_solid(eV/cell)', 'ET(eV/cell)']] = df_thermo_mono_disp[['G_solid(eV/cell)','ET(eV/cell)']].apply(lambda x: 0.5 * x)

    # get thermodynamic water data from preprocessed file.
    df_thermo_gas = analyzer.get_gas_thermodynamic_results_nist(path='data/NIST_JANAF_WATER.dat',
                                                                U_H2O=U_dft['water'],
                                                                units='eV/fu')
    # set up the gas calculator
    gas_calc = thermodynamicGasCalculator(
        H_relative_to_stp=df_thermo_gas[['T(K)', 'H-H(Tr)']],
        S_relative_to_stp=df_thermo_gas[['T(K)', 'S']],
        U_DFT=U_dft['water']  # eV/fu apparently.
    )

    # for this to work I have to re-run crystal calc's at temp intervals of 100K
    # as the data does not currently match with the gas data - INTERPOLATE????.
    df_thermo_gas['dMu(T)'] = df_thermo_gas['T(K)'].apply(lambda x: gas_calc.dMu(x))

    T = df_thermo_mono_disp['T(K)'].values[1::2]
    Fmono = df_thermo_mono_disp['ET(eV/cell)'].values[1::2]
    Fmapi = df_thermo_mapi_disp['ET(eV/cell)'].values[1::2]
    dFT = Fmono-Fmapi
    U_water = (U_dft['water'])*np.ones(len(T))
    dMuT = df_thermo_gas[df_thermo_gas['T(K)'].isin(range(100,800,100))]['dMu(T)'].values

    critical_p = critical_p_line(T, dFT, dMuT)

    p_sat = [SaturationVapourPressure(x) for x in T]

    T_p = zip(T, critical_p, p_sat)
    critical_RH = [100.0*x[1]/x[2] for x in T_p]

    #plt.figure(1)
    #plt.plot(T, dGT-U_water-dMuT)
    #plt.close()

    #THIS LOOKS LIKE A WINNER WINNER CHICKEN DINNER!
    # Note - it's a bit blocky but I can sort that with some more T values.
    plt.figure(2)
    plt.ylim(0,100)
    plt.xlim(200,700)
    plt.plot(T, critical_RH)
    plt.show()