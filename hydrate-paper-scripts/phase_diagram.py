"""
5/12/2017        crystal_results_file.close()

        dftp = pd.DataFrame(df_TP)
        dfg = pd.DataFrame(df_G)
        dfet = pd.DataFrame(df_ET)

        df_thermo_solid = pd.concat([dftp, dfg, dfet], axis=1, join='inner')
        #df_thermo_solid = pd.concat([dftp, dfg], axis=1, join='inner')
        self.df_thermo_solid = df_thermo_solid.astype(float)
        #print df_thermo_solid
        return df_thermo_solid.astype(float)

Create a phase diagram from the calculated crystal data.
"""


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from thermodynamicanalyzer import ThermodynamicAnalyzer
from thermodynamicGas import thermodynamicGasCalculator
from vapourPressureCalculator import SaturationVapourPressureCalculator

if __name__=="__main__":
    # instantiate a thermodynamic analyzer instance
    analyzer = ThermodynamicAnalyzer()

    df_thermo_mapi_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')
    df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan.out')
    df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan.out')
    df_thermo_water = analyzer.get_gas_thermodynamic_results_crystal(path='data/water_opt.out')

    dg = 0.5*df_thermo_mono['G_solid(eV/cell)'] - df_thermo_mapi['G_solid(eV/cell)'] - df_thermo_water['G_solid(eV/cell)']
    dg_kj = 0.5*df_thermo_mono['G_solid(kJ/mol)'] - df_thermo_mapi['G_solid(kJ/mol)'] - df_thermo_water['G_solid(kJ/mol)']


    dg_disp = 0.5*df_thermo_mono['G_solid(eV/cell)'] - df_thermo_mapi_disp['G_solid(eV/cell)'] - df_thermo_water['G_solid(eV/cell)']


    # =============================================================================
    #     DISPERSION VS GAMMA THERMODYNAMIC FUNCTIONS FROM CRYSTAL.
    # =============================================================================
    df_thermo_mapi_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')
    df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan.out')

    plt.figure(1)
    ax = df_thermo_mapi[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label=r'PBESol - $\Gamma$ - MAPbI$_3$')
    df_thermo_mapi_disp[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label='PBESol - Dispersion - MAPbI$_3$', ax=ax)
    plt.xlabel('T(K)')
    plt.ylabel('G(eV/cell)')
    plt.savefig('mapi_disp_vs_gamma.png')
    plt.show()


    # =============================================================================
    #     DISPERSION VS GAMMA THERMODYNAMIC FUNCTIONS FROM CRYSTAL.
    # =============================================================================
    df_thermo_mono_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan_1x3x1.out')
    df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan.out')

    # divide by 3 to give the unit cell total energy - didn't have to do this for MAPI so have I ran a supercell for the unti cell calcs?
    df_thermo_mono_disp[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']] = df_thermo_mono_disp[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']]/3

    plt.figure(2)
    ax = df_thermo_mono[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label=r'PBESol - $\Gamma$ - Mono')
    df_thermo_mono_disp[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label='PBESol - Dispersion - Mono', ax=ax)
    plt.xlabel('T(K)')
    plt.ylabel('G(eV/cell)')
    plt.savefig('mono_disp_vs_gamma.png')
    #plt.show()

    # divide by two to get the Gibbs free energy per formula unit for the water vapour
    #print 0.5*df_thermo_mono_disp['G_solid(eV/cell)']-df_thermo_mapi_disp['G_solid(eV/cell)']-df_thermo_water_vapour['G_gas(eV/cell)']
    plt.figure(3)
    plt.plot(df_thermo_mono_disp['T(K)'].values, (
    0.5 * df_thermo_mono_disp['G_solid(eV/cell)'] - df_thermo_mapi_disp['G_solid(eV/cell)'] - df_thermo_water[
        'G_solid(eV/cell)']).values)
    plt.show()