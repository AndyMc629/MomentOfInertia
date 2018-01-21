#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 20:03:11 2017

@author: andrew

"""
import pandas as pd

from thermodynamicanalyzer import ThermodynamicAnalyzer
from thermodynamicGas import thermodynamicGasCalculator
from vapourPressureCalculator import SaturationVapourPressureCalculator
import matplotlib.pyplot as plt
import seaborn as sns

conv=27.2114 # 1 hartree = 27.2114 eV

# the U_DFT for water below was found by deleting all the atoms inside a monohydrate unit cell and leaving the waters
# and running the optimisation with all of the same parameters.
U_dft = {'water': -0.5*1.5221135574129E+02*conv,
         'mono': -0.5*2.3951108085020E+03*conv,
         'mapi': -0.5*7.983651216778E+02*conv,
         'di': -86466.55,
         'pbi2': -2.154603426126E+02*conv}
         #
         #'mono': -0.5*7.9836512167779E+02*conv,
         #'mapi': -3.2304127634186E+02*conv}

# the U_DFT for water given below was calculated using PBESol but in a molecule (not in a cell) - this must be the
# cause of the large energy difference and then the screwy results that come from that.

#U_dft = {'water':-2040.240234654193}
if __name__ == "__main__":
    
    
    # instantiate a thermodynamic analyzer instance
    analyzer = ThermodynamicAnalyzer()
    
    # get mono data
    #df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan.out')
    df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan_1x3x1.out')
    # divide by 3 to give the unit cell total energy - didn't have to do this for MAPI so have I ran a supercell for the unti cell calcs?
    df_thermo_mono[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']] = df_thermo_mono[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']]/3

    #get mapi data
    #df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan.out')
    df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')

    # get thermodynamic water data from preprocessed file.
    df_thermo_gas = analyzer.get_gas_thermodynamic_results_nist(path='data/NIST_JANAF_WATER.dat', 
                                                U_H2O=U_dft['water'],
                                                units='eV/fu')

    # set up the gas calculator
    gas_calc = thermodynamicGasCalculator(
            H_relative_to_stp=df_thermo_gas[['T(K)','H-H(Tr)']],
            S_relative_to_stp=df_thermo_gas[['T(K)','S']],
            U_DFT=U_dft['water'] #eV/fu apparently.
            )
    

    # for this to work I have to re-run crystal calc's at temp intervals of 100K
    # as the data does not currently match with the gas data - INTERPOLATE????.
    df_thermo_gas['G_H2O(eV/fu)'] = df_thermo_gas['T(K)'].apply(lambda x: gas_calc.mu_gas(T=x,P=1))
    
    #plt.figure()
    #df_thermo_mapi.plot(x='T(K)', y='G_solid(eV/cell)', label='mapi')
    #df_thermo_mono.plot(x='T(K)', y='G_solid(eV/cell)', label='mono')
    #df_thermo_gas.plot(x='T(K)', y='G_H2O(eV/fu)', label='water')
    #plt.show()
    
    #=============================================================================
    # The chemical potential calculations for water I have performed so far
    # do not include the intermediate temperatures that I have calculated
    # phonon energies for in crystal - therefore I need to do some interpolation
    #=============================================================================

    water_vapour_regression = pd.ols(x=df_thermo_gas['T(K)'][0:7], y=df_thermo_gas['G_H2O(eV/fu)'][0:7])
    
    df_thermo_water_vapour = df_thermo_mono.copy()
    df_thermo_water_vapour.drop(
            ['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)'],
            inplace=True,
            axis=1)
    g_water_vapour_interpolated = water_vapour_regression.predict(
            beta=water_vapour_regression.beta, 
            x=df_thermo_mono['T(K)'])
    
    df_thermo_water_vapour['G_gas(eV/cell)'] = g_water_vapour_interpolated
    df_thermo_water_vapour.rename(columns={'P(MPA)': 'atm'}, inplace=True)
    df_thermo_water_vapour['atm']=1.0


    df_thermo_water_all = []
    #for pressure in [1E-7, 1E-5, 1E-3, 1E-1, 1, 1E2]:#np.arange(1.0,6.0,1.0): #atm
    for pressure in [1E-5, 1E-3, 1E-1, 1]:#, 1E2]:
    #for pressure in [3.5077E-5, 3.5077E-3, 3.5077E-2]:
        # for this to work I have to re-run crystal calc's at temp intervals of 100K
        # as the data does not currently match with the gas data - INTERPOLATE????.
        df_thermo_gas['G_H2O(eV/fu)'] = df_thermo_gas['T(K)'].apply(lambda x: gas_calc.mu_gas(T=x,P=pressure))
        
        water_vapour_regression = pd.ols(x=df_thermo_gas['T(K)'][0:7], y=df_thermo_gas['G_H2O(eV/fu)'][0:7])
    
        df_thermo_water_vapour = df_thermo_mono.copy()
    
        df_thermo_water_vapour.drop(
                ['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)'],
                inplace=True,
                axis=1)
        
        g_water_vapour_interpolated = water_vapour_regression.predict(
                beta=water_vapour_regression.beta, 
                x=df_thermo_mono['T(K)'])
        
        df_thermo_water_vapour['G_gas(eV/cell)'] = g_water_vapour_interpolated
        df_thermo_water_vapour.rename(columns={'P(MPA)': 'atm'}, inplace=True)
        df_thermo_water_vapour['atm']=pressure
        df_thermo_water_all.append(df_thermo_water_vapour)
        
    """
    
    """
    df_thermo_water_all = pd.concat(df_thermo_water_all)
    
    df_final_result_phase_diagram = df_thermo_water_all.copy()
    
    df_final_result_phase_diagram['dG(eV)'] = df_thermo_water_all.apply(
            lambda row:
                0.5*df_thermo_mono[df_thermo_mono['T(K)']==row['T(K)']]['G_solid(eV/cell)'].values[0] -
                df_thermo_mapi[df_thermo_mapi['T(K)']==row['T(K)']]['G_solid(eV/cell)'].values[0]
                - row['G_gas(eV/cell)'],
                axis=1)
        
        
    vapcalc = SaturationVapourPressureCalculator()
    df = df_final_result_phase_diagram
    plt.figure(1)
    for P in df['atm'].value_counts().index:
        if P>1E-5:
            plt.plot(df[df['atm']==P]['T(K)'], df[df['atm']==P]['dG(eV)'], label=str(P))
    plt.axhline(0)
    plt.legend()
    plt.show()
    
    #df = df[df['T(K)'].isin(range(300, 600, 100))]
    #df['P_critical_atm'] = df['T(K)'].apply(lambda x: gas_calc.critical_P_calculator_mono_with_vib(x, df_thermo_mono, df_thermo_mapi))
    #df['RH_critical_percent']=df.apply(lambda row: 100*vapcalc.convert_pressure_to_rh(row['T(K)'], row['P_critical_atm']), axis=1)
    
    #=============================================================================
    #     HEATMAP WITH SEABORN
    #=============================================================================
    #df2 = df.sort_values('atm', ascending=False)
    
    df2 = df[df['T(K)']<350]
    df2.rename(columns={'atm': 'P(atm)'}, inplace=True)
    df2 = df2.pivot('P(atm)', 'T(K)', 'dG(eV)')
    
    fig, ax = plt.subplots(figsize=(10,10))
    #sns.heatmap(data=df2, ax=ax, mask=df2>0, annot=True, center=0).invert_yaxis()
    sns.heatmap(data=df2, ax=ax, annot=True, center=0, cmap="RdYlBu").invert_yaxis()
    plt.savefig('phase_diag_heatmap.png')

    test = df[df['T(K)']<350]['T(K)'].apply(lambda x: vapcalc.calculate_water_vapour_pressure(T=x))

    # =============================================================================
    #     DISPERSION VS GAMMA THERMODYNAMIC FUNCTIONS FROM CRYSTAL.
    # =============================================================================
    df_thermo_mapi_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan_2x2x2.out')
    df_thermo_mapi = analyzer.get_solid_thermodynamic_results_crystal(path='data/mapi_freq_tscan.out')

    plt.figure(2)
    ax = df_thermo_mapi[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label=r'PBESol - $\Gamma$ - MAPbI$_3$')
    df_thermo_mapi_disp[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label='PBESol - Dispersion - MAPbI$_3$', ax=ax)
    plt.xlabel('T(K)', fontsize=16)
    plt.ylabel('G(eV/cell)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.savefig('mapi_disp_vs_gamma.png')
    plt.show()

    # =============================================================================
    #     DISPERSION VS GAMMA THERMODYNAMIC FUNCTIONS FROM CRYSTAL.
    # =============================================================================
    df_thermo_mono_disp = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan_1x3x1.out')
    df_thermo_mono = analyzer.get_solid_thermodynamic_results_crystal(path='data/mono_freq_tscan.out')

    # divide by 3 to give the unit cell total energy - didn't have to do this for MAPI so have I ran a supercell for the unti cell calcs?
    df_thermo_mono_disp[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']] = df_thermo_mono_disp[['G_solid(Hartree/cell)', 'G_solid(eV/cell)', 'G_solid(kJ/mol)']]/3

    plt.figure(3)
    ax = df_thermo_mono[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label=r'PBESol - $\Gamma$ - Mono')
    df_thermo_mono_disp[['T(K)', 'G_solid(eV/cell)']].plot(x='T(K)', y='G_solid(eV/cell)', label='PBESol - Dispersion - Mono', ax=ax)
    plt.xlabel('T(K)', fontsize=16)
    plt.ylabel('G(eV/cell)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.savefig('mono_disp_vs_gamma.png')
    plt.close()

    # divide by two to get the Gibbs free energy per formula unit for the water vapour
    print 0.5*df_thermo_mono_disp['G_solid(eV/cell)']-df_thermo_mapi_disp['G_solid(eV/cell)']-df_thermo_water_vapour['G_gas(eV/cell)']

    plt.figure(4)
    plt.plot(df_thermo_mono_disp['T(K)'].values,
             (0.5 * df_thermo_mono_disp['G_solid(eV/cell)']-\
              df_thermo_mapi_disp['G_solid(eV/cell)']-\
              df_thermo_water_vapour['G_gas(eV/cell)']).values,
             label=r'$\Delta G_{Mono}$ - PBESol - Dispersion')

    plt.plot(df_thermo_mono_disp['T(K)'].values,
             (0.5 * df_thermo_mono['G_solid(eV/cell)']-\
              df_thermo_mapi['G_solid(eV/cell)']-\
              df_thermo_water_vapour['G_gas(eV/cell)']).values,
             label=r'$\Delta G_{Mono}$ - PBESol - $\Gamma$')
    plt.axhline(0, color='red', linestyle='--')
    plt.legend(loc='upper left')
    plt.xlabel('T(K)', fontsize=16)
    plt.ylabel(r'$\Delta G_{Mono}$', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.savefig('deltaG_mono_dispersion_deltaG_dispersion_vs_gamma.png')
    plt.show()



    plt.figure(5)
    comparison_pressure_atm = 0.1
    plt.plot(df_thermo_mono_disp['T(K)'].values,
             (0.5 * df_thermo_mono_disp['G_solid(eV/cell)']-\
              df_thermo_mapi_disp['G_solid(eV/cell)']-\
              df_thermo_water_vapour['G_gas(eV/cell)']).values,
             label=r'$\Delta G_{Mono}$ - PBESol - Dispersion - $P_{H_2O}$ = 1 atm')

    plt.plot(df_thermo_mono_disp['T(K)'].values,
             (0.5 * df_thermo_mono['G_solid(eV/cell)']-\
              df_thermo_mapi['G_solid(eV/cell)']-\
              df_thermo_water_all[df_thermo_water_all['atm']==comparison_pressure_atm]['G_gas(eV/cell)']).values,
             label=r'$\Delta G_{Mono}$ - PBESol - Dispersion - $P_{H_2O}$ = '+str(comparison_pressure_atm)+' atm')
    plt.axhline(0, color='red', linestyle='--')
    plt.legend(loc='upper left', fontsize=12)
    plt.xlabel('T(K)', fontsize=16)
    plt.ylabel(r'$\Delta G_{Mono}$', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid()
    plt.tight_layout()
    plt.savefig('deltaG_mono_dispersion_two_pressures.png')
    plt.show()



    # =============================================================================
    #
    # =============================================================================
    df_thermo_crystal_water = analyzer.get_gas_thermodynamic_results_crystal(path='data/water_freq_tscan.out')
    # divide by two to get the Gibbs free energy per formula unit for the water vapour
    print 0.5*df_thermo_mono_disp['G_solid(eV/cell)']-df_thermo_mapi_disp['G_solid(eV/cell)']-0.5*df_thermo_crystal_water['G_solid(eV/cell)']