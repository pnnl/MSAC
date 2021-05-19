import pandas as pd
import numpy as np


def calculate_total_mz(adduct_tuple, mass):
    '''
    calculate adduct mz using multiplier, charge, and mass
    '''
    input_mz_multiplier = adduct_tuple[0]
    adduct_charge = adduct_tuple[1]
    adduct_mz = adduct_tuple[2]
    #  calculating as adduct mz + ((input multiplier*input mass)/charge)
    total_mz = adduct_mz + input_mz_multiplier*mass/np.absolute(adduct_charge)
    return total_mz


def calculate_all_mz(df, mass_file, mass_col):
    input_masses = pd.read_csv(mass_file)

    #  create a lookup table for the adduct information:
    #  name, input mass mulitpler, charge, m/z
    d = {adduct: [mult, charge, mass] for adduct, (mult, (charge, mass))
         in zip(df['adduct'], zip(df['input_mass_multiplier'],
                zip(df['charge'], df['m/z'])))}

    masses_to_calc = input_masses[mass_col]
    original_cols = input_masses.columns
    all_masses = []
    for adduct in d.keys():
        #  put pos or neg on adduct name
        #  to designate electrospray ionization state
        if int(d[adduct][1]) > 0:
            adduct_name = adduct + '_ESIpos'
        else:
            adduct_name = adduct + '_ESIneg'
        input_masses[adduct_name] = [calculate_total_mz(d[adduct], mass)
                                     for mass in masses_to_calc]
    
    input_masses = pd.melt(input_masses, id_vars=original_cols, var_name='adduct', value_name='adduct mass')

    return input_masses
