import pandas as pd
import numpy as np

def calculate_total_mz(df, mass_file, mass_col):
    input_masses = pd.read_csv(mass_file)
    
    # create a lookup table for the adduct information: name, input mass mulitpler, charge, m/z
    d = {adduct:[mult, charge, mass] for adduct, (mult, (charge, mass)) in zip(df['adduct'], zip(df['input_mass_multiplier'], zip(df['charge'], df['m/z'])))}

    
    masses_to_calc = input_masses[mass_col]
    all_masses = []
    for adduct in d.keys():
        # calculating as: adduct mz + ((input multiplier*input mass)/charge)
        input_masses[adduct] = [d[adduct][2] + ((d[adduct][0]*mass)/np.absolute(d[adduct][1])) for mass in masses_to_calc]
    
    return input_masses