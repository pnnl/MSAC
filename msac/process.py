from os import path
import numpy as np
import pandas as pd
import pkg_resources

from msac import calculate_adduct_mz, calculate_input_mz

def process_file(input_mass, mass_col = None, no_mass_formula_col = None, adduct_file = None, include_neutral_loss = False, outname = None, coverage_cutoff = 1.0, restrict = None):
    # calculate adduct mz
    if adduct_file:
        df = calculate_adduct_mz.calculate_adduct_mz(adduct_file, None)
        print("Using supplied adduct file {}. Coverage cutoff not used.".format(adduct_file))
    else:
        if include_neutral_loss:
            ADDUCT_FILE = pkg_resources.resource_filename('msac',
                                                      'example_data/adduct_list_full.csv')
        else: # don't include neutral losses, adducts only
            ADDUCT_FILE = pkg_resources.resource_filename('msac',
                                                      'example_data/adduct_only_list.csv')
        df = calculate_adduct_mz.calculate_adduct_mz(ADDUCT_FILE, coverage_cutoff)

    # load masses if not dataframe
    if not isinstance(input_mass, pd.DataFrame):
        try:
            input_masses = pd.read_csv(input_mass)
        except Exception as e:
            print("Please provide input masses as either a dataframe with a 'Mass' column or a valid filename of a .csv")
            print("Error was", e)
    else:
        input_masses = input_mass
        input_mass = 'input_masses.csv'

    # if no_mass is specified, calculate monoisotopic mass from formula
    if no_mass_formula_col:
        input_masses = calculate_input_mz.calculate_mass_from_formula(input_masses, no_mass_formula_col)

    # calculate input mass mz for each adduct and add adduct mz
    output = calculate_input_mz.calculate_all_mz(df, input_masses,
                                                      mass_col, restrict)


    if outname:
        output_name = outname
    else:
        inname, ext = path.splitext(input_mass)
        output_name = inname + '_adducts.csv'
    output.to_csv(output_name, index=False)
    return output