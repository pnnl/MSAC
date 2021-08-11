from os import path
import numpy as np
import pandas as pd
import pkg_resources

import calculate_adduct_mz
import calculate_input_mz

def process_file(input_mass_file, mass_col = None, adduct_file = None, outname = None, coverage_cutoff = None, restrict = None):
    # calculate adduct mz
    if adduct_file:
        df = calculate_adduct_mz.calculate_adduct_mz(adduct_file, None)
        print("Using supplied adduct file {}. Coverage cutoff not used.".format(adduct_file))
    else:
        ADDUCT_FILE = pkg_resources.resource_filename('msac',
                                                      'example_data/adduct_list_full.csv')
        df = calculate_adduct_mz.calculate_adduct_mz(ADDUCT_FILE, coverage_cutoff)

    # if formula given, remove adducts if they can't be lost

    # calculate input mass mz for each adduct and add adduct mz
    output = calculate_input_mz.calculate_all_mz(df, input_masses,
                                                      mass_col, restrict)
    print(output.columns)
    print(output.head())

    if outname:
        output_name = outname
    else:
        inname, ext = path.splitext(input_masses)
        output_name = inname + '_adducts.csv'
    output.to_csv(output_name, index=False)