import argparse
from os import path
import numpy as np
import pandas as pd
import pkg_resources

import msac


def main():
    p = {}

    parser = argparse.ArgumentParser(description='Calculates the m/z of potential adducts for a given compound m/z')
    parser.add_argument('-v', '--version', action='version', version=msac.__version__, help='print version and exit')

    parser.add_argument('input_masses', help=".csv with 'mass' column")
    parser.add_argument('-f', '--adduct_file', help="path to a .csv with an 'adduct' col and a 'charge' col. Defaults to 'example_data/adduct_list.csv' ")
    parser.add_argument('-o', '--outname', help='an output filename (.csv) for the calculated adducts')
    parser.add_argument('-m', '--mass_col', default='mass', type=str, help="if the mass column isn't called 'mass'")
    parser.add_argument('-c', '--coverage_cutoff', default='1', type=float, help='if using the default adducts, picks top adducts as percentile or (if >1) number of adducts to calculate')
    parser.add_argument('-r', '--restrict', default=None, type=str, help="name of additional column in the mass file specifying molecular formula (ex formula C8H8)")


    args = parser.parse_args()

    # calculate adduct mz
    if args.adduct_file:
        df = msac.calculate_adduct_mz.calculate_adduct_mz(args.adduct_file, None)
        print("Using supplied adduct file {}. Coverage cutoff not used.".format(args.adduct_file))
    else:
        ADDUCT_FILE = pkg_resources.resource_filename('msac',
                                                      'example_data/adduct_list_full.csv')
        df = msac.calculate_adduct_mz.calculate_adduct_mz(ADDUCT_FILE, args.coverage_cutoff)

    # if formula given, remove adducts if they can't be lost

    # calculate input mass mz for each adduct and add adduct mz
    output = msac.calculate_input_mz.calculate_all_mz(df, args.input_masses,
                                                      args.mass_col, args.restrict)
    print(output.columns)
    print(output.head())

    if args.outname:
        output_name = args.outname
    else:
        inname, ext = path.splitext(args.input_masses)
        output_name = inname + '_adducts.csv'
    output.to_csv(output_name, index=False)
