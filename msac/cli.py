import argparse
from os import path
import numpy as np
import pandas as pd
import msac


def main():
    p = {}

    parser = argparse.ArgumentParser(description='Calculates the m/z of potential adducts for a given compound m/z')
    parser.add_argument('-v', '--version', action='version', version=msac.__version__, help='print version and exit')

    parser.add_argument('input_masses', help=".csv with 'mass' column")
    parser.add_argument('-f','--adduct_file', help="path to a .csv with an 'adduct' col and a 'charge' col. Defaults to 'example_data/adduct_list.csv' ")
    parser.add_argument('-o', '--outname', help='an output filename (.csv) for the calculated adducts')

    args = parser.parse_args()

    # calculate adduct mz
    if args.adduct_file:
        df = msac.calculate_adduct_mz.calculate_adduct_mz(args.adduct_file) 
    else:
        df = msac.calculate_adduct_mz.calculate_adduct_mz('example_data/adduct_list.csv')

    # calculate input mass mz for each adduct and add adduct mz
    output = msac.calculate_input_mz.calculate_total_mz(df, args.input_masses)

    if args.outname:
        output_name = args.outname
    else:
        inname, ext = path.splitext(args.input_masses)
        output_name = inname + '_adducts.csv'
    output.to_csv(output_name, index=False)

    

    