import argparse
from os import path


from msac import process, __version__


def main():
    p = {}

    parser = argparse.ArgumentParser(description='Calculates the m/z of potential adducts for a given compound m/z')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    parser.add_argument('input_masses', help=".csv with 'mass' column")
    parser.add_argument('-a', '--adduct_file', help="path to a .csv with an 'adduct' col and a 'charge' col. Defaults to 'example_data/adduct_list.csv' ")
    parser.add_argument('-n', '--neutral_losses_included', action='store_true', help="add this flag to use the neutral loss/adduct combination list")
    parser.add_argument('-o', '--outname', help='an output filename (.csv) for the calculated adducts')
    parser.add_argument('-m', '--mass_col', default='mass', type=str, help="if the mass column isn't called 'mass'")
    parser.add_argument('-f', '--formula_col', default=None, type=str, help="specify formula column if you don't have exact mass in your dataframe but have formula")
    parser.add_argument('-c', '--coverage_cutoff', default=1.0, type=float, help='if using the default adducts, picks top adducts as percentile or (if >1) number of adducts to calculate')
    parser.add_argument('-r', '--restrict', default=None, type=str, help="name of additional column in the mass file specifying molecular formula (ex formula C8H8)")


    args = parser.parse_args()

    process.process_file(args.input_masses, args.mass_col, args.formula_col, args.adduct_file, args.neutral_losses_included, args.outname, args.coverage_cutoff, args.restrict)
