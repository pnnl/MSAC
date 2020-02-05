'''
Author: Madison Blumer
Purpose: Calculate mass for adducts of a given M.
Notes:
    Molmass gets its mass values from NIST
'''
import pandas as pd
import numpy as np 
import re
from molmass import Formula
import argparse

abbrev_to_formula = {'ACN':'CH3CN', 'DMSO':'C2H6OS', 'FA':'CH2O2', 
                    'HAc':'CH3COOH', 'TFA':'C2HF3O2',
                    'IsoProp':'CH3CHOHCH3', 'MeOH':'CH3OH'}

MASS_ELECTRON_DALTON = 0.00054857990924 # Fiehn https://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/

def get_ions(s, atom_dict):
    '''
    identifies the ions present in the adduct set and splits the adduct names
    '''
    ions = re.split('(\+|\-|\.)', s)
    ions = ions[1:]
    ions = split_coeff(ions)
    for ion in ions:
        for token in ion:
            if token not in atom_dict and not token.isdigit():
                atom_dict[token] = 1
    return ions

def split_coeff(ions):
    '''
    splits coefficients off ion names
    '''
    coeff = []
    all_coeff = []
    for ion in ions:
        if ion[0].isdigit():
            coeff = [ion[0], ion[1:]]
        else:
            coeff = [ion]
        all_coeff.append(coeff)
    return all_coeff

def get_adduct_masses(atom_dict, mass_dict, all_atoms):
    # calculate masses for each adduct, before accounting for electrons
    adduct_mass = {}
    for form in all_atoms:
        # all atoms is list of lists of lists
        # all_atoms = list of adducts, adduct = list of seperated adduct (even (incl 0) indices are +/-, odd are list of atoms), 
        #       atom = list containing any multipler plus the ion tag
        mass = 0.0
        adduct = ''
        for ind in range(len(form)):
            # remakes the full adduct tag out of the listed version as it loops
            adduct = adduct + ''.join(form[ind])
            # if the current part of the adduct is an atom
            if ind % 2 != 0:
                if len(form[ind]) == 1:
                    # ['Na'], ['H']
                    m = mass_dict[form[ind][0]]
                else:
                    # ['2','H'], ['3','Na']
                    m = float(form[ind][0]) * mass_dict[form[ind][1]]
                # mass multipler is specified by the previous listed adduct part, either +/-
                if form[ind-1][0] == '+':
                    mult = 1
                elif form[ind-1][0] == '-':
                    mult = -1
                else:
                    mult = 0
                    print("issue", form[ind-1])
                mass = mass + mult*m
        adduct_mass[adduct] = mass
    return adduct_mass

def get_atom_masses(df):
    # obtain masses for all possible atoms in the list. will also change things that are in the abbreviation dict to formula
    atom_dict = {}
    mass_dict = {}
    all_atoms = [get_ions(s, atom_dict) for s in df['adduct']]
    for atom in atom_dict.keys():
        if atom == '+' or atom == '-':
            continue
        if atom in abbrev_to_formula:
            f = Formula(abbrev_to_formula[atom])
        else:
            f = Formula(atom)
        if atom not in mass_dict:
            mass_dict[atom] = f.isotope.mass #monoisotopic mass
    return atom_dict, mass_dict, all_atoms

def get_compound_masses(mass, df):
    # columns: has 'adduct', 'charge', 'input_mass_multiplier', 'm/z'

    df['expected m/z'] = [(mult*mass)/np.absolute(charge) for mult,charge in zip(df['input_mass_multiplier'],df['charge'])]

    df['expected m/z'] = [mass + adduct_mass for mass,adduct_mass in zip(df['expected m/z'], df['m/z'])]

    output_df = pd.concat([df['expected m/z'], df['adduct'], df['m/z'], df['charge']], 
                          axis=1, names=['expected_m/z','adduct','adduct_m/z','charge'])

    return output_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given an input mass and a .csv of adducts with charges, returns a DataFrame of the possible adducts and their m/z values. ")
    parser.add_argument('input_masses', help=".csv with 'mass' column")
    parser.add_argument('-f','--adduct_file', help="path to a .csv with an 'adduct' col and a 'charge' col. Defaults to 'example_input.csv' ")

    args = parser.parse_args()

    if args.adduct_file:
        df = pd.read_csv(args.adduct_file, dtype={'charge':float})
        try:
            df['adduct'][0]
        except:
            print("Please title the adduct name column as 'adduct'. Adducts in form 2M+H, M-H+Na")
        try:
            df['charge'][0]
        except:
            print("Please title the charge column 'charge'. Charges in form 2, -2, 1, -1")  
    else:
        df = pd.read_csv('example_input.csv')
    
    # add multiplier column for calculating electron gain/loss
    df['electron_multiplier'] = [-s for s in df['charge']]
    # multiplier for input mass--if it's a 2M versus and M adduct
    df['input_mass_multiplier'] = [float(s[0]) if s[0].isdigit() else 1 for s in df['adduct']]

    # calculate masses of atoms in the addutcs, then the adducts themselves
    atom_dict, mass_dict, all_atoms = get_atom_masses(df)
    
    adduct_mass = get_adduct_masses(atom_dict, mass_dict, all_atoms)

    # get the calculated mass for the adduct from the dict made above
    df['mass_noelectron'] = [adduct_mass[s[1:]] if s[0] == 'M' else adduct_mass[s[2:]] for s in df['adduct']]

    # get accurate mass that INCLUDES electron lost/gained with charge
    df['mass'] = [x + mult*MASS_ELECTRON_DALTON for x, mult in zip(df['mass_noelectron'],df['electron_multiplier'])]

    df['m/z'] = [mass/np.absolute(charge) for mass, charge in zip(df['mass'], df['charge'])]


    # Calculate adducts + input masses
    d = {adduct:(mult, mass) for adduct, (mult, mass) in zip(df['adduct'], zip(df['input_mass_multiplier'], df['m/z']))}

    input_masses = pd.read_csv(args.input_masses)
    masses_to_calc = input_masses['mass']
    all_masses = []
    for adduct in d.keys():
        m = [d[adduct][1] + d[adduct][0]*mass for mass in masses_to_calc]
        all_masses.append(m)
    output_df = pd.DataFrame(all_masses, index=list(d.keys()), columns=masses_to_calc).transpose()

    #output_df = get_compound_masses(float(args.mass), df)

    output_df.to_csv('masses_recalc.csv', index=True)
    
