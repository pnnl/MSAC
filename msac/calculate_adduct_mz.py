"""
check : Calculate the combined input+ adduct mass.
author: @m-blumer
"""

# Imports
import pandas as pd
import numpy as np
import re
import argparse
from os import path
from molmass import Formula

abbrev_to_formula = {'ACN': 'CH3CN', 'DMSO': 'C2H6OS', 'FA': 'CH2O2',
                     'HAc': 'CH3COOH', 'TFA': 'C2HF3O2',
                     'IsoProp': 'CH3CHOHCH3', 'MeOH': 'CH3OH'}

#  from Fiehn Lab
#  https://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
MASS_ELECTRON_DALTON = 0.00054857990924

# Functions
def get_ions(s, atom_dict):
    """Identifies the ions present in the adduct set and splits the adduct names
    Parameters
    ----------
    adduct_tuple : tuple of floats
        input mass multiplier, adduct charge, and adduct mz 
    mass : int or float
        Total mass of input molecule
    Returns
    -------
    float
        Returns the total mz for the adduct/molecule pair.
    """
    ions = re.split(r'(\+|\-|\.)', s)
    ions = ions[1:]
    ions = split_coeff(ions)
    for ion in ions:
        for token in ion:
            if token not in atom_dict and not token.isdigit():
                atom_dict[token] = 1
    return ions


def split_coeff(ions):
    '''
    split coefficients off ion names--can't be processed by molmass
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
    '''
    calculate masses for each adduct, before accounting for electrons
    '''
    adduct_mass = {}
    for form in all_atoms:
        #  all atoms is list of lists of lists
        #  all_atoms = list of adducts,
        #       adduct = list of seperated adduct
        #           (even (incl 0) indices are +/-, odd are list of atoms),
        #       atom = list containing any multipler plus the ion tag
        mass = 0.0
        adduct = ''
        for ind in range(len(form)):
            #  remakes the full adduct tag out of the listed version
            adduct = adduct + ''.join(form[ind])
            #  if the current part of the adduct is an atom
            if ind % 2 != 0:
                if len(form[ind]) == 1:
                    # ['Na'], ['H']
                    m = mass_dict[form[ind][0]]
                else:
                    # ['2','H'], ['3','Na']
                    m = float(form[ind][0]) * mass_dict[form[ind][1]]
                #  mass multipler is specified
                #  by whether the adduct is added or subtracted, either +/-
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
    '''
    obtain masses for all possible atoms in the list.
    will also change things that are in the abbreviation dict to formula
    '''
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
            #  use monoisotopic mass
            mass_dict[atom] = f.isotope.mass
    return atom_dict, mass_dict, all_atoms


def limit_by_percent_coverage(df, cutoff):
    # calcaulate cumulative sum from percent coverage
    cumsumdf = df.percent_coverage.cumsum().to_frame()
    cumsumdf.reset_index(inplace=True)
    cumsumdf.columns = ['adduct','cumsum']
    cumsumdf['index'] = [str(num) for num in cumsumdf.index]
    totaldf = cumsumdf.merge(percoveragedf, on='adduct')
    if cutoff <= 1.0:
        cutoffdf = totaldf[totaldf['cumsum'] <= cutoff] # things at or above the requested quantile
    else:
        cutoffdf = totaldf[totaldf.index < cutoff]
    if cutoffdf.empty:
        cutoffdf = totaldf.sort_values('cumsum').head(1)
    return cutoffdf


def calculate_adduct_mz(fname, cutoff):
    df = pd.read_csv(fname, dtype={'charge': float})
    try:
        df['adduct'][0]
    except KeyError:
        print("Please title the adduct name column as 'adduct'.\
               Adducts in form 2M+H, M-H+Na")
    try:
        df['charge'][0]
    except KeyError:
        print("Please title the charge column 'charge'.\
               Charges in form 2, -2, 1, -1")

    # strip off leading/trainling spaces, if present
    df['adduct'] = [adduct.lstrip('[').split(']')[0] for adduct in df.adduct]
    
    # subset data frame if using default & cutoff value specified
    if cutoff:
        df = limit_by_percent_coverage(df, cutoff)

    # add multiplier column for calculating electron gain/loss
    df['electron_multiplier'] = [-s for s in df['charge']]
    # multiplier for input mass--if it's a 2M versus and M adduct
    df['input_mass_multiplier'] = [float(s[0]) if s[0].isdigit() else 1
                                   for s in df['adduct']]

    # calculate masses of atoms in the addutcs, then the adducts themselves
    atom_dict, mass_dict, all_atoms = get_atom_masses(df)

    adduct_mass = get_adduct_masses(atom_dict, mass_dict, all_atoms)

    # get the calculated mass for the adduct from the dict made above
    df['mass_noelectron'] = [adduct_mass[s[1:]]
                             if s[0] == 'M' else adduct_mass[s[2:]]
                             for s in df['adduct']]

    # get accurate mass that INCLUDES electron lost/gained with charge
    df['mass'] = [x + mult*MASS_ELECTRON_DALTON for x, mult
                  in zip(df['mass_noelectron'], df['electron_multiplier'])]

    df['m/z'] = [mass/np.absolute(charge) for mass, charge
                 in zip(df['mass'], df['charge'])]

    return df
