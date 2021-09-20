"""
check : Calculate the combined input+ adduct mass.
author: @m-blumer
"""
# Imports
import pandas as pd
import numpy as np
from molmass import Formula
from msac import check

# Functions
def calculate_mass_from_formula(input_df, formula_col):
    input_df['mass'] = [Formula(mol).isotope.mass for mol in input_df[formula_col]]
    return input_df

def calculate_total_mz(adduct_tuple, mass):
    """Calculate adduct mz using multiplier, charge, and mass
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
    input_mz_multiplier = adduct_tuple[0]
    adduct_charge = adduct_tuple[1]
    adduct_mz = adduct_tuple[2]
    #  calculating as adduct mz + ((input multiplier*input mass)/charge)
    total_mz = adduct_mz + input_mz_multiplier*mass/np.absolute(adduct_charge)
    return total_mz


def calculate_all_mz(df, input_masses, mass_col, formula_col):
    """Calculate adduct mz using multiplier, charge, and mass
    Parameters
    ----------
    df : DataFrame
        Contains an entry representing each adduct of interest, 
        including adduct name, charge, m/z, and input mass multipler
    mass_file : str
        File containing masses of the input molecules.
    mass_col : str
        Name of mass column in mass_file.
    formula_col : str
        Name of formula column; None if formula not given
    Returns
    -------
    DataFrame
        Returns a table of caluclated masses across all adducts for each input molecule.
    """
    input_masses = input_masses.copy() # to make sure we aren't modifying the original dataframe
    #  create a lookup table for the adduct information:
    #  name, input mass mulitpler, charge, m/z
    d = {adduct: [mult, charge, mass] for adduct, (mult, (charge, mass))
         in zip(df['adduct'], zip(df['input_mass_multiplier'],
                zip(df['charge'], df['m/z'])))}

    if formula_col is not None:
        input_masses['parent_atoms'] = [check.formula_to_dict(formula) for formula in input_masses[formula_col]]

    masses_to_calc = input_masses[mass_col]
    original_cols = input_masses.columns
    all_masses = []

    for adduct in d.keys():
        if formula_col is not None:
            input_masses['parent_atoms'] = [check.formula_to_dict(formula) for formula in input_masses[formula_col]]
            input_masses[adduct] = [calculate_total_mz(d[adduct], mass)
                                    if check.adduct_in_parent(adduct, parent_atoms)
                                    else np.nan
                                    for mass, parent_atoms in zip(input_masses[mass_col], input_masses['parent_atoms'])]
        else:
            input_masses[adduct] = [calculate_total_mz(d[adduct], mass)
                                    for mass in input_masses[mass_col]]
    input_masses = pd.melt(input_masses, id_vars=original_cols, var_name='adduct', value_name='adduct mass')

    return input_masses
 