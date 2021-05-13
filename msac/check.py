"""
check : Verify whether a molecule can produce listed adducts.

author: @christinehc
"""
# Imports
import re

from collections import Counter


# Functions
def get_adduct_charge(adduct):
    """Find adduct charge given an adduct name.

    Note: works on most adducts, may need tweaking to include more.

    Args:
    ----
        adduct : str
            Name of adduct as a string (e.g. 'M+H')

    Returns
    -------
        int
            Adduct charge (e.g. 1)

    """
    charge = 0
    neg_ions = ["Cl", "Br", "I", "OH", "HCO2", "H2"]
    h2_pattern = re.compile("-H-H(?![0-9|CHNOPSF])")
    h2_match = h2_pattern.search(str(adduct))
    if h2_match:
        adduct = adduct.replace("-H-H", "-H2")

    # only count -H toward negative charge
    for neg in adduct.split("-")[1:]:
        # if adduct is M-
        if not any(neg):
            charge -= 1
        else:
            if neg[0].isdigit() and neg[1:].split("+")[0] == "H":
                charge -= int(neg[0])
            elif neg.split("+")[0] == "H":
                charge -= 1
            elif neg.split("+")[0] in neg_ions:
                charge -= 1

    # count all +el toward positive charge
    for pos in adduct.split("+")[1:]:
        # if adduct is M+
        if not any(pos):
            charge += 1
        else:
            ions = pos[1:].split("-")[0]
            if pos[0].isdigit() and ions in neg_ions:
                charge -= int(pos[0])
            elif pos in neg_ions:
                charge -= 1
            elif pos[0].isdigit() and ions not in neg_ions:
                charge += int(pos[0])
            else:
                charge += 1

    # # account for coefficients
    # if any(adduct.split("M")[0]):
    #     charge = charge / int(adduct.split("M")[0])

    return charge


def formula_to_dict(formula, coefficient=1):
    """Parse molecular formula into dictionary of each atom count.

    e.g. Parses C6H6 into {'C': 6, 'H': 6}

    Notes:
    - Assumes that a coefficient at the beginning of the formula
      is a multiplicative factor (e.g. 2NaCl = 2 Na + 2 Cl atoms).
    - Fractions or decimals are ignored; only whole integers may
      be used.

    Parameters
    ----------
    formula : str
        Chemical formula in string format.
    coefficient : int or float
        Factor by which to divide the total atom counts (default: 1).


    Returns
    -------
    dict
        Dictionary {element: count} from the input formula, where:
            element : str
                Chemical element (e.g. H, C, Na, ...)
            count : int
                Number of occurrences of the element as specified by
                the input chemical formula

    """
    # remove any preceding or ending -/+ signs for multiatom ions
    if len(formula) > 4:
        if (formula[-1] == '-'
            or formula[-1] == '+') and formula[-2].isdigit():
            formula = formula[:-2]
        elif (formula[-2] == '-'
              or formula[-2] == '+') and formula[-1].isdigit():
            formula = formula[:-2]
    formula = formula.strip("+-")

    # check if there is a multiplier
    check_multiplier = re.search(r"(^[-]{0,1}[\d+]{1,1000}[\w]+)", formula)

    # if multiplier exists, save value and remove
    if check_multiplier:
        # count digits in multiplier portion of formula
        digits = 1
        while check_multiplier.group()[:digits].isdigit():
            digits += 1
        multiplier = int(check_multiplier.group()[:digits - 1])
        formula = "-" + check_multiplier.group()[digits - 1:]
    else:
        multiplier = 1

    # split formula into atom, number groups
    regex = r"[\(]{0,1}[A-Z]{1}[a-z]{0,1}[-\+]{0,1}[\d]{0,10}[\]\)]{0,1}[\d+]{0,1000}"
    split = re.findall(regex, formula)

    # find parentheses for molecular units


    # parse string for atom symbols and add to dictionary
    atom_dict = dict()
    parentheses = False
    for i, s in enumerate(split):
        # for (molecule)_subscript, search for and apply multiplier
        if not parentheses:
            unit = 1
        if '(' in s:
            remaining = "".join(split[i:])
            unit = int(re.search(r"(?<=\))[\d]{0,1}", remaining).group())
            parentheses = True
        if ')' in s:
            parentheses = False
            # unit = int(s.split(')')[1])
            s = s.split(')')[0]

        # parse remaining atom, subscript pairs
        if re.search(r"[A-Za-z]+[-\+]{1}[\d]{0,10}", s):
            pattern = re.search(r"[A-Za-z]+[-\+]{1}[\d]{0,10}", s)
        else:
            pattern = re.search(r"[A-Za-z]+", s)
        idx = pattern.end()
        atom = pattern.group()

        # initialize dict entry if none exists
        if atom not in atom_dict.keys():
            atom_dict[atom] = 0

        # formula subscript
        if s[idx:]:
            atom_dict[atom] += int(s[idx:]) * unit * multiplier / coefficient
        else:
            atom_dict[atom] += 1 * unit * multiplier / coefficient
    return atom_dict


def adduct_in_parent(adduct, parent):
    """Check if an adduct can form from a parent structure.

    Because adducts form based on losses in the parent molecule, the
    function compares the adduct formula against the parent ion
    formula and returns False if the parent ion does not contain the
    requisite atoms.

    Args:
    -----
        adduct : str
            Name of adduct as a string (e.g. "M+H")
        parent : str
            Chemical formula of parent compound (e.g. "C6H6")

    Returns
    -------
        bool
            Returns True for possible or False for impossible adduct.

    """
    # account for coefficient (e.g. 2M-2H only needs to lose -H)
    try:
        coefficient = int(re.search(r"[\d]+(?=M\+|M\-)", adduct).group())
    except AttributeError:  # if no coefficient
        coefficient = 1.0

    # split adduct into additions and losses
    ions = re.findall(r"[+-][\w]+", adduct)

    # parse only losses and expand out all atom symbols for adducts
    losses = [ion for ion in ions if ion[0] == "-"]
    atoms = [formula_to_dict(f, coefficient=coefficient) for f in list(losses)]
    adduct_atoms = sum((Counter(dict(x)) for x in atoms), Counter())

    # count number of atoms that appear in parent molecule
    parent_atoms = formula_to_dict(parent)

    # check that all adduct atoms are contained within parent formula
    for atom in list(adduct_atoms):
        if (atom not in parent_atoms) or (adduct_atoms[atom]
                                          > parent_atoms[atom]):
            return False
    return True
