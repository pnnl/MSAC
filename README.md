# Molecular Spectrometry Adduct Calculator (MSAC)

The Molecular Spectrometry Adduct Calculator (MSAC) calculates the m/z of potential adducts from given compound m/z. This can help in reading a mass spectrometry spectra.
The user provides a .csv with one column containing monoisotopic masses titled 'mass'; other columns will be preserved. Optionally, the user can also supply a more targeted list of adducts to calculate: a .csv with columns 'adduct' and 'charge', where adducts are written without brackets such as M+Na or M-H.

## Getting Started

We recommend you install `msac` into its own conda environment using the instructions below.

### Prerequisites

Recommended installation is through anaconda, so dependencies are compatible with each other.

Dependencies are pandas, numpy, and molmass.

Using conda:

```bash
conda env create -f msac-env.yml
conda activate msac-env
```

A requirements.txt file is also available.

### Installing

In your activated conda environment, install the calculator. The easiest way to install MSAC is via `pip`:

```bash
pip install msac
```

## Running the calculator

MSAC is provided with two adduct lists: a default list of 13 common adducts and an extensive list of adducts found in NIST17.
To calculate the m/z for these adducts based on your input masses, pass MSAC a .csv with the compound masses in a column titled 'mass'. Other columns present in the file will be preserved.

```bash
msac input.csv
```

If your input file's monoisotopic mass column is not called 'mass', you can pass in the title as a string.

```bash
msac input.csv -m "Monoisotopic Mass"
```

You can also specify an output file name; default is {input_name}_adducts.csv.

```bash
msac input.csv -o my_output_name.csv
```

If you want to use your own list of adducts, create a csv with a column called 'adduct' and an column of 'charge'. Check example_data/adduct_list_full.csv for an example.

The default adduct list included in the distribution is used here.

```bash
msac input.csv -f adduct_only_list.csv
```

If you want to include potential neutral losses, add the -n flag.

```bash
msac input.csv -n
```

If you're using the default adduct list or the neutral loss list with -n, you can limit the number of adducts used based on its prevalence in NIST/GNPS/MassBank. Using a value 1.0 or less will give you all adducts that cover more than the specified percent of NIST/GNPS/MassBank. Using an integer value 2+ gives that number of adducts.

```bash
msac input.csv -c 0.75 -o "75th_percentile_adducts.csv"
msac input.csv -c 5 -o "top_5_most_common_adducts.csv"
```

If you have a column of formulas in your input file, you can have MSAC calculate the input masses for you and/or have it output 'NaN' values when a neutral loss cannot be done due to the number of atoms in the molecule.

```bash
# have MSAC calculate the mass using --formula_col
msac input.csv -f "Formula"
# have MSAC restrict losses to those numerically possible using --restrict
msac input.csv -r "Formula"
```

## Authors

* **Madison Blumer**

## License

This project is licensed under the BSD License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thank you to PurpleBooth for a great [README template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2).
