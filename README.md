# Molecular Spectrometry Adduct Calculator (MSAC)

The Molecular Spectrometry Adduct Calculator (MSAC) calculates the m/z of potential adducts from given compound m/z. This can help in reading a mass spectrometry spectra.
The user provides a .csv with one column containing monoisotopic masses titled 'mass'; other columns will be preserved. Optionally, the user can also supply a more targeted list of adducts to calculate: a .csv with columns 'adduct' and 'charge', where adducts are written without brackets such as M+Na or M-H. 

## Getting Started

These instructions will get you to run the calculator as a python package.

### Prerequisites

Dependencies are pandas, numpy, and molmass.


Using conda: 
```
conda create -n msac-env python pandas numpy
conda activate msac-env
pip install molmass
```

### Installing

In your activated conda environment, install the calculator.

```
# clone/install
git clone https://github.com/pnnl/msac.git
pip install msac/

# direct
pip install git+https://gpnnl/m-blumer/msac
```

## Running the calculator

MSAC is provided with a list of 300+ adducts and their charges. To calculate the m/z for these adducts based on your input masses, pass MSAC a .csv with the compound masses in a column titled 'mass'. Other columns present in the file will be preserved.
```
msac input.csv
```
If your input file's monoisotopic mass column is not called 'mass', you can pass in the title as a string.
```
msac input.csv -m "Monoisotopic Mass"
```
You can also specify an output file name; default is {input_name}_adducts.csv.
```
msac input.csv -o my_output_name.csv
```
If you want to use your own list of adducts, create a csv with a column called 'adduct' and an column of 'charge'. Check example_data/adduct_list.csv for an example.
```
msac input.csv -f my_adduct_list.csv
```

## Authors

* **Madison Blumer** 

## License

This project is licensed under the BSD License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thank you to PurpleBooth for a great [README template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2).
