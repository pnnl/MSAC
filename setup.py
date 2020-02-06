from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()
    required = None

pkgs = find_packages(exclude=('examples', 'docs', 'resources'))

setup(
    name='msac',
    version='0.1.0',
    description='Molecular Spectrometry Adduct Calculator: Calculates m/z of potential adducts given compound m/z',
    long_description=readme,
    author='Madison Blumer',
    author_email='madison.blumer@pnnl.gov',
    url='https://github.com/m-blumer/msac',
    license=license,
    packages=pkgs,
    package_data={'': ['example_data/adduct_list.csv']},
    include_package_data=True,
    install_requires=required,
    entry_points={
        'console_scripts': ['msac = msac.cli:main']
    }
)