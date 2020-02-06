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
    name='Molecular Spectrometry Adduct Calculator',
    version='0.1.0',
    description='Calculates m/z of potential adducts given compound m/z',
    long_description=readme,
    author='Madison Blumer',
    author_email='madison.blumer@pnnl.gov',
    url='https://github.com/m-blumer/msac',
    license=license,
    packages=pkgs,
    install_requires=required,
    #entry_points={
    #    'console_scripts': ['darkchem = darkchem.cli:main']
    #}
)