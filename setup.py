"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='HInt-ppi',
    version='0.7.3',
    description=(
        'A tool to find homologous interactions and speed up AlphaFold-based structural modeling.'
    ),
    author='Quentin Rouger',
    author_email='quentin.rouger@univ-rennes.fr',
    license='GPL-3.0 license',
    url='https://github.com/Qrouger/HInt',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'alphapulldown==2.5.0',
        'matplotlib',
        'nvidia-ml-py',
        'ihm',
        'tokamax==0.0.11',
        'jax[cuda]==0.9.1',
        'pandas',
        'pydantic',
        'packaging',
        'opt_einsum',
        'zstandard==0.23.0', 
        'typeguard==2.13.3',
        'triton==3.1.0',
        'jax_triton==0.2.0',
        'rdkit==2025.9.4',
        'gemmi'
    ],
    entry_points={'console_scripts': ['HInt=HInt.HInt:main',],}
)
