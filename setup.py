"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='HInt-ppi',
    version='0.6.2',
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
        'alphapulldown',
        'matplotlib',
        'nvidia-ml-py',
        'ihm',
        'scipy==1.16.0',
        'jax[cuda12]==0.5.3',
        'numpy==1.26.4',
        'pandas',
        'pydantic',
        'packaging',
        'opt_einsum',
        'rdkit==2024.3.5', 
        'zstandard==0.23.0', 
        'jaxtyping==0.2.34', 
        'typeguard==2.13.3',
        'jax_triton==0.2.0',
        'triton==3.1.0',
        'torch>=1.6',
        'gemmi'
    ],
    entry_points={'console_scripts': ['HInt=HInt.HInt:main',],}
)
