"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='HInt-ppi',
    version='0.8.3',
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
        'matplotlib',
        'nvidia-ml-py',
        'ihm',
        'tokamax==0.0.12',
        'flax==0.12.9',
        'qwix==0.1.8',
        "jax[cuda12]==0.10.2; sys_platform != 'darwin'",
        "jax-mps==0.10.9; sys_platform == 'darwin' and platform_machine == 'arm64'",
        'pandas',
        'pydantic',
        'packaging',
        'opt_einsum',
        'zstandard==0.23.0', 
        'typeguard==2.13.3',
        'triton==3.1.0',
        'jax_triton==0.2.0',
        'rdkit==2025.9.4',
        'gemmi',
        'scipy==1.18.1',
        'alphafold-colabfold==2.3.20'
    ],
    entry_points={'console_scripts': ['HInt=HInt.HInt:main',],}
)
