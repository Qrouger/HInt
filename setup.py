"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='HInt',
    version='0.1',
    description=(
        'Interaction homolog search tool'
    ),
    author='Quentin Rouger',
    author_email='quentin.rouger@univ-rennes.fr',
    license='GPL-3.0 license',
    url='https://github.com/Qrouger/HInt',
    include_package_data=True,
    packages=find_packages(),
    install_requires=[
        'alphapulldown',
 #       'seaborn',
 #       'urllib3',
        'matplotlib',
#        'scipy',
        'torchdata==0.9.0',
        'pandas',
        'pydantic',
        'packaging',
        'opt_einsum',
        'torch-geometric'
    ],
    entry_points={'console_scripts': ['PPIFold=PPIFold.PPIFold:main',],}
)
