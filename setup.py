"""Install script for setuptools."""
from setuptools import find_packages
from setuptools import setup

setup(
    name='HInt-ppi',
    version='0.4.2',
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
        'alphapulldown==2.1.8',
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
        "nvidia-cublas-cu12==12.8.5.5",
        "nvidia-cuda-cccl-cu12==12.9.27",
        "nvidia-cuda-cupti-cu12==12.9.79",
        "nvidia-cuda-nvcc-cu12==12.9.86",
        "nvidia-cuda-nvrtc-cu12==12.9.86",
        "nvidia-cuda-runtime-cu12==12.9.79",
        "nvidia-cudnn-cu12==9.20.0.48",
        "nvidia-cufft-cu12==11.4.1.4",
        "nvidia-cusolver-cu12==11.7.5.82",
        "nvidia-cusparse-cu12==12.5.10.65",
        "nvidia-nccl-cu12==2.29.7",
        "nvidia-nvjitlink-cu12==12.9.86",
        "nvidia-nvshmem-cu12==3.5.21",
    ],
    entry_points={'console_scripts': ['HInt=HInt.HInt:main',],}
)
