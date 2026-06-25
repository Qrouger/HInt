#!/bin/bash

set -e
set -o pipefail


# === CONFIG ===
ENV_NAME="HInt"
CONDA_BIN="${CONDA_BIN:-conda}"
ENV_PATH="$($CONDA_BIN info --base)/envs/$ENV_NAME"

echo "=== Initializing conda ==="
if [[ "$CONDA_BIN" == *micromamba* ]] || [[ "$CONDA_BIN" == *mamba* ]]; then
    eval "$($CONDA_BIN shell hook --shell=bash)"
else
    source "$($CONDA_BIN info --base)/etc/profile.d/conda.sh"
fi


if $CONDA_BIN env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "=== Removing existing conda environment ==="
    $CONDA_BIN remove -n "$ENV_NAME" --all -y

elif [ -d "$ENV_PATH" ]; then
    echo "=== Removing leftover environment directory ==="
    rm -rf "$ENV_PATH"

fi

echo "=== Creating environment ==="
$CONDA_BIN create -n $ENV_NAME -y \
    -c conda-forge -c bioconda \
    python=3.11 \
    pdbfixer=1.9 \
    mafft \
    kalign2 \
    hhsuite \
    hmmer \
    mmseqs2 \
    git \
    setuptools=81

echo "=== Installing HInt ==="

$CONDA_BIN run -n $ENV_NAME pip install -U hint-ppi

#$CONDA_BIN run -n $ENV_NAME pip uninstall -y colabfold || true

$CONDA_BIN run -n $ENV_NAME pip install  --no-deps \
"colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"

$CONDA_BIN run -n $ENV_NAME pip install numpy==1.26.4


echo "=== Installing AlphaFold3 ==="
if [ -d "alphafold3" ]; then
    rm -rf alphafold3
fi

$CONDA_BIN run -n $ENV_NAME git clone https://github.com/KosinskiLab/alphafold3

cd alphafold3

$CONDA_BIN run -n $ENV_NAME git checkout 6ad1a65994c2111d291a386cdc048d8c9bfae4af


$CONDA_BIN run -n $ENV_NAME pip install . --no-deps

$CONDA_BIN run -n $ENV_NAME build_data || echo "WARNING: build_data failed"

$CONDA_BIN run -n $ENV_NAME $CONDA_BIN install nvidia/label/cuda-12.4.1::cuda -c nvidia/label/cuda-12.4.1 -y

cd ..

$CONDA_BIN install -y -c nvidia/label/cuda-12.4.1 cuda

echo "=== DeepLocPro ==="
if [ ! -d "deeplocpro" ]; then
   $CONDA_BIN run -n $ENV_NAME git clone https://github.com/Jaimomar99/deeplocpro
fi
cd deeplocpro

$CONDA_BIN run -n $ENV_NAME pip install -q . torch==2.6.0
$CONDA_BIN run -n $ENV_NAME pip install -q triton==3.1.0

cd ..

echo "=== Installation completed successfully ==="

