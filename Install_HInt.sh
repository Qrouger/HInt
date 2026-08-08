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
    python=3.12 \
    pdbfixer=1.10 \
    mafft \
    kalign2 \
    hhsuite \
    hmmer \
    mmseqs2 \
    git \
    setuptools=81

echo "=== Installing HInt ==="

$CONDA_BIN run -n $ENV_NAME pip install -U fast-hint-ppi


$CONDA_BIN run -n $ENV_NAME pip install  --no-deps \
"colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"

echo "=== Installing AlphaFold3 ==="
if [ -d "alphafold3" ]; then
    rm -rf alphafold3
fi

$CONDA_BIN run -n $ENV_NAME git clone https://github.com/KosinskiLab/alphafold3

cd alphafold3

$CONDA_BIN run -n $ENV_NAME git checkout 86b9ea3feacc8934e6e2a581c49eb4c37a2a3d20

$CONDA_BIN run -n $ENV_NAME pip install .

$CONDA_BIN run -n $ENV_NAME build_data || echo "WARNING: build_data failed"

cd ..

echo "=== DeepLocPro ==="
if [ ! -d "deeplocpro" ]; then
   $CONDA_BIN run -n $ENV_NAME git clone https://github.com/Jaimomar99/deeplocpro
fi
cd deeplocpro

$CONDA_BIN run -n $ENV_NAME pip install -q . torch==2.6.0
$CONDA_BIN run -n $ENV_NAME pip install -q triton==3.1.0

cd ..

$CONDA_BIN run -n $ENV_NAME pip install nvidia-cudnn-cu12==9.25.0.15


echo "=== Installation completed successfully ==="

