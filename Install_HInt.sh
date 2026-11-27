#!/bin/bash

set -e
set -o pipefail


# === CONFIG ===
ENV_NAME="HInt"
CONDA_BIN="${CONDA_BIN:-conda}"
CONDA_DIR="${CONDA_DIR:-/lustre/fsn1/worksf/projects/rech/conda}"   # envs + caches, need to be modify
ENV_PATH="$CONDA_DIR/envs/$ENV_NAME"                                             # HInt environment only

echo "Conda : $CONDA_BIN"
echo "Environment : $ENV_PATH"

echo "=== Initializing conda ==="
if [[ "$CONDA_BIN" == *micromamba* ]] || [[ "$CONDA_BIN" == *mamba* ]]; then
    eval "$($CONDA_BIN shell hook --shell=bash)"
else
    source "$($CONDA_BIN info --base)/etc/profile.d/conda.sh"
fi

# Environments and package cache in CONDA_DIR (saved in ~/.condarc, also used for pi_score)
$CONDA_BIN config --add envs_dirs "$CONDA_DIR/envs"
$CONDA_BIN config --add pkgs_dirs "$CONDA_DIR/pkgs"
export PIP_CACHE_DIR="$CONDA_DIR/.cache/pip"


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
    libstdcxx-ng=13 \
    setuptools=81

echo "=== Installing HInt ==="

$CONDA_BIN run -n $ENV_NAME pip install "hint-ppi @ git+https://github.com/Qrouger/HInt.git@HInt-slurm"

$CONDA_BIN run -n $ENV_NAME pip install -v "alphapulldown @ git+https://github.com/Qrouger/AlphaPulldown.git"

$CONDA_BIN run -n $ENV_NAME pip install  --no-deps \
"colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"

echo "=== Installing AlphaFold3 ==="

AF3_DIR="${WORK:?WORK is not set}/alphafold3"

if [ -d "$AF3_DIR" ]; then
    rm -rf "$AF3_DIR"
fi

$CONDA_BIN run -n $ENV_NAME git clone https://github.com/KosinskiLab/alphafold3 "$AF3_DIR"

cd "$AF3_DIR"

$CONDA_BIN run -n $ENV_NAME git checkout 98785987139da9781ba271b159c93efaaf23a06f

$CONDA_BIN run -n $ENV_NAME pip install . --no-deps

$CONDA_BIN run -n $ENV_NAME build_data || echo "WARNING: build_data failed"

cd ..

echo "=== DeepLocPro ==="

DLP_DIR="${WORK:?WORK is not set}/deeplocpro"

if [ ! -d "$DLP_DIR" ]; then
   $CONDA_BIN run -n $ENV_NAME git clone https://github.com/Jaimomar99/deeplocpro "$DLP_DIR"
fi
cd "$DLP_DIR"

$CONDA_BIN run -n $ENV_NAME pip install -q . torch==2.6.0
$CONDA_BIN run -n $ENV_NAME pip install -q triton==3.1.0
$CONDA_BIN run -n $ENV_NAME pip install -q nvidia-cudnn-cu12==9.26.0.51


cd ..

echo "=== Installation completed successfully ==="
