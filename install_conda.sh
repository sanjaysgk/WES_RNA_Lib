#!/bin/bash
# ===============================================
# Automatic Conda installation and environment setup
# For WES & RNA-seq Snakemake workflow
# ===============================================

set -e  # exit on any error

# Where to install Miniconda
CONDA_DIR="$PWD/miniconda3"

# Check if conda already exists
if [ -d "$CONDA_DIR" ]; then
    echo "Conda already installed at $CONDA_DIR"
else
    echo "Downloading Miniconda..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh

    echo "Installing Miniconda silently..."
    bash Miniconda3.sh -b -p "$CONDA_DIR"

    echo "Removing installer..."
    rm -f Miniconda3.sh

    echo "Initializing conda for current shell..."
    source "$CONDA_DIR/etc/profile.d/conda.sh"
fi

# Ensure conda is available
source "$CONDA_DIR/etc/profile.d/conda.sh"

# List of environment YAMLs
ENV_YAMLS=("workflow/envs/get_genome.yaml" 
           "workflow/envs/simulate_reads.yaml" 
           "workflow/envs/validate_genome.yaml"
           "workflow/envs/WES_data.yml"
           "workflow/envs/RNA_data.yml")

# Create environments
for env_file in "${ENV_YAMLS[@]}"; do
    env_name=$(basename "$env_file" .yaml)
    echo "Creating environment: $env_name"
    
    if conda env list | grep -q "$env_name"; then
        echo "Environment $env_name already exists, skipping..."
    else
        conda env create -f "$env_file"
    fi
done

echo "All environments are ready."
echo "Activate an environment using: conda activate <env_name>"
#     $ conda deactivate
