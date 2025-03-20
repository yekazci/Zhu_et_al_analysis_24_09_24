#!/bin/bash

# Submit this script inside bash_scripts folder.

# SLURM directives:

#SBATCH -J pando_eGRN
#SBATCH --export=ALL
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=yusufenes.kazci@mdc-berlin.de
#SBATCH --output=slurm_output_%j.txt
#SBATCH --error=slurm_error_%j.txt
#SBATCH --chdir=./                # Set the working directory to the current, which is default.


Rscript ../r_scripts/PANDO_Infer_32_cpu_CHUNKS.R

echo "completed."
