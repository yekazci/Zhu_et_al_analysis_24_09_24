#!/bin/bash

# SLURM directives:


#SBATCH -J PANDO_GRN
#SBATCH --export=ALL
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=yusufenes.kazci@mdc-berlin.de
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt


Rscript PANDO_Infer_GRN-SLURM_16_cpu_200G.R

echo message from shell: FINISHED!
