#!/bin/bash

# FILENAME: demultiplex.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH --job-name mitotype
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/demult.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/demult.out

#Dylan Ryals 06 MAY 2023
#last edited


date

module load anaconda

####

#install package:
#     conda create -n gbs
#     conda activate gbs
#     conda install -c conda-forge -c bioconda ultraplex

conda activate gbs












####
echo "DONE"
date





    

