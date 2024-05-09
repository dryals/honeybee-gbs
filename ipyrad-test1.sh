#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00
#SBATCH --job-name gbs_pipeline
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/test.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/test.out

#Dylan Ryals 09 MAY 2023
#last edited

date

#iPyrad

#testing with one lane of 96 samples

module load anaconda use.own
conda activate ipyrad

# #first-time activation
#     conda create -n ipyrad
#     conda activate ipyrad
#     conda install ipyrad -c conda-forge -c bioconda
# 
# #run R to output barcodes
#     Rscript --vanilla --silent barcodes.R
# 
# #copy parameter file into scratch directory
cp params-1test.txt $CLUSTER_SCRATCH/gbs/ipyrad

#barcodes
cd barcodes



echo "starting ipyrad..."
    cd $CLUSTER_SCRATCH/gbs/ipyrad

    
    ipyrad -p params-1test.txt -s 1 -c $SLURM_NTASKS -d -f --MPI



####
echo "DONE"
date





    

