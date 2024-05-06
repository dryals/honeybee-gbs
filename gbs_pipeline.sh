#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --job-name gbs_pipeline
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out

#Dylan Ryals 06 MAY 2023
#last edited


date

#iPyrad

module load anaconda
conda activate ipyrad

#first-time activation
#     conda create -n ipyrad
#     conda activate ipyrad
#     conda install ipyrad -c conda-forge -c bioconda
    
#initalize parameter file 
# cd $CLUSTER_SCRATCH/gbs
# mkdir -p ipyrad
# ipyrad -n test-gbs

#copy over parameter file
cp params-test-gbs.txt $CLUSTER_SCRATCH/gbs/ipyrad
cd $CLUSTER_SCRATCH/gbs/ipyrad

ipyrad -p params-test-gbs.txt -s 1


####
echo "DONE"
date





    

