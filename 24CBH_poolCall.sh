#!/bin/bash

# FILENAME: 24CBH_poolCall.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=7-00:00:00
#SBATCH --partition cpu
#SBATCH --job-name poolCall
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pool.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pool.out

#Dylan Ryals 25 AUG 2025
#last edited 

#determining sequence depth and ref allele count for pooled samples

#testing on 2023 data

date
#####

module load biocontainers samtools


#grab some sites

cd /scratch/negishi/dryals/gbs/24CBH/analysis

samtools mpileup -b test.samps \
    -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
    -C 50 -q 20 -Q 20 -l test.sites \
    -o pileuptest.out


#use R script to translate into beethoven format

    
######
echo "DONE"
date








    

