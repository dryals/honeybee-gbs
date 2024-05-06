#!/bin/bash

# FILENAME: demultiplex.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --job-name demultiplex
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/demult.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/demult.out

#Dylan Ryals 06 MAY 2023
#last edited


date

#stacks?

module load biocontainers stacks

#run r script ...
    #TO DO


process_shortreads \
-1 $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1/Bag13_p1_L1_1.fq.gz \
-2 $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1/Bag13_p1_L1_2.fq.gz \
-b barcodes/P1barcodes.txt \
-o $CLUSTER_SCRATCH/gbs/bag13/samples \
-inline_inline \
-i gzfastq \
-y gzfastq \
-q -c -r

process_shortreads \
-1 $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1/Bag13_p1_L1_1.fq.gz \
-2 $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1/Bag13_p1_L1_2.fq.gz \
-b barcodes/P1barcodes.txt \
-o $CLUSTER_SCRATCH/gbs/bag13/samples \
-inline_inline \
-i gzfastq \
-y gzfastq \
-q -c -r


####
echo "DONE"
date





    

