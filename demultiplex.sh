#!/bin/bash

# FILENAME: demultiplex.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --job-name mitotype
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/demult.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/demult.out

#Dylan Ryals 06 MAY 2023
#last edited


date

#module load anaconda

#install package:
#     conda create -n gbs
#     conda activate gbs
#     conda install -c conda-forge -c bioconda ultraplex

#     conda activate gbs


module load anaconda r

#prepare barcodes files
    #run r script

# ultraplex
conda activate gbs

ultraplex -i $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1/Bag13_p1_L1_1.fq.gz \
-b barcodes/P1barcodes.txt \
-d $CLUSTER_SCRATCH/gbs/bag13/samples \
-t $SLURM_NTASKS
    
    
    
# #stacks?
# process_radtags -p $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1 \
# -o $CLUSTER_SCRATCH/gbs/bag13/samples \
# -b barcodes/P1barcodes.txt -r -c -q -i gzfastq


####
echo "DONE"
date





    

