#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --job-name gbs_pipeline
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out

#Dylan Ryals 06 MAY 2024
#last edited 15 MAY 2024

date

#iPyrad

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
cp params-bag13p1.txt $CLUSTER_SCRATCH/gbs/ipyrad
cp params-bag13p2.txt $CLUSTER_SCRATCH/gbs/ipyrad

# 
# #rename fastqs: _R1_ and _R2_ required in filename!!!
#     cd $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1
#     mv Bag13_p1_L1_1.fq.gz Bag13_p1_L1_R1_.fastq.gz
#     mv Bag13_p1_L1_2.fq.gz Bag13_p1_L1_R2_.fastq.gz
#     cd $CLUSTER_SCRATCH/gbs/bag13/Bag13_p2
#     mv Bag13_p2_L1_1.fq.gz Bag13_p2_L1_R1_.fastq.gz
#     mv Bag13_p2_L1_2.fq.gz Bag13_p2_L1_R2_.fastq.gz


echo "starting ipyrad..."
#     cd $CLUSTER_SCRATCH/gbs/ipyrad
#     #first step: demultiplexing
#     ipyrad -p params-bag13p1.txt -s 1 -c $SLURM_NTASKS -d -f --MPI
#     ipyrad -p params-bag13p2.txt -s 1 -c $SLURM_NTASKS -d -f --MPI
#     
#     ipyrad -m bag13 params-bag13p1.txt params-bag13p2.txt
# 
#     #attempt the rest of the steps
#     ipyrad -p params-bag13.txt -s 234567 -c $SLURM_NTASKS -d -f --MPI
    ipyrad -p params-bag13.txt -s 567 -c $SLURM_NTASKS -d -f --MPI

    
        #with 20 tasks and 10GB per, this takes 24hrs to get to finish step 4 ...
    
####
echo "DONE"
date








    

