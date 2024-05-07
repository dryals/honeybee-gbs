#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --job-name gbs_pipeline
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out

#Dylan Ryals 06 MAY 2023
#last edited


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
# cp params-test-gbs.txt $CLUSTER_SCRATCH/gbs/ipyrad
# 
# #rename fastqs: _R1_ and _R2_ required in filename
#     cd $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1
#     mv Bag13_p1_L1_1.fq.gz Bag13_p1_L1_R1_.fastq.gz
#     mv Bag13_p1_L1_2.fq.gz Bag13_p1_L1_R2_.fastq.gz
#  

echo "starting ipyrad..."
    cd $CLUSTER_SCRATCH/gbs/ipyrad
    #first step: demultiplexing
    #ipyrad -p params-test-gbs.txt -s 1 -c $SLURM_NTASKS -d -f
        #this takes around 1.5hr for one lane of 96 samples with 32 cores
        
    #ipyrad -p params-test-gbs.txt -s 2 -c $SLURM_NTASKS -d -f
    
    ipyrad -p params-test-gbs.txt -s 3 -c $SLURM_NTASKS -d -f
        #s3 requires more than 32 cores ...
    


####
echo "DONE"
date





    

