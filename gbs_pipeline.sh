#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --nodes=1 
#SBATCH --ntasks=32
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

#first-time activation
#     conda create -n ipyrad
#     conda activate ipyrad
#     conda install ipyrad -c conda-forge -c bioconda
    
#initalize parameter file 
# cd $CLUSTER_SCRATCH/gbs
# mkdir -p ipyrad
# ipyrad -n test-gbs

#run R to output barcodes
    #...

#copy over parameter file
cp params-test-gbs.txt $CLUSTER_SCRATCH/gbs/ipyrad

#rename fastqs
#     cd $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1
#     mv Bag13_p1_L1_1.fq.gz Bag13_p1_L1_R1_.fastq.gz
#     mv Bag13_p1_L1_2.fq.gz Bag13_p1_L1_R2_.fastq.gz
    
echo "starting ipyrad..."
    cd $CLUSTER_SCRATCH/gbs/ipyrad
    ipyrad -p params-test-gbs.txt -s 1 -c $SLURM_NTASKS -d -f


####
echo "DONE"
date





    

