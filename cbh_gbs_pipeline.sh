#!/bin/bash

# FILENAME: cbh_gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:00:00
#SBATCH --job-name cbh_gbs
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/cbh.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/cbh.out

#Dylan Ryals 05 DEC 2024
#last edited 

date

#iPyrad
#this takes around 30 hours on 120 cores...

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

P=1 #plate number

echo "setup plate $P ... "




#copy parameter file into scratch directory
    #TODO: create param files with vim or something
    mkdir -p $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
    cp params/params-23CBH_${P}.txt $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}


#rename fastqs: _R1_ and _R2_ required in filename!!!
    #TODO: ensure all dirs and fastq are named the same...
#     cd data/CBH2023
#     if[! dir 23CBH_${P}]
#         mv CBH_${P} 23CBH_${P}
#         cd 23CBH_${P}
#         mv CBH_${P}_1.fq.gz 23CBH_${P}_R1_.fastq.gz
#         mv CBH_${P}_2.fq.gz 23CBH_${P}_R2_.fastq.gz
#     fi
#     
    

echo "starting ipyrad..."
    cd $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
    #first step: demultiplexing
    ipyrad -p params-23CBH_${P}.txt -s 1 -c $SLURM_NTASKS -d -f --MPI
        #8 cores * 10GB work, not minimum
            #time: 
    
#OLD
    
#     ipyrad -p params-bag13p2.txt -s 1 -c $SLURM_NTASKS -d -f --MPI
#     
#     ipyrad -m bag13 params-bag13p1.txt params-bag13p2.txt
# 
#     #attempt the rest of the steps
#     ipyrad -p params-bag13.txt -s 234567 -c $SLURM_NTASKS -d -f --MPI
    
    #drop samples that failed to assemble
#    ipyrad -p params-bag13.txt -b bag13-final - 23-II18w09
    
    #output final results
#    ipyrad -p params-bag13-final.txt -s 7 -c $SLURM_NTASKS -d -f --MPI

    
####
echo "DONE"
date








    

