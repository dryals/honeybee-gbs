#!/bin/bash

# FILENAME: gbs_pipeline.sh

#SBATCH -A bharpur
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-00:00:00
#SBATCH --job-name gbs_pipeline
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pipeline.out

#Dylan Ryals 06 MAY 2023
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
cp params-test-gbs.txt $CLUSTER_SCRATCH/gbs/ipyrad
# 
# #rename fastqs: _R1_ and _R2_ required in filename!!!
#     cd $CLUSTER_SCRATCH/gbs/bag13/Bag13_p1
#     mv Bag13_p1_L1_1.fq.gz Bag13_p1_L1_R1_.fastq.gz
#     mv Bag13_p1_L1_2.fq.gz Bag13_p1_L1_R2_.fastq.gz
#  

echo "starting ipyrad..."
    cd $CLUSTER_SCRATCH/gbs/ipyrad
    #first step: demultiplexing
    #ipyrad -p params-test-gbs.txt -s 1 -c $SLURM_NTASKS -d -f
        #this takes around 1.5hr with 32 cores
        
    #ipyrad -p params-test-gbs.txt -s 2 -c $SLURM_NTASKS -d -f
        #this is pretty quick
        #read filter could be more stringent
    
    #ipyrad -p params-test-gbs.txt -s 3 -c $SLURM_NTASKS -d -f --MPI
        #s3 runs with 6GB mem per task
        #takes ??hrs on 24 cores
            #"building clusters" hangs at 0% for a while then jumps

    #s4 takes 7 min
            
    ipyrad -p params-test-gbs.txt -s 567 -c $SLURM_NTASKS -d -f --MPI
    
    #some steps are memory limited (3) while others are not (1,2,4...?)
        #may need a better way of managing jobs to optimize ... 
        

    #split to remove dead samp  
    ipyrad -p params-data1 -b test-branch - 23-II18w09


####
echo "DONE"
date








    

