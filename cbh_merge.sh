#!/bin/bash

# FILENAME: cbh_merge.sh

#SBATCH -A bharpur
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-10:00:00
#SBATCH --job-name cbh_merge
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/merge.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/merge.out

#Dylan Ryals 09 DEC 2024
#last edited 

date

#iPyrad


module load anaconda use.own
conda activate ipyrad

echo "-------------"

# echo "merging plates..."
#     #gather completed files
#     cd $CLUSTER_SCRATCH/gbs/23CBH
#         mkdir -p merged
#         cd merged
#         #create var with all plate names
#         ls ../*/*.json
#         #...
#         echo -n "" > mergep.txt
#         for i in 1 2 3 7 12
#         do
#             echo -n "../23CBH_${i}/params-23CBH_${i}.txt " >> mergep.txt
#         done 
#             echo "" >> mergep.txt
#         mp=$( cat mergep.txt )
#         
#     #create merged param file
#     ipyrad -m merged $mp
#         #edit if needed...
        
echo "launching ipyrad..."
#     ipyrad -p params-merged.txt -s 34567 -c $SLURM_NTASKS -d -f --MPI
    #try 18 cores qtih 10GB ea ...
    
    cd $CLUSTER_SCRATCH/gbs/23CBH/testmerge
    ipyrad -p params-testmerge.txt -s 34567 -c $SLURM_NTASKS -d -f --MPI
        #s3 alone may take 7 days.. I'm not sure ...
            #try more cores, else may require several merging steps ...
            
    
    
    # 
    #     #attempt the rest of the steps
    #     ipyrad -p params-bag13.txt -s 234567 -c $SLURM_NTASKS -d -f --MPI
        
        #drop samples that failed to assemble
    #    ipyrad -p params-bag13.txt -b bag13-final - 23-II18w09
        
        #output final results
    #    ipyrad -p params-bag13-final.txt -s 7 -c $SLURM_NTASKS -d -f --MPI

echo "-------------"
echo "done"
date








    

