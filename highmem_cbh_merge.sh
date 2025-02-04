#!/bin/bash

# FILENAME: highmem_cbh_merge.sh

#SBATCH -A highmem
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=30G
#SBATCH --time=1-00:00:00
#SBATCH --job-name hm_cbh_merge
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/hm_merge.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/hm_merge.out

#Dylan Ryals 09 DEC 2024
#last edited 4 FEB 2025

date

#iPyrad


module load anaconda use.own
conda activate ipyrad

echo "-------------"

# echo "merging plates..."
#     #gather completed files
#     cd $CLUSTER_SCRATCH/gbs/23CBH
#         mkdir -p varcalltest
#         cd varcalltest
#         #create var with all plate names
#         ls ../*/*.json
#         #...
#         echo -n "" > mergep.txt
#         for i in 1 2 3 4 5 6 7 8 9 10 11 12 17 18 19 20 21 22 23 30 31
#         #for i in 22 30
#         do
#             echo -n "../23CBH_${i}/params-23CBH_${i}.txt " >> mergep.txt
#         done 
#             echo "" >> mergep.txt
#         mp=$( cat mergep.txt )
#         
#     #create merged param file
#     ipyrad -m varcalltest $mp
#     
#     #edit if needed...
#         
echo "launching ipyrad..."

    #WARNING: ensure the correct param file is used! edit if needed after merging...

    cd $CLUSTER_SCRATCH/gbs/23CBH/varcalltest
#     ipyrad -p params-varcalltest.txt -s 6 -c $SLURM_NTASKS -d -f --MPI
# 
#     rcall=$( cat ~/ryals/honeybee-gbs/data/toremove.txt | tr '\n' ' ' )
# 
#     #branch, remove samples, and output vcf
#     ipyrad -p params-varcalltest.txt -b varcallfinal2 - $rcall
#     
#     #edit param file ...

    ipyrad -p params-varcallfinal2.txt -s 7 -c $SLURM_NTASKS -d -f --MPI
    
    #trying wiht 32GB * 8 cores, maybe the architecture makes a difference?
        #this works, but takes barely > 1 day
    #trying 16GB * 12 cores ... hopefully this is enough ram and cores 
        #16GB was not enough... rude
    #trying 30GB * 10 cores
        

echo "-------------"
echo "done"
date








    

