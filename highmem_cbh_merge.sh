#!/bin/bash

# FILENAME: highmem_cbh_merge.sh

#SBATCH -A highmem
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=48G
#SBATCH --time=1-00:00:00
#SBATCH --job-name hm_cbh_merge
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/hm_merge.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/hm_merge.out

#Dylan Ryals 09 DEC 2024
#last edited 01 APR 2025

date

#iPyrad

module load anaconda use.own
conda activate ipyrad

echo "-------------"

# echo "merging plates..."
#     #gather completed files
#     cd $CLUSTER_SCRATCH/gbs/23CBH
#         mkdir -p varcall-update
#         cd varcall-update
#         
#         #create var with all plate names
#         ls ../*/*.json
#         #...
#         echo -n "" > mergep.txt
#         for i in {1..31}
#         do
#             echo -n "../23CBH_${i}/params-23CBH_${i}.txt " >> mergep.txt
#         done 
#             echo "" >> mergep.txt
#             
#             
#         mp=$( cat mergep.txt )
#         
#     #create merged param file
#     ipyrad -m varcall-update $mp
#     
#     #edit if needed...
#      
     
echo "launching ipyrad..."

    #WARNING: ensure the correct param file is used! edit if needed after merging...
    

    cd $CLUSTER_SCRATCH/gbs/23CBH/varcall-update
    
#     ipyrad -p params-varcall-update.txt -s 6 -c $SLURM_NTASKS -d -f --MPI
#     
#     #hanging on 'building database' with 6GB * 30 cores? Try highmem?
#     #OOM with 12G, trying 16G
#     #finally succedded with 20G * 14 tasks
#     
#     #manual analysis with cbh_analysis...
#     toremove=$( cat toremove.txt )
#     
#     #adjust line 21 in param file to refelct new 90% locus threshold
#     
#     
#     #branch, remove samples, and output vcf
#     ipyrad -p params-varcall-update.txt -b varcall-update-final - $toremove
#  
#     ipyrad -p params-varcall-update-final.txt -s 7 -c $SLURM_NTASKS -d -f --MPI
#     
# 
# #trying branching 
#     toadd=$( cat ~/ryals/honeybee-gbs/data/split-1.txt )
#     ipyrad -p params-varcall-update.txt -b varcall-update-split1 $toadd
#     
#     
#     toadd=$( cat ~/ryals/honeybee-gbs/data/split-2.txt )
#     ipyrad -p params-varcall-update.txt -b varcall-update-split2 $toadd
#     
#     #edit both param files for new line 21
    
    #ipyrad -p params-varcall-update-split1.txt -s 7 -c $SLURM_NTASKS -d -f --MPI
    
    ipyrad -p params-varcall-update-split2.txt -s 7 -c $SLURM_NTASKS -d -f --MPI


    #split in half, trying 48G * 10 cores ...
        #may try more mem, less cores later... hopefully this is fast
        #runs in about 9hrs
    
    
    #trying s7 64G * 8 tasks on highmem
        #may need to re-run everything with stricter filtering to run this in < 24hrs
        # ... or break in half and merge vcfs later???
        #fail
    #trying 80G * 6 tasks
        #... definitely will need to break into separate tasks ...
    
    
    
    #old notes

    
    #s7 froze at 33% completion, trying with 
    #s7 also errors at 8G * 24 cores: ipyparallel.error.EngineError: Engine b'98b61bdc-b5e555fd93b43e8c2ec98921' died while running task '03985a02-0d72077978c48d97fcac8e00_1845710_144'
    #trying 10G * 20cores ... fail
    #tryin 12G * 16 cores ... fail
    #14G * 15 cores ... fail
    #20G * 11 cores ... fail
    #30G * 8 cores ... fail
    #50G * 4 cores ...
        #try adjusting filters to decrease memory usage on s7... or GATK...
        #try highmem: 64G * 2 cores ???? ...at that point, just use GATK seriously
   
    #s6 works with 6GB by 30cores
        #try adjusting params 11 and 12 to increase usable data (but decrease quality?)
        #try adjusting params 14 for quicker s6
        

echo "-------------"
echo "done"
date







    

