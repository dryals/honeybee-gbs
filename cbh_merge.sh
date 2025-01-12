#!/bin/bash

# FILENAME: cbh_merge.sh

#SBATCH -A bharpur
#SBATCH --ntasks=30
#SBATCH --mem-per-cpu=6G
#SBATCH --time=7-00:00:00
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
   
    #branch, remove samples, and output vcf
    ipyrad -p params-varcalltest.txt -b varcallfinal - 23CBH193_4 23CBH149_2 23CBH323_4 23CBH224_1 \
        23CBH183_5 23CBH347_2 23CBH121_1 23CBH196_6 23CBH335_6 23CBH344_8
   
    ipyrad -p params-varcallfinal.txt -s 7 -c $SLURM_NTASKS -d -f --MPI
   
   
    #trying 6GB by 30cores
        #try adjusting params 11 and 12 to increase usable data (but decrease quality?)
        #try adjusting params 14 for quicker s6
        
        #remove '23CBH193_4', '23CBH149_2', '23CBH323_4', '23CBH224_1', '23CBH183_5', '23CBH347_2', '23CBH121_1', '23CBH196_6', '23CBH335_6', '23CBH344_8'
   
  
  #old
    
    #18 cores * 11Gb seems to work ... 
        #try changing params to make consensus calls (s5) faster?
            #dcrease 19, 
        #try smaller batches in parallel??
        
        # ... try assembling with somehting else after ipyrad filters and demultiplexes ...
    
        #only 5 plates
        #s3 alone may take 4.5 days with only 4 tasks... 
            #try 20 tasks * 6GB
                # ~ 1.4 days 
                #might not be enough mem?
                #not neccesarily scaling with number of cores.. .may need to break bulk
            #also try 64 cores on highmem for 1d .. although this might not complete in time...
            #finally try re-filtering with stricter settings (v2) to reduce data
            #... or several merging steps to get past s3 in reasonable time ...
        #I'm now thinking only merging at s5 ... 
            

echo "-------------"
echo "done"
date








    

