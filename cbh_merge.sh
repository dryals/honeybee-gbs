#!/bin/bash

# FILENAME: cbh_merge.sh

#SBATCH -A bharpur
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-00:00:00
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

echo "merging plates..."
    #gather completed files
    cd $CLUSTER_SCRATCH/gbs/23CBH
        mkdir -p testmerge4
        cd testmerge4
        #create var with all plate names
        ls ../*/*.json
        #...
        echo -n "" > mergep.txt
        #for i in 1 2 3 4 5 6 7 8 9 10 11 12 17 18 19 20 21 22 23 30 31
        for i in 22 30
        do
            echo -n "../23CBH_${i}/params-23CBH_${i}.txt " >> mergep.txt
        done 
            echo "" >> mergep.txt
        mp=$( cat mergep.txt )
        
    #create merged param file
    ipyrad -m testmerge4 $mp
    
    #edit if needed...
        
echo "launching ipyrad..."
#     ipyrad -p params-merged.txt -s 34567 -c $SLURM_NTASKS -d -f --MPI
    #cd $CLUSTER_SCRATCH/gbs/23CBH/varcalltest
    ipyrad -p params-testmerge4.txt -s 567 -c $SLURM_NTASKS -d -f --MPI
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








    

