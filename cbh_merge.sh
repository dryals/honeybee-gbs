#!/bin/bash

# FILENAME: cbh_merge.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-10:00:00
#SBATCH --job-name array_cbh_gbs
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/merge.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/merge.out

#Dylan Ryals 09 DEC 2024
#last edited 

date

#iPyrad


module load anaconda use.own
conda activate ipyrad

#echo
echo "merging plates..."

    cd $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
    #first step: demultiplexing
    ipyrad -p params-23CBH_${P}.txt -s 12 -c $SLURM_NTASKS -d -f --MPI
        #s1 8 cores * 10GB work, not minimum
            #time: 2:55
        #s1-5 8 cores * 10 GB (48 CPU): basically 1 day
            #i'm not sure which steps are better in parallel or together...
            #test if we can get by with less GB and tasks may take fine-tuning per step...
        #testing s12 with 6GB per task
            #this take 32 cores, I can run about 3 at a time...
            #s1-2 takes basically 4.5 hrs
            #failures in parallel mode will be real hard to track...
            #consider just the first 2 steps (or something), then the rest merged...
        #find min memorgy req for first step
            #run max in parallel, run multiple samples in sequence
            #merge everything to run from step2 onwards more efficiently (i think?)
        #try dierting output to separate logfiles?
                

    #TODO: attempt to merge, may need to re-run some with -f to get all on the same page ... 
 
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
echo "DONE set $SLURM_ARRAY_TASK_ID" >> $log
date








    

