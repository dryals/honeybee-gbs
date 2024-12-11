#!/bin/bash

# FILENAME: array_cbh_gbs.sh

#SBATCH -A bharpur
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-10:00:00
#SBATCH --job-name array_cbh_gbs
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/dump_%a.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/dump_%a.out

#Dylan Ryals 06 DEC 2024
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

#output logfile to track all jobs at once
log=/home/dryals/ryals/honeybee-gbs/outputs/array.out

#additional setup if first task
     #TODO: somehow verify task 1 has done this before continuing 
        #just manually reset log for now...
        #echo -n "" > $log
# if [  $SLURM_ARRAY_TASK_ID == 1 ]; then 
#     date > $log
# 
#     #sort which plates to process
#         echo -n "" > $CLUSTER_SCRATCH/gbs/23CBH/todo.txt
#         #for i in 1 2 3 4 5 6 7 8 10 11 12 17 18
#         #for i in 19 20 21 22 23 30 31 9
#         #for i in 4 5 6 8 9 10 11 17 18 19 20 21 22 23 30 31
#         for i in 22 30
#         do
#             echo $i >> $CLUSTER_SCRATCH/gbs/23CBH/todo.txt
#         done
# fi
    

    k=1 #number of plates per job
    nstart=$((( $SLURM_ARRAY_TASK_ID - 1 ) * $k + 1))
    nend=$(( $nstart + $k - 1 ))
    #echoing commands to stdout and logfile so progress can be tracked both places
    echo "task $SLURM_ARRAY_TASK_ID processing n = $nstart - $nend" >> $log
        echo "task $SLURM_ARRAY_TASK_ID processing n = $nstart - $nend"
 
#main processing loop
cat $CLUSTER_SCRATCH/gbs/23CBH/todo.txt | sed -n "${nstart},${nend} p" | while read P
do
        
    echo "    starting plate ${P}" >> $log
        echo -n "        " >> $log
        date >> $log
        echo "    starting plate ${P}" 


    #copy parameter file into scratch directory
        #new working dir for this plate
        mkdir -p $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
        #edit param file to use plate name and save to dir
        cd ~/ryals/honeybee-gbs
        param=$( cat params/params-23CBH_PLATE.txt )
        echo "${param//PLATE/"$P"}" > $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}/params-23CBH_${P}.txt
        
    #rename dirs and fastqs
        #_R1_ and _R2_ required in filename!!!
        cd data/CBH2023
        #exit if file not found
        ls *CBH_${P} &> /dev/null || ( echo "dir for ${P} does not exist" ; exit )
        #rename dirs without year code
        if [ -d "CBH_${P}" ]; then
            echo "renaming dir"
            mv CBH_${P} 23CBH_${P}
            cd 23CBH_${P}
            mv CBH_${P}_1.fq.gz 23CBH_${P}_R1_.fastq.gz
            mv CBH_${P}_2.fq.gz 23CBH_${P}_R2_.fastq.gz
            cd ..
        fi
        #rename dirs with year code
        if [ -f "23CBH_${P}/23CBH_${P}_1.fq.gz" ]; then
            echo "renaming fastq"
            cd 23CBH_${P}
            mv 23CBH_${P}_1.fq.gz 23CBH_${P}_R1_.fastq.gz
            mv 23CBH_${P}_2.fq.gz 23CBH_${P}_R2_.fastq.gz
            cd ..
        fi    

    #ipyrad
        cd $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
        #first step: demultiplexing
        ipyrad -p params-23CBH_${P}.txt -s 2 -c $SLURM_NTASKS -d -f --MPI
            #s1 8 cores * 10GB work, not minimum
                #time: 2:55
            #s1-5 8 cores * 10 GB (48 CPU): basically 1 day
                #i'm not sure which steps are better in parallel or together...
                #test if we can get by with less GB and tasks may take fine-tuning per step...
            #testing s12 with 6GB per task
                #this take 32 cores, I can run about 3 at a time...
                #s1-2 takes basically 4.5 hrs
            #try decreasing memory and just running s1 for the rest of plates
                #then merge everything at s2 and go from there with max cores... easiest to manage?

                
    echo "    finished plate ${P}" >> $log
        echo -n "        " >> $log
        date >> $log
 
 done
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
echo "DONE task $SLURM_ARRAY_TASK_ID" >> $log
    echo -n "    " >> $log
    date >> $log
date








    

