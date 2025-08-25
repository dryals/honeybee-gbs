#!/bin/bash

# FILENAME: 24CBH_array.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=7-00:00:00
#SBATCH --partition cpu
#SBATCH --job-name 24CBH_array
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/dump_%a.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/dump_%a.out

#Dylan Ryals 25 AUG 2025
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
# 
# #output logfile to track all jobs at once
 log=/home/dryals/ryals/honeybee-gbs/outputs/array.out
 
#     echo -n "" > $log
#     date > $log
#     

#additional setup if first task
     #TODO: somehow verify task 1 has done this before continuing 
        #just manually reset log for now...
        #echo -n "" > $log
# if [  $SLURM_ARRAY_TASK_ID == 1 ]; then 
#     date > $log
# 
    #sort which plates to process
        echo -n "" > $CLUSTER_SCRATCH/gbs/24CBH/todo.txt
        #for i in 1 2 3 4 5 6 7 8 10 11 12 17 18
        #for i in 19 20 21 22 23 30 31 9
        #for i in 1 2 3 4 5 6 7 8 9 10 11 12 17 18 19 20 21 22 23 30 31
        #for i in 13 14 15 16 24 25 26 27 28 29 
        for i in 1 2 3 4 5
        do
            echo $i >> $CLUSTER_SCRATCH/gbs/24CBH/todo.txt
        done
# fi
    

    k=2 #number of plates per job
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
        mkdir -p $CLUSTER_SCRATCH/gbs/24CBH/24CBH_${P}
        #edit param file to use plate name and save to dir
        cd ~/ryals/honeybee-gbs
        param=$( cat params/params-24CBH_PLATE.txt )
        echo "${param//PLATE/"$P"}" > $CLUSTER_SCRATCH/gbs/24CBH/24CBH_${P}/params-24CBH_${P}.txt
        
    #rename dirs and fastqs
        #_R1_ and _R2_ required in filename!!!
        cd data/24CBH
        #rename files without underscore
        if [ -d "24CBH${P}" ]; then
            echo "renaming dir"
            mv 24CBH${P} 24CBH_${P}
            cd 24CBH_${P}
            mv 24CBH${P}_L1_1.fq.gz 24CBH_${P}_R1_.fastq.gz
            mv 24CBH${P}_L1_2.fq.gz 24CBH_${P}_R2_.fastq.gz
            cd ..
        fi
        #exit if file not found
        #ls *CBH_${P} &> /dev/null || ( echo "dir for ${P} does not exist" ; exit )
        #rename dirs without year code
#         if [ -d "CBH_${P}" ]; then
#             echo "renaming dir"
#             mv CBH_${P} 24CBH_${P}
#             cd 24CBH_${P}
#             mv CBH_${P}_1.fq.gz 24CBH_${P}_R1_.fastq.gz
#             mv CBH_${P}_2.fq.gz 24CBH_${P}_R2_.fastq.gz
#             cd ..
#         fi
#         #rename dirs with year code
#         if [ -f "24CBH_${P}/24CBH_${P}_1.fq.gz" ]; then
#             echo "renaming fastq"
#             cd 24CBH_${P}
#             mv 24CBH_${P}_1.fq.gz 24CBH_${P}_R1_.fastq.gz
#             mv 24CBH_${P}_2.fq.gz 24CBH_${P}_R2_.fastq.gz
#             cd ..
#         fi    

    #ipyrad
        cd $CLUSTER_SCRATCH/gbs/24CBH/24CBH_${P}
        ipyrad -p params-24CBH_${P}.txt -s 123 -c $SLURM_NTASKS -d -f --MPI
        
        
            #s1 8 cores * 10GB work, not minimum
                #time: 2:55
            #s1-5 8 cores * 10 GB (48 CPU): basically 1 day
                #i'm not sure which steps are better in parallel or together...
                #test if we can get by with less GB and tasks may take fine-tuning per step...
            #testing s1-4 with 6GB per task
                #this take 32 cores, I can run about 3 at a time...
                #s1-2 takes basically 4.5 hrs
            #try decreasing memory and just running s1 for the rest of plates
                #then merge everything at s2 and go from there with max cores... easiest to manage?
                
            #extending parallel processing wtihin plates to s5
                #does this really require 6gb?
                #time for one plate: 
                #consider GATK if this becomes unmanageable 

                    #
                #8 cores * 8Gb and old params: ~ 650 min consesus calling
                #10 cores *6Gb and new params: ~4hrs s5

                
    echo "    finished plate ${P}" >> $log
        echo -n "        " >> $log
        date >> $log
 
 done

    
####
echo "DONE task $SLURM_ARRAY_TASK_ID" >> $log
    echo -n "    " >> $log
    date >> $log
date








    

