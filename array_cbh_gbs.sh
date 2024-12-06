#!/bin/bash

# FILENAME: array_cbh_gbs.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1-10:00:00
#SBATCH --job-name array_cbh_gbs
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/dump.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/dump.out

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

log=/home/dryals/ryals/honeybee-gbs/outputs/array.out
echo "here we go" > $log

#sort which plates to process

    for i in 2 4 5 6 7 8 10 11 12 17 18
    do
        echo $i >> todo.txt
    done
    
    SLURM_ARRAY_TASK_ID=3

    k=4 #number of mitotypes per job
    nstart=$((( $SLURM_ARRAY_TASK_ID - 1 ) * $k + 1))
    nend=$(( $nstart + $k - 1 ))
    echo "$nstart $nend"
 
#main processing loop
cat todo.txt | sed -n "${nstart},${nend} p" | while read P
do
        
    echo "starting plate ${P}" >> $log
        echo -n "    " >> $log
        date >> $log



    #copy parameter file into scratch directory
        #TODO: create param files with vim or something
        #new working dir for this plate
        mkdir -p $CLUSTER_SCRATCH/gbs/23CBH/23CBH_${P}
        #edit param file to use plate name and save to dir
        param=$( cat params/params-23CBH_PLATE.txt)
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

    echo "starting ipyrad plate $P ..." >> $log
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
                
    echo "finished plate ${P}" >> $log
        echo -n "    " >> $log
        date >> $log
 
 done
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








    

