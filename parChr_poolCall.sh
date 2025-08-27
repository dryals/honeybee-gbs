#!/bin/bash

# FILENAME: parChr_poolCall.sh

#SBATCH -A bharpur
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition cpu
#SBATCH --job-name par_poolCall
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pool_dump.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pool_dump.out

#Dylan Ryals 27 AUG 2025
#last edited 27 AUG 2025

#parallelized by chromosome
#determining sequence depth and ref allele count for pooled samples


log=/home/dryals/ryals/honeybee-gbs/outputs/par_pool.out

#####

module load biocontainers samtools

    cd /scratch/negishi/dryals/gbs/24CBH/analysis

    chrShort=$( echo $SLURM_ARRAY_TASK_ID )
    chrLong=$( sed "${chrShort}q;d" chrsrename.txt | awk '{print $2}' )
    
    mkdir -p chrPar/chr${chrShort}
    cd chrPar/chr${chrShort}
    
    grep $chrLong ../../23CBH-t.sites > 23CBH_chr${chrShort}.sites
    
    echo "starting chr $chrShort" >> $log
    
    samtools mpileup -b ../../24CBHpool.bamlist \
    -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
    -l 23CBH_chr${chrShort}.sites \
    -C 50 -q 20 -Q 20 -d 200 \
    -a -o 24CBHchr${chrShort}.mpileup
    
    echo "FINISHED chr $chrShort" >> $log


#ON LOCAL MACHINE

#     PPP=/home/dylan/Documents/bees/harpurlab/project/gensel/poPoolation2/popoolation2_1201
# 
#     java -ea -Xmx7g -jar $PPP/mpileup2sync.jar \
#     --input pileuptest4.out --output test4.sync \
#     --fastq-type sanger --min-qual 20 --threads 8

    
######
echo "DONE"
date








    

