#!/bin/bash

# FILENAME: 24CBH_poolCall.sh

#SBATCH -A bharpur
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1-00:00:00
#SBATCH --partition cpu
#SBATCH --job-name poolCall
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/pool.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/pool.out

#Dylan Ryals 25 AUG 2025
#last edited 27 AUG 2025

#determining sequence depth and ref allele count for pooled samples

#testing on 2023 data

date
#####

module load biocontainers samtools bcftools


# #grab some sites
#     #a quick test...
# 
#     cd /scratch/negishi/dryals/gbs/24CBH/analysis
# 
#     samtools mpileup -b test2.samps \
#         -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
#         -C 50 -q 20 -Q 20 -l test.sites \
#         -a -o pileuptest2.out
#         
#     #use R script to translate into beethoven format
#     
#     
# #testing on a single plate with all sites identified in 2023
#     awk '{print $1, $2}' /scratch/negishi/dryals/gbs/23CBH/analysis/23CBH-updated.frq > 23CBH.sites
#     
#     #translate numeric chr to codes (annoying)
#         awk '{print $2, $1}' ~/ryals/admixPipeline/chrsrename.txt  > chrsrename.txt
#         awk '{print $1}' 23CBH.sites > tmp.1
#         
#         awk '
#             NR==FNR { map[$1] = $2; next }
#             $0 in map { $0 = map[$0] }
#             { print }
#             ' chrsrename.txt tmp.1 > tmp.2
#         
#         awk '{print $2}' 23CBH.sites > tmp.3
#         paste tmp.2 tmp.3 > 23CBH-t.sites
#         rm tmp.*
#     
#     #try on 4 samples and time
#     date
#     
#         samtools mpileup -b test2.samps \
#         -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
#         -C 50 -q 20 -Q 20 -l 23CBH-t.sites \
#         -a -o pileuptest3.out
#     
#     date
#      #12 sec for 4 samps; 3sec/samp ; 96 samps in 5min ; 400 samps in 20min
#         #96 samples takes longer ... 28min
#     
# #whole plate
#     ls $CLUSTER_SCRATCH/gbs/24CBH/24CBH_3/24CBH_3_refmapping/*.bam > test3.samps
# 
#     date
# 
#     samtools mpileup -b test3.samps \
#     -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
#     -C 50 -q 20 -Q 20 -l 23CBH-t.sites \
#     -o pileuptest4.out
# 
#     date
#     
#parallel processig by sample probably, pasting all togehter at end


#all 2024 data together, will likely take ages to run
    echo "mpileup..."
    cd /scratch/negishi/dryals/gbs/24CBH/analysis
    
    #list all samples
    ls $CLUSTER_SCRATCH/gbs/24CBH/24CBH_*/*_refmapping/*.bam > 24CBH.bamlist
    
    #remove singlebee samples, these should be added to previous year run
    grep -vE "*23CBH[0-9]{3}_[0-9]-mapped*" 24CBH.bamlist > 24CBHpool.bamlist
    
    samtools mpileup -b 24CBHpool.bamlist \
    -f /depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna \
    -l 23CBH-t.sites \
    -C 50 -q 20 -Q 20 -d 200 \
    -a -o 24CBH.mpileup


#ON LOCAL MACHINE

#     PPP=/home/dylan/Documents/bees/harpurlab/project/gensel/poPoolation2/popoolation2_1201
# 
#     java -ea -Xmx7g -jar $PPP/mpileup2sync.jar \
#     --input pileuptest4.out --output test4.sync \
#     --fastq-type sanger --min-qual 20 --threads 8

    

    
######
echo "DONE"
date








    

