#!/bin/bash

# FILENAME: cbh_analysis.sh

#SBATCH -A bharpur
#SBATCH --ntasks=16
#SBATCH --time=05:00:00
#SBATCH --job-name gbs_analysis
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/analysis_cbh.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/analysis_cbh.out

#Dylan Ryals 17 JAN 2025
#last edited 

date


module load biocontainers vcftools bcftools plink anaconda r

#set some global variables 
chrsLong=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz
rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )


#pull a small file for testing
    cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
    echo "sorting input..."
    #sort the input
        cat 23CBH.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > 23CBH-sorted.vcf
        conda activate ipyrad 
        bgzip 23CBH-sorted.vcf
        conda deactivate
        
        bcftools index -c 23CBH-sorted.vcf.gz 
        bcftools annotate 23CBH-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force \
            -Ob -o 23CBH.bcf.gz
        bcftools index -c 23CBH.bcf.gz
        
        
    echo "filtering input..."
    #filter the input
        #remove low-qual samples and drones
        bcftools view 23CBH.bcf.gz -S -q 0.01:minor -e 'F_MISSING>0.10' --threads $SLURM_NTASKS \
            -Ob -o 23CBH-filter.bcf.gz
        bcftools index -c 23CBH-filter.bcf.gz
        
        #pull one chr for testing
        bcftools view 23CBH-filter.bcf.gz -r 11 -Ov -o 23CBH-small.vcf
        
#         #depth and coverage stats
#         bcftools query -l bag13-filter.bcf.gz > bag13-filter.names
#         bcftools query -f'%CHROM\t%POS\t%DP\n' bag13-filter.bcf.gz > bag13-filter.depth
#         #rscript ...

    

####
echo "DONE"
date








    

