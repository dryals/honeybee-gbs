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



# #select samples to remove
# cd $CLUSTER_SCRATCH/gbs/23CBH
# 
# 
# head -n 1 23CBH_1/*consens/s5* > s5summary.txt
# sed -s 1d  */*consens/s5* >> s5summary.txt

#move results into working directory
    #...

#sort and filter vcf
    #this is probably too many steps :/
    cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
#     echo "sorting input..."
#     #sort the input
#         cat 23CBH.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > 23CBH-sorted.vcf
#         conda activate ipyrad 
#         bgzip 23CBH-sorted.vcf
#         conda deactivate
#         
#         bcftools index -c 23CBH-sorted.vcf.gz 
#         bcftools annotate 23CBH-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force \
#             -Ob -o 23CBH.bcf.gz
#         bcftools index -c 23CBH.bcf.gz
#         
#         
    echo "filtering input..."
    #filter the input
        #remove missing, keep all alleles (no MAF filter)
        bcftools view 23CBH.bcf.gz -M2 -q 0.001:minor -e 'F_MISSING>0.10' \
            -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 \
            --threads $SLURM_NTASKS -Ob -o 23CBH-filter.bcf.gz
            
        bcftools index -c 23CBH-filter.bcf.gz
        
        #maf filter
        bcftools view 23CBH-filter.bcf.gz -q 0.01:minor --threads $SLURM_NTASKS -Ob -o 23CBH-maf.bcf.gz
        bcftools index -c 23CBH-maf.bcf.gz
        
        
        #pull vcf
        #bcftools view 23CBH-filter.bcf.gz -Ov -o 23CBH-filter.vcf
        bcftools view 23CBH-maf.bcf.gz -Ov -o 23CBH-maf.vcf
        
        #pull sample names
        #grep "#CHROM" -m 1 23CBH.vcf > header.txt
        
        
    echo "calculataing allele freqs..."
    bcftools view 23CBH-maf.bcf.gz -S ~/ryals/honeybee-gbs/data/balanceSet.txt -Ou | \
            bcftools +fill-tags -Ob -o 23CBH-balanced.bcf.gz
            
    bcftools index -c 23CBH-balanced.bcf.gz
            
    bcftools view 23CBH-balanced.bcf.gz -Ou | bcftools +fill-tags -Ou | \
        bcftools query -f'%CHROM\t%POS\t%AF\n' -o 23CBH.frq
   
# #predict queen and average worker genotypes
    #see R script...
    
    #TODO: this isn't getting great results, double-check and debug script. 
    Rscript --vanilla --silent /home/dryals/ryals/honeybee-gbs/queencaller.R
    
    
    
#     #split chrs
#     mkdir -p chrs
#     cd chrs
#     for i in {1..16}
#     do
#         mkdir -p chr${i}
#         bcftools view ../23CBH-filter.bcf.gz -r $i -Ov -o chr${i}/23CBH-filter-${i}.vcf
#     
#     done
    

#load into plink
    #create worker
    plink --vcf 23CBH-filter.vcf --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
            --set-missing-var-ids @:# --maf 0.01 --mind 0.5 --geno 0.2 \
            --threads $SLURM_NTASKS --out plink/workers
            
    cd plink
    plink --bfile workers --read-freq balworkers.frq --make-rel square --out workers
    
    #queens
    plink --file qgt --make-bed --out plink/qraw
    cd plink
    plink --bfile qraw --make-bed --geno 0.2 --mind 0.5 --maf 0.05  \
        --make-rel square --out qfilter
    
        #optional: read balanced worker set to calc af
    plink --bcf 23CBH-balanced.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
            --set-missing-var-ids @:# \
            --threads $SLURM_NTASKS --out plink/balworkers
    cd plink
    plink --bfile balworkers --freq --out balworkers
    plink --bfile qraw --make-bed --read-freq balworkers.frq --geno 0.2 --mind 0.5 --maf 0.05 \
        --make-rel square --out qfilter
    
    
    
    plink --bfile qfilter --threads $SLURM_NTASKS --maf 0.05 --pca 100 --out pca
    
    #queens and workers together
#     cat qraw.fam workers.fam > testqw.fam
#     cp workers.fam > 
#     plink testqw
    
    
    #output for blup
        #there has GOT to be a better way ...
    plink --bfile qfilter --recode A --out test
        awk '{$1=$3=$4=$5=$6=""; print $0}' test.raw | tail -n +2 |\
            sed -r 's/[NA]+/5/g' | sed 's/ \{2,\}/!/g' | tr -d ' ' | tr '!' ' ' > blup.raw
            
    

#     q
#         #depth and coverage stats
#         bcftools query -l bag13-filter.bcf.gz > bag13-filter.names
#         bcftools query -f'%CHROM\t%POS\t%DP\n' bag13-filter.bcf.gz > bag13-filter.depth
#         #rscript ...

    

####
echo "DONE"
date








    

