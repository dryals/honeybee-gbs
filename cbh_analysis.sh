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

# 
# 
# #select samples to remove
# cd $CLUSTER_SCRATCH/gbs/23CBH
# #head -n 1 23CBH_1/*consens/s5* > s5summary.txt
# sed -s 1d  */*consens/s5* >> s5summary.txt
# 
# #move results into working directory
#     #cp 23CBH_1/varcall-update-split*/*.vcf analysis
# 
# #sort and filter vcfs
#     #this is probably too many steps :/
#     cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
#     echo "sorting input..."
#     #sort the input
#         for i in 1 2 
#         do
#             echo "    working on split ${i} ..."
#             cat varcall-update-split${i}.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' \
#                 > split${i}-sorted.vcf
#             conda activate ipyrad 
#             bgzip split${i}-sorted.vcf
#             conda deactivate
#             
#             bcftools index -c split${i}-sorted.vcf.gz 
#             bcftools annotate split${i}-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force \
#                 -Ob -o split${i}.bcf.gz
#             bcftools index -c split${i}.bcf.gz
#         done
#     
#     echo "merging..."
#         cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
#         
#         bcftools merge split1.bcf.gz split2.bcf.gz -0 --threads $SLURM_NTASKS -Ob -o 23CBH-updated.bcf.gz
#         bcftools index -c 23CBH-updated.bcf.gz
#         
    echo "filtering input..."
    #filter the input
        #remove missing, keep all alleles
#         bcftools view 23CBH-updated.bcf.gz -M2 -q 0.01:minor -e 'F_MISSING>0.10' \
#             -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 \
#             --threads $SLURM_NTASKS -Ob -o 23CBH-updated-filter.bcf.gz 
#             
#         bcftools index -c 23CBH-updated-filter.bcf.gz
        cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
        
        bcftools view 23CBH-updated.bcf.gz -M2 -m2 -e 'F_MISSING>0.10' \
            -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 \
            --threads $SLURM_NTASKS -Ob -o 23CBH-updated-geno.bcf.gz 
            
        bcftools index -c 23CBH-updated-geno.bcf.gz
        
        bcftools view 23CBH-updated-geno.bcf.gz -q 0.001:minor \
            --threads $SLURM_NTASKS -Oz -o 23CBH-ap.vcf.gz 
            
        bcftools index -c 23CBH-ap.vcf.gz
        
        #TODO: remove bad inds!!!!!!@!!!!!!!!!!!!
        
#         
#         #how many sites?
#         bcftools view 23CBH-updated-filter.bcf.gz | grep -v "#" | wc -l
#             #21k sites ... not so great :/
#             #try less stringent filtering (perhaps)... more testing required...
#         
# #         #maf filter
# #         bcftools view 23CBH-filter.bcf.gz -q 0.01:minor --threads $SLURM_NTASKS -Ob -o 23CBH-maf.bcf.gz
# #         bcftools index -c 23CBH-maf.bcf.gz
#         
#         
#         #pull vcf
#         bcftools view 23CBH-updated-filter.bcf.gz -Ov -o 23CBH-updated-filter.vcf
#         
#         
#         #pull sample names
#         grep "#CHROM" -m 1 23CBH-updated-filter.vcf > header.txt
#

#KING relatedness
    cd $CLUSTER_SCRATCH/gbs/23CBH/analysis/plink
    
    module load biocontainers plink
    
    plink --bcf ../23CBH-updated-geno.bcf.gz --make-bed \
        --geno 0.1 --mind 0.4 --maf 0.005 \
        --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
        --set-missing-var-ids @:# \
        --out 23CBH
        
    module purge
    module load biocontainers plink2

    plink2 --bfile 23CBH --make-king square --out 23CBH

    #module purge
    
    awk -v OFS='_' '{print $1, $2}' 23CBH.fam > keep.samps
    
    cd ..
    
    bcftools view 23CBH-updated-geno.bcf.gz -q 0.01:minor \
            -S ap/goodworkers.txt \
            --threads $SLURM_NTASKS -Ov -o 23CBH-ap.vcf 



#     echo "calculataing allele freqs..."
#     # see R script ... 
#     bcftools view 23CBH-updated-filter.bcf.gz -S ~/ryals/honeybee-gbs/data/balanceSet.txt -Ou | \
#             bcftools +fill-tags -Ob -o 23CBH-balanced.bcf.gz
#             
#     bcftools index -c 23CBH-balanced.bcf.gz
#             
#     bcftools view 23CBH-balanced.bcf.gz -Ou | bcftools +fill-tags -Ou | \
#         bcftools query -f'%CHROM\t%POS\t%AF\n' -o 23CBH-updated.frq


#...using alpha peel!
    #conda environment
#         module purge
#         module load anaconda
#         
#         conda create -n AlphaPeel
#         conda activate AlphaPeel
#         conda install pip pip setuptools wheel
#         pip install AlphaPeel
#         pip install cgi
#         pip install Cython=0.29.1
#         pip install AlphaPlinkPython
        
        module load anaconda
        cd $CLUSTER_SCRATCH/gbs/23CBH/analysis/ap
        conda activate AlphaPeel
        
        AlphaPeel -genotypes APgeno.txt -pedigree APped.txt \
            -out ap2
            
        


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
    

# #load into plink
#     #create worker
#     plink --vcf 23CBH-filter.vcf --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
#             --set-missing-var-ids @:# --maf 0.01 --mind 0.5 --geno 0.2 \
#             --threads $SLURM_NTASKS --out plink/workers
#             
#     cd plink
#     plink --bfile workers --read-freq balworkers.frq --make-rel square --out workers
#     
#     #queens
#     plink --file qgt --make-bed --out plink/qraw
#     cd plink
#     plink --bfile qraw --make-bed --geno 0.2 --mind 0.5 --maf 0.05  \
#         --make-rel square --out qfilter
#     
#         #optional: read balanced worker set to calc af
#     plink --bcf 23CBH-balanced.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort \
#             --set-missing-var-ids @:# \
#             --threads $SLURM_NTASKS --out plink/balworkers
#     cd plink
#     plink --bfile balworkers --freq --out balworkers
#     plink --bfile qraw --make-bed --read-freq balworkers.frq --geno 0.2 --mind 0.5 --maf 0.05 \
#         --make-rel square --out qfilter
#     
#     
#     plink --bfile qfilter --threads $SLURM_NTASKS --maf 0.05 --pca 100 --out pca
#     
#     #queens and workers together
#     cd $CLUSTER_SCRATCH/gbs/23CBH/analysis
#     cat qgt.ped wgt.ped > qw.ped
#     cp qgt.map qw.map
#     
#      plink --file qw --make-bed --out plink/qwraw
#      
#      cd plink
#      plink --bfile qwraw --make-bed --read-freq balworkers.frq --geno 0.2 --mind 0.5 --maf 0.05 \
#         --make-rel square --out qwfilter
#     
#     
#     
#     #output for blup
#         #there has GOT to be a better way ...
#     plink --bfile qfilter --recode A --out test
#         awk '{$1=$3=$4=$5=$6=""; print $0}' test.raw | tail -n +2 |\
#             sed -r 's/[NA]+/5/g' | sed 's/ \{2,\}/!/g' | tr -d ' ' | tr '!' ' ' > blup.raw
#             
#     

#     q
#         #depth and coverage stats
#         bcftools query -l bag13-filter.bcf.gz > bag13-filter.names
#         bcftools query -f'%CHROM\t%POS\t%DP\n' bag13-filter.bcf.gz > bag13-filter.depth
#         #rscript ...

    

####
echo "DONE"
date








    

