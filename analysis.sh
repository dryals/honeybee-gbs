#!/bin/bash

# FILENAME: analysis.sh

#SBATCH -A bharpur
#SBATCH --ntasks=4
#SBATCH --time=1-00:00:00
#SBATCH --job-name gbs_analysis
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/analysis.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/analysis.out

#Dylan Ryals 15 MAY 2023
#last edited

date

#iPyradS

#testing with one lane of 96 samples

module load biocontainers vcftools bcftools plink anaconda



#format and fix the input file
    #set some global variables 
    chrsLong=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
    refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz
    rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )
        
    bcftools sort test-branch.vcf.gz -Ob -o test.bcf.gz
    
#sort the input
    cat test-branch.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > test-sorted.vcf
    conda activate ipyrad 
    bgzip test-sorted.vcf
    
    bcftools index -c test-sorted.vcf.gz &
    bcftools annotate test-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force -Ob -o test.bcf.gz
    bcftools index -c test.bcf.gz &
    
#filter the input
    bcftools view test.bcf.gz -q 0.05:minor -e 'F_MISSING>0.1' --threads $SLURM_NTASKS -Ob -o test-filter.bcf.gz
    bcftools index -c test-filter.bcf.gz

#grab reference file 
    bcftools index -c reference.bcf.gz


#merge in the references
    #save vcf sites
    bcftools query test-filter.bcf.gz -f'%CHROM\t%POS\n' -o gbs.sites
    #filter 
    bcftools view reference.bcf.gz -T gbs.sites --threads $SLURM_NTASKS -Ob -o reference-filter.bcf.gz
        #index
        bcftools index -c reference-filter.bcf.gz
    #merge
    bcftools merge reference-filter.bcf.gz test-filter.bcf.gz -m snps -Ou | bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o admix.bcf.gz
    
    bcftools index -c admix.bcf.gz

#find AIMs??? or nah???
    mkdir -p aim
    cd aim
    for pop in A C M O
    do
        bcftools view ../reference-filter.bcf.gz -S ~/ryals/admixPipeline/references/${pop}.txt -Ou | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%AF\n' -o ${pop}.frq
        
        awk '{print $3}' ${pop}.frq > ${pop}.tmp
    done

    paste A.frq C.tmp M.tmp O.tmp > ref.popfrq
    rm *.tmp *.frq



    
#plink
    plink --bcf admix.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --threads $SLURM_NTASKS --silent --out admix
    
        cd /home/dryals/ryals/honeybee-gbs
        #create pop file
        R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
        
#run admixture
    cd $CLUSTER_SCRATCH/gbs/analysis
    ADMIX=/depot/bharpur/apps/admixture/admixture
    #remember to change this for the correct k value!!
    $ADMIX admix.bed 4 -j4 --cv=20 --supervised > admix.out
        
#run NGSremix
    #this takes a hot minute
    NGSremix=/depot/bharpur/apps/NGSremix/src/NGSremix
    #make sure to edit 'select' to remove the references... this can be done programmatically
    $NGSremix -plink admix -f admix.4.P -q admix.4.Q -P 1 -o test.rel -select 75-169
    

    
    



####
echo "DONE"
date








    

