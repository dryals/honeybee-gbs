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


# %coverage ???
    #23-T3w04-mapped-sorted.bam as an example file
    cd $CLUSTER_SCRATCH/gbs/ipyrad/test-gbs_refmapping
    samtools depth -a 23-II42w06-mapped-sorted.bam > ../../analysis/II42w06.depth
    
    cd $CLUSTER_SCRATCH/gbs/analysis
    
    echo "$( grep -Pc "\t0" II42w06.depth ) / $( wc -l II42w06.depth | awk '{print $1}')" | bc -l
    

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
    
    bcftools index -c test-sorted.vcf.gz 
    bcftools annotate test-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force -Ob -o test.bcf.gz
    bcftools index -c test.bcf.gz 
    
    
#filter the input
    bcftools view test.bcf.gz -q 0.05:minor -e 'F_MISSING>0.1' --threads $SLURM_NTASKS -Ob -o test-filter.bcf.gz
    bcftools index -c test-filter.bcf.gz
    
    #depth and coverage stats
    bcftools query -l test-filter.bcf.gz > test-filter.names
    bcftools query -f'%CHROM\t%POS\t%DP\n' test-filter.bcf.gz > test-filter.depth
    #rscript ...

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

#find AIMs
    #get popfrq
    mkdir -p aim
    cd aim
    for pop in A C M O
    do
        bcftools view ../reference-filter.bcf.gz -S ~/ryals/admixPipeline/references/${pop}.txt -Ou | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%AF\n' -o ${pop}.frq
        
        awk '{print $3}' ${pop}.frq > ${pop}.tmp
    done

    paste A.frq C.tmp M.tmp O.tmp > ref.popfrq
    rm *.tmp *.frq
    
    #calc AIM
    R --vanilla --no-save --no-echo --silent < ~/ryals/honeybee-gbs/aimIa_v2.R
    
    #format
    sort -k3 -gr test.ia > test-sorted.ia
    awk 'OFS=":" {print$1, $2}' test-sorted.ia | head -n 5000 > plink_aim.txt
    

#plink
    cd ..
    
    #sanitize variant ids 
    bcftools annotate admix.bcf.gz -I '.' --threads $SLURM_NTASKS -Ob -o admix2.bcf.gz
    
    plink --bcf admix2.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --extract aim/plink_aim.txt --threads $SLURM_NTASKS --silent --out admix
    
#run admixture
    #create pop file
    cd ~/ryals/honeybee-gbs/
    R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
    sleep 3
    cd -
    
    ADMIX=/depot/bharpur/apps/admixture/admixture
    $ADMIX admix.bed 4 -j${SLURM_NTASKS} --cv=20 --supervised > admix.out
        
#run NGSremix
    #this takes a hot minute
    NGSremix=/depot/bharpur/apps/NGSremix/src/NGSremix
    #make sure to edit 'select' to remove the references... this can be done programmatically
    $NGSremix -plink admix -f admix.4.P -q admix.4.Q -P 1 -o test.rel -select 75-169
    

    
#### additional analysis ###
    
#LD decay on Indiana honeybees 
    #bell
    cd $CLUSTER_SCRATCH/gbs/analysis
    

plink --bcf $CLUSTER_SCRATCH/pipeline/allsamp.filter.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --keep ~/ryals/honeybee-gbs/data/plink_indy.txt --thin 0.1 --r2 gz --ld-window 100 --ld-window-kb 2000 --ld-window-r2 0 --silent --threads $SLURM_NTASKS --out indy2 &

    #Rscript...
        


####
echo "DONE"
date








    

