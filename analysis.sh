#!/bin/bash

# FILENAME: analysis.sh

#SBATCH -A bharpur
#SBATCH --ntasks=8
#SBATCH --time=03:00:00
#SBATCH --job-name gbs_analysis
#SBATCH --output=/home/dryals/ryals/honeybee-gbs/outputs/analysis.out
#SBATCH --error=/home/dryals/ryals/honeybee-gbs/outputs/analysis.out

#Dylan Ryals 15 MAY 2023
#last edited 22 May 2023

date

#iPyradS

#testing with one lane of 96 samples

module load biocontainers vcftools bcftools plink anaconda r


# %coverage ???
#     #23-T3w04-mapped-sorted.bam as an example file
#     cd $CLUSTER_SCRATCH/gbs/ipyrad/test-gbs_refmapping
#     samtools depth -a 23-II42w06-mapped-sorted.bam > ../../analysis/II42w06.depth
#     
#     cd $CLUSTER_SCRATCH/gbs/analysis
#     
#     echo "$( grep -Pc "\t0" II42w06.depth ) / $( wc -l II42w06.depth | awk '{print $1}')" | bc -l
#     

#format and fix the input file
    #set some global variables 
    chrsLong=$( cat /home/dryals/ryals/diversity/chrfilter.txt | tr '\n' ',' )
    refs=/depot/bharpur/data/popgenomes/HarpurPNAS/output_snp.vcf.gz
    rename=/home/dryals/ryals/admixPipeline/chrsrename.txt
    chrsShort=$( awk '{print $2}' $rename | tr '\n' ' ' )
    
#     cd $CLUSTER_SCRATCH/gbs/ipyrad/bag13-final_outfiles
#     cp bag13-final.vcf ../../analysis
#     
     cd $CLUSTER_SCRATCH/gbs/analysis
   
# echo "sorting input..."
#sort the input
#     cat bag13-final.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > bag13-sorted.vcf
#     conda activate ipyrad 
#     bgzip bag13-sorted.vcf
#     conda deactivate
    
#     bcftools index -c bag13-sorted.vcf.gz 
#     bcftools annotate bag13-sorted.vcf.gz --rename-chrs $rename --threads $SLURM_NTASKS --force -Ob -o bag13.bcf.gz
#     bcftools index -c bag13.bcf.gz 
#     
#     
# echo "filtering input..."
# #filter the input
#     #remove low-qual samples and drones
#     bcftools view bag13.bcf.gz -S ^bag13-remove.txt -q 0.01:minor -e 'F_MISSING>0.10' --threads $SLURM_NTASKS -Ob -o bag13-filter.bcf.gz
#     bcftools index -c bag13-filter.bcf.gz
#     
#     #depth and coverage stats
#     bcftools query -l bag13-filter.bcf.gz > bag13-filter.names
#     bcftools query -f'%CHROM\t%POS\t%DP\n' bag13-filter.bcf.gz > bag13-filter.depth
#     #rscript ...
# 
# #grab reference file 
# #     bcftools index -c reference.bcf.gz
# 
# 
# echo "merging with reference..."
# #merge in the references
#     #save vcf sites
#     bcftools query bag13-filter.bcf.gz -f'%CHROM\t%POS\n' -o gbs.sites
#     #filter 
#     bcftools view reference.bcf.gz -T gbs.sites --threads $SLURM_NTASKS -Ob -o reference-filter.bcf.gz
#         #index
#         bcftools index -c reference-filter.bcf.gz
#     #merge
#     bcftools merge reference-filter.bcf.gz bag13-filter.bcf.gz -m snps -Ou | bcftools norm -m +snps -Ou | bcftools view -M2 -m2 --threads $SLURM_NTASKS -Ob -o admix.bcf.gz
#     
#     bcftools index -c admix.bcf.gz
# 
# echo "finding AIMs..."
# #find AIMs
#     #get popfrq
#     mkdir -p aim
#     cd aim
#     for pop in A C M O
#     do
#         bcftools view ../admix.bcf.gz -S ~/ryals/admixPipeline/references/${pop}.txt -Ou | bcftools +fill-tags | bcftools query -f'%CHROM\t%POS\t%AF\n' -o ${pop}.frq
#         
#         awk '{print $3}' ${pop}.frq > ${pop}.tmp
#     done
# 
#     paste A.frq C.tmp M.tmp O.tmp > ref.popfrq
#     rm *.tmp *.frq
#     
#     #calc AIM
#     cd ~/ryals/honeybee-gbs
#     R --vanilla --no-save --no-echo --silent < aimIa_v2.R
#     
#     cd $CLUSTER_SCRATCH/gbs/analysis/aim
#     
#     #format
#     sort -k3 -gr bag13.ia > bag13-sorted.ia
#     awk 'OFS=":" {print$1, $2}' bag13-sorted.ia | head -n 5000 > plink_aim.txt
#     
# 
# #plink
# echo "running plink..."
#     cd $CLUSTER_SCRATCH/gbs/analysis
#     #sanitize variant ids 
#     bcftools annotate admix.bcf.gz -I '.' --threads $SLURM_NTASKS -Ob -o admix2.bcf.gz
#     
#     #consider ld pruning...
#     plink --bcf admix2.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --extract aim/plink_aim.txt --threads $SLURM_NTASKS --silent --out admix
#     
# echo "estimating admixture"
# #run admixture
#     #create pop file
#     cd ~/ryals/honeybee-gbs/
#     R --vanilla --no-save --no-echo --silent < makeAdmixPop.R
#     sleep 3
#     cd -
#     
#     ADMIX=/depot/bharpur/apps/admixture/admixture
#     $ADMIX admix.bed 4 -j${SLURM_NTASKS} --cv=20 --supervised > admix.out
#     
# echo "estimating relatedness..."
# #run NGSremix
#     NGSremix=/depot/bharpur/apps/NGSremix/src/NGSremix
#     #calculate sample start and end
#     end=$( wc -l admix.pop | awk '{print $1}' )
#     samps=$( grep -c "-" admix.pop )
#     start=$( echo "$end - $samps + 1" | bc )
#     
#     $NGSremix -plink admix -f admix.4.P -q admix.4.Q -P 1 -o bag13.rel -select ${start}-${end} &
    
#GRM
    #filter a test file
    bcftools view bag13-filter.bcf.gz -r 11 --threads $SLURM_NTASKS \
    -Oz -o bag13-filter-short.vcf.gz
    
    #R script
    
    
    
    
#### additional analysis ###
    
# #LD decay on Indiana honeybees 
#     #bell
#     cd $CLUSTER_SCRATCH/gbs/analysis
#     
# 
# plink --bcf $CLUSTER_SCRATCH/pipeline/allsamp.filter.bcf.gz --make-bed --allow-extra-chr --chr-set 16 no-xy -chr $chrsShort --set-missing-var-ids @:# --keep ~/ryals/honeybee-gbs/data/plink_indy.txt --thin 0.1 --r2 gz --ld-window 100 --ld-window-kb 2000 --ld-window-r2 0 --silent --threads $SLURM_NTASKS --out indy2 &

    #Rscript...
        


####
echo "DONE"
date








    

