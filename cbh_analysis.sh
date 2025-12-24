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
        
        #TODO: remove bad inds!!!!!!!!!!!!!!!!!!
        
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
            
    #output site list for ap
    grep -v "#" 23CBH-ap.vcf | awk '{print $1, $2}' > 23CBH-ap.sites
    
    cat 23CBH-ap.sites | tr " " "\t" > 23CBH-ap-tr.sites
    
    #output full sample list for queencaller
    bcftools view 23CBH-updated-geno.bcf.gz \
        -T 23CBH-ap-tr.sites \
        --threads $SLURM_NTASKS -Ov -o 23CBH-queencaller.vcf
    
    
    #find total map lengths
    cp /depot/bharpur/data/ref_genomes/AMEL/resources/recombination/AMELMarley.renamed.fixed.txt \
        ./marey.txt
        
        cat marey.txt | tr "\"" " " > marey2.txt
        
#         #done in R bacause awk is a pain
#         R
#         #print length of ea chrs
#         marey = read.delim("marey2.txt", header = T, sep = "\t")
#         library(dplyr)
#         chr.big = marey %>% group_by(map) %>% arrange(desc(gen)) %>% slice(1) %>% select(map, gen) %>% 
#             mutate(gen = gen / 100)
#         write.table(file="chrlength.txt", chr.big, quote = F, col.names = F, row.names = F)
#         #print start and stop indicies
#         sites = read.delim("23CBH-ap.sites", header = F, sep = "")
#         
#         chridx = data.frame(chr = unique(sites$V1),
#                             start = NA,
#                             stop = NA
#                             )
#         chridx$start[1] = 1
#         chridx$stop[1] = which.max( sites$V2[sites$V1 == 1] )
#         for( i in 2:16){
#         
#         chridx$start[i] = chridx$stop[(i-1)] + 1
#         chridx$stop[i] = chridx$stop[(i-1)] + which.max( sites$V2[sites$V1 == i] )
#         
#         }
#         write.table(chridx, file = "chridx.txt", col.names = F, row.names = F, quote = F)
#         #write out map files
#         for (i in 1:16){
#             sites.tmp = sites %>% filter(V1 == i)
#                 colnames(sites.tmp) = c("chr","bp")
#                 sites.tmp$name = paste0("snp", 1:nrow(sites.tmp))
#                 
#             write.table(sites.tmp %>% select(chr, name, bp), 
#                 file = paste0("ap/chr", i, "/map.txt"),
#                 row.names = F, col.names = F, quote = F
#             )
#         }
#         
#         quit(save = "no")
        

        




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

#download devel branch 
    cd ~/ryals
    git clone --recurse-submodules -b devel https://github.com/AlphaGenes/AlphaPeel.git
    cd AlphaPeel
    conda create -n alphaPeel3
    conda activate alphaPeel3
    conda install pip
    python3 -m pip install --upgrade build
    python3 -m build
    python3 -m pip install dist/alphapeel*.whl
    
    
    
        
        module load anaconda
        cd $CLUSTER_SCRATCH/gbs/23CBH/analysis/ap
        conda activate alphaPeel3
        
        #TODO
            #reach out to Jana for help!
            #use Purdue metafounder
            #include suspect colonies as unplaced indiv
            #does AP generate good parental genotypes? use instead of queencaller
            #hybrid peeling? phasing?
                #input a map file?
                #use phasing information for queen caller????
            #use imputation? run from raw seq data?
        

         #loop through chrs ... would be better in an array job

         for chr in {1..16}
         do
            #take chr parameters
            start=$( sed "${chr}q;d" ../chridx.txt | awk '{print $2}')
            stop=$( sed "${chr}q;d" ../chridx.txt | awk '{print $3}' )
            len=$( sed "${chr}q;d" ../chrlength.txt | awk '{print $2}' )
            
            echo "$chr : $start , $stop , $len"
            
            mkdir -p chr${chr}
            
            #run ap for chr
            AlphaPeel -geno_file APgeno-reduced.txt \
                -method multi \
                -ped_file APped-reduced.txt \
                -out_file chr${chr}/ap \
                -rec_length $len -start_snp $start -stop_snp $stop \
                -n_cycle 8 \
                -n_thread $SLURM_NTASKS \
                -n_io_thread 8 \
                -no_dosage \
                -est_alt_allele_prob \
                -update_alt_allele_prob \
                -alt_allele_prob \
                -geno_threshold 0.34 -geno \
                -map_file chr${chr}/map.txt

         done
         
    #combine AF and GT into dingle files
        #AF
        tail -n +2 chr1/ap.alt_allele_prob.txt > all_founderAF2.txt
        #GT
        cat chr1/ap.geno_0.34.txt > all_GT.txt
        for chr in {2..16}
        do
            #AF
            tail -n +2 chr${chr}/ap.alt_allele_prob.txt >> all_founderAF2.txt
            #GT
            paste all_GT.txt \
            <( awk '{OFS=" "; $1=""; gsub(/[[:space:]]+/, " "); print $0}' chr${chr}/ap.geno_0.34.txt ) > tmp
            cat tmp > all_GT.txt
        done
        

        
        
        #compare against naive afs
        cd ..
    bcftools view 23CBH-ap.vcf | bcftools +fill-tags | \
        bcftools query -f'%CHROM\t%POS\t%AF\n' -o 23CBH-ap.frq
        
        
        


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








    

