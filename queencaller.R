#for parallel processing
# args = commandArgs(trailingOnly=TRUE)
# chrno = args[1]

suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))

select = dplyr::select

setwd("/scratch/negishi/dryals/gbs/23CBH/analysis")

#read in sample info
samples = read.delim("header.txt", header = F) %>% t() %>% 
  as.data.frame() %>% 
  rename(sample_id = 1) %>% 
  filter(grepl("23CBH", sample_id)) %>% 
  mutate(queen_id = gsub("_[0-9]*", "", sample_id),
         colony_id = gsub("23CBH", "", queen_id)) 

#read in allele frequencies ...
af = read.delim("23CBH.frq", header = F) 
  colnames(af) = c("chr", "pos", "f")
  #af = af %>% filter(chr == chrno)

#read in genotypes
  #vcf.raw = read.vcfR(paste0("chrs/chr",chrno,"/23CBH-filter-", chrno, ".vcf"))
  vcf.raw = read.vcfR("23CBH-filter.vcf")
  
  
  #convert to dosage
  gt = extract.gt(vcf.raw)
  
  gt2 = gt
  gt2[gt == "0/0"] = "0"
  gt2[gt2 == "1/1"] = "2"
  gt2[gt2 == "1/0" | gt2 == "0/1"] = "1"
  gt2 = as.data.frame(gt2)
  
  rm(gt, vcf.raw)
  
  
  
#randomly sample for testing
  # set.seed(123)
  # fewsites = sample(1:nrow(gt2), 4000) %>% sort()
  # 
  # gt2 = gt2[fewsites, ]
  # af = af[fewsites, ]

  
callqueen = function(CHR){
  #filter 
  chrrows = af$chr == CHR
  af.chr = af[chrrows,]
  gt2.chr = gt2[chrrows,]
  
  #create object for queen and worker group gt
  
  qgt = matrix(nrow = nrow(af.chr), 
               ncol = length(unique(samples$queen_id))) %>% 
    as.data.frame
  colnames(qgt) = unique(samples$queen_id)
  
  
  
  #loop through sites
  for (i in (1:nrow(qgt))){
    p = af.chr$f[i]
    #queen genotype likelihood
    qgl = data.frame("w2" = c(p, p/2, 0), 
                     "w1" = c(1-p, 1/2, p), 
                     "w0" = c(0, (1-p)/2, 1-p))
    rownames(qgl) = c("q2", "q1", "q0")
    
    #loop through colonies
    for (col in colnames(qgt)){
      #pull workers
      cworkers.cols = grepl(col, colnames(gt2.chr))
      #estimate queen gt and error accounting for missing genotypes
      cworkers.gt = as.numeric(gt2.chr[i, cworkers.cols])
      cworkers.gt = cworkers.gt[!is.na(cworkers.gt)]
      nworker = length(cworkers.gt)
      #require 5 non-missing
      if(nworker < 6){
        qgt[i, col] = NA
        next
      }
      cworkers.counts = c(sum(cworkers.gt == 2), 
                          sum(cworkers.gt == 1), 
                          sum(cworkers.gt == 0))
      #multinomial
      #lik of q2
      lq2 = dmultinom(cworkers.counts, nworker, as.numeric(qgl[1,]))
      #likelihood of q1
      lq1 = dmultinom(cworkers.counts, nworker, as.numeric(qgl[2,]))
      #lik of q0
      lq0 = dmultinom(cworkers.counts, nworker, as.numeric(qgl[3,]))
      
      qgt[i, col] = c(2, 1, 0)[which.max(c(lq2, lq1, lq0))]
      
    }
  }
  return(qgt)
}  

#run in parallel
  registerDoParallel(cores = 16)
  
  #starttime = Sys.time()
  
    finalqgt = foreach(i=1:16, .combine=rbind) %dopar%
      callqueen(i)
  
  #Sys.time() - starttime


#write object 
  #write queen for plink .ped
  pt1 = data.frame(IID = names(finalqgt)) %>% 
    left_join(samples %>% select(IID = queen_id, FID = queen_id), multiple = 'first') %>% 
    mutate(FID = gsub("-", "", FID)) %>% 
    select(FID, IID) %>% 
    mutate(father = 0, mother = 0, sex = 2, pheno = 0)
  
  pt2 = t(finalqgt)
  pt2[pt2 == 0] = "TT"
  pt2[pt2 == 1] = "AT"
  pt2[pt2 == 2] = "AA"
  pt2[is.na(pt2)] = "00"
  
  qgt.out = cbind(pt1, pt2)
  
  write.table(qgt.out, "qgt.ped", 
              row.names = F, col.names = F, quote = F, sep = " ")
  #write .map file
  map = af %>% 
    mutate(variant = paste0(chr, ":", pos),
           cm = 0) %>% 
    select(chr, variant, cm, pos)
  
  write.table(map, "qgt.map", 
              row.names = F, col.names = F, quote = F, sep = " ")



#mean(af$f[af$chr == i])


  
  
# #write out in brief form
# write.table(af %>% select(chr, pos), 
#             file = paste0("chrs/chr",chrno,"/", chrno, ".sites"),
#             col.names = F, quote = F, row.names = F)
# 
# write.table(qgt, 
#             file = paste0("chrs/chr",chrno,"/", chrno, ".geno"),
#             col.names = F, quote = F, row.names = F, sep = " ")
# 
# write.table(colnames(qgt), 
#             file = paste0("chrs/chr",chrno,"/", chrno, ".names"),
#             col.names = F, quote = F, row.names = F, sep = " ")
# 
# 
#   
