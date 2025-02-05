#for parallel processing
args = commandArgs(trailingOnly=TRUE)
chrno = args[1]

suppressMessages(library(tidyverse))

setwd("/scratch/negishi/dryals/gbs/23CBH/analysis")


#read in allele frequencies ...
af = read.delim("23CBH.frq", header = F) 
colnames(af) = c("chr", "pos", "f")

#read in genotypes
  vcf.raw = read.vcfR("data/23CBH-filter.vcf")
  
  #convert to dosage
  gt = extract.gt(vcf.raw)
  
  gt2 = gt
  gt2[gt == "0/0"] = "0"
  gt2[gt2 == "1/1"] = "2"
  gt2[gt2 == "1/0" | gt2 == "0/1"] = "1"
  gt2 = as.data.frame(gt2)
  
  rm(gt, vcf.raw)
  


#queen and average worker gt
#####
#create object for queen and worker group gt
qgt = matrix(nrow = nrow(af), ncol = length(unique(samples$queen_id))) %>% 
  as.data.frame
colnames(qgt) = unique(samples$queen_id)
wgt = qgt

#loop through sites
#could also be coded by creating one new obj for each step...
#takes like 20min with 7k sites
#TODO: possible to call one queen allele if not both???
for (i in (1:nrow(qgt))){
  
  p = af$f[i]
  #queen genotype likelihood
  qgl = data.frame("w2" = c(p, p/2, 0), 
                   "w1" = c(1-p, 1/2, p), 
                   "w0" = c(0, (1-p)/2, 1-p))
  rownames(qgl) = c("q2", "q1", "q0")
  
  #loop through colonies
  for (c in colnames(qgt)){
    #pull workers
    cworkers.cols = grepl(c, colnames(gt2))
    #estimate queen gt and error accounting for missing genotypes
    cworkers.gt = as.numeric(gt2[i, cworkers.cols])
    cworkers.gt = cworkers.gt[!is.na(cworkers.gt)]
    nworker = length(cworkers.gt)
    #require 5 non-missing
    if(nworker < 6){
      qgt[i, c] = NA
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
    
    qgt[i, c] = c(2, 1, 0)[which.max(c(lq2, lq1, lq0))]
    
    #TODO: calculate average worker gt
  }
}