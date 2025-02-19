#for parallel processing
# args = commandArgs(trailingOnly=TRUE)
# chrno = args[1]

suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(matrixcalc))

select = dplyr::select

setwd("/scratch/negishi/dryals/gbs/23CBH/analysis")

#read in sample info
samples = read.delim("header.txt", header = F) %>% t() %>% 
  as.data.frame() %>% 
  rename(sample_id = 1) %>% 
  filter(grepl("23CBH", sample_id)) %>% 
  mutate(queen_id = gsub("_[0-9]*", "", sample_id),
         colony_id = gsub("23CBH", "", queen_id)) 

  # samples %>% group_by(colony_id) %>% summarise(n = n()) %>% 
  # filter(n < 7) %>%  ungroup %>% arrange(n) 

#read in allele frequencies ...
af = read.delim("23CBH.frq", header = F) 
  colnames(af) = c("chr", "pos", "f")
  #af = af %>% filter(chr == chrno)

#read in genotypes
  #vcf.raw = read.vcfR(paste0("chrs/chr",chrno,"/23CBH-filter-", chrno, ".vcf"))
  vcf.raw = read.vcfR("23CBH-maf.vcf")
  
  
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
  # fewsites = sample(1:nrow(gt2), 5000) %>% sort()
  # 
  # oldaf = af
  # oldgt = gt2
  # 
  # gt2 = gt2[fewsites, ]
  # af = af[fewsites, ]

  
callqueen = function(CHR){
  #filter to selected chr
  chrrows = af$chr == CHR
  af.chr = af[chrrows,]
  gt2.chr = gt2[chrrows,]
  
  #create object for queen gt
  qgt = matrix(nrow = nrow(af.chr), 
               ncol = length(unique(samples$queen_id))) %>% 
    as.data.frame
  colnames(qgt) = unique(samples$queen_id)

  #loop through sites
  for (i in (1:nrow(qgt))){
    p = af.chr$f[i]
    
    #contingent worker genotypes probability for each queen genotype
      #worker probabilities in order: 2,1,0
    cwgp = rbind(c(p, 1-p, 0),
                c(p/2, 1/2, (1-p)/2),
                c(0, p, 1-p))
    #independent queen genotype probabilities (HW)
      #2,1,0
    iqgp = c(p**2, 2*p*(1-p), (1-p)**2)

    #loop through colonies
    for (col in colnames(qgt)){
      #pull genotypes
      cworkers.gt = as.numeric(gt2.chr[i, grepl(col, colnames(gt2.chr))])
      #remove NAs and count
      cworkers.gt = cworkers.gt[!is.na(cworkers.gt)]
      nworker = length(cworkers.gt)
      #require 5 non-missing
      if(nworker < 6){
        qgt[i, col] = NA
        next
      }
      
      #simulation for testing
      # rqg = c(0,1) #'real' queen genotype
      # cworkers.gt = sample(rqg, nworker, replace = T) + 
      #   sample(c(1,0), nworker, prob = c(p, 1-p), replace = T)
      # 
      # 
      #sum genotypes
      cworkers.counts = c(sum(cworkers.gt == 2), 
                          sum(cworkers.gt == 1), 
                          sum(cworkers.gt == 0))
      #cworkers.counts
      
      #bayes for each queen genotype
        #independent probability of observing worker count
          iwgp = dmultinom(cworkers.counts, nworker, iqgp)
      
        #probability of q2
        pq2 = (dmultinom(cworkers.counts, nworker, cwgp[1,]) * iqgp[1]) / 
                iwgp
        #prob of q1
        pq1 = (dmultinom(cworkers.counts, nworker, cwgp[2,]) * iqgp[2]) /
                iwgp
        #prob of q0
        pq0 = (dmultinom(cworkers.counts, nworker, cwgp[3,]) * iqgp[3]) / 
                iwgp
        
        #TODO: probs do not add to 1 ... this is bad ... 
        #sum(pq2, pq0, pq1)
        
      #write max prob
      qgt[i, col] = c(2, 1, 0)[which.max(c(pq2, pq1, pq0))]
    }
  }
  return(qgt)
}  


#callqueen(16)
#run in parallel
  registerDoParallel(cores = 16)
  
  #starttime = Sys.time()
  
    finalqgt = foreach(chr=1:16, .combine=rbind) %dopar%
      callqueen(chr)
  
    
    # sum(finalqgt1 == finalqgt2 |
    #       is.na(finalqgt1 & is.na(finalqgt2)))
    # 
    # prod(dim(finalqgt1))
    
  #Sys.time() - starttime


#write queen object 
  #write queen for plink .ped
  pt1 = data.frame(FID = names(finalqgt), IID = names(finalqgt)) %>% 
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
  
#write raw worker object
  #write worker for plink .ped
  pt1 = data.frame(FID = colnames(gt2), IID = colnames(gt2)) %>% 
    mutate(father = 0, mother = 0, sex = 2, pheno = 0)
  
  pt2 = t(gt2)
  pt2[pt2 == 0] = "TT"
  pt2[pt2 == 1] = "AT"
  pt2[pt2 == 2] = "AA"
  pt2[is.na(pt2)] = "00"
  
  wgt.out = cbind(pt1, pt2)
  
  write.table(wgt.out, "wgt.ped", 
              row.names = F, col.names = F, quote = F, sep = " ")
  #write .map file
  map = af %>% 
    mutate(variant = paste0(chr, ":", pos),
           cm = 0) %>% 
    select(chr, variant, cm, pos)
  
  write.table(map, "wgt.map", 
              row.names = F, col.names = F, quote = F, sep = " ")
  

  
#workers
  #create data matrix
  wgt = matrix(nrow = nrow(af), 
               ncol = length(unique(samples$queen_id))) %>% 
    as.data.frame
  colnames(wgt) = paste0(unique(samples$queen_id), "_w")
  
  #as numeric
  gt3 = apply(gt2, 2, as.numeric)
  
  #loop through colonies
  ucols = unique(samples$queen_id)
  for(i in 1:length(ucols)){
    #TODO: ensure not too many NAs
    wgt[,i] = rowMeans(gt3[,grepl(ucols[i], colnames(gt3))], na.rm = T)
    
  }  
  
  #center by subtracting 2*f
  wgt.c = apply(wgt, 2, function(x){
    x - (2*af$f)
  })
  
  sum(is.na(wgt.c))
  
  #assume 0 where NA???
    #try imputation??
  wgt.c[is.na(wgt.c)] = 0
  
  #vanraiden scaling parameter
  k = 2 * sum(af$f * (1-af$f))
  
  #G matrix
  G.w = (t(wgt.c) %*% wgt.c) / k
  
  #heatmap(G.w)
  
  #fix symmetry
  for(i in 1:dim(G.w)[1]){
    for(j in 1:i){
      G.w[j,i] = G.w[i,j]
    }
  }
  is.positive.definite(G.w)

  #TODO:
    #address diagonal (not inbreeding?)
    #compare new to old w matrix
    #combine with queens
    #GWAS ... ??? 
  
  
  
