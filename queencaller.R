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
af.raw = read.delim("23CBH-updated.frq", header = F) 
  colnames(af.raw) = c("chr", "pos", "f")
  #af = af %>% filter(chr == chrno)
  
  #remove zeros
  sum(af.raw$f == 0)
  
  af = af.raw[af.raw$f != 0,]


  #read in genotypes
  #vcf.raw = read.vcfR(paste0("chrs/chr",chrno,"/23CBH-filter-", chrno, ".vcf"))
  vcf.raw = read.vcfR("23CBH-updated-filter.vcf")
  
  
  #convert to dosage
  gt = extract.gt(vcf.raw)
  
  gt2 = gt
  gt2[gt == "0/0"] = "0"
  gt2[gt2 == "1/1"] = "2"
  gt2[gt2 == "1/0" | gt2 == "0/1"] = "1"
  gt2 = as.data.frame(gt2)
  
  rm(gt, vcf.raw)
  
  #remove zeros
    gt2 = gt2[af.raw$f != 0,]
  
  
  
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
      
      #TODO: use pedigree information?
      
      #pull genotypes
      cworkers.gt = as.numeric(gt2.chr[i, grepl(col, colnames(gt2.chr))])
      #remove NAs and count
      cworkers.gt = cworkers.gt[!is.na(cworkers.gt)]
      nworker = length(cworkers.gt)
      #require 4 non-missing
      if(nworker < 5){
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
 #TODO: average inbreeding???
    
       
#combine with queens
    
    qwgt = cbind(finalqgt, wgt)
    
    #switch to ref allele instead of alt
    qwgt.ref = 2-qwgt
    
    qwgt.out = qwgt.ref
      qwgt.out[is.na(qwgt.out)] = 5
      qwgt.out = round(qwgt.out,4)
      
        
    # #write out simplified form
    # write.table(t(qwgt.out), file = "queenworker.geno",
    #             sep = " ", col.names = F, quote = F)
    
    #different format to combine wtih poolseq data
      #TODO: ensure these are the correct positions!
    gencomb = cbind(af$chr, af$pos, qwgt.out)
      names(gencomb)[1:2] = c("chr", "pos")
    
    write.table(gencomb, file = "23CBH_qw_ref.geno",
                sep = " ", col.names = T, quote = F, row.names = F)
    
    
    
#     
#     #rm(finalqgt, wgt)
# 
# 
# #write queen object 
#   #write queen for plink .ped
#   pt1 = data.frame(FID = names(finalqgt), IID = names(finalqgt)) %>% 
#     mutate(father = 0, mother = 0, sex = 2, pheno = 0)
#   
#   pt2 = t(finalqgt)
#   pt2[pt2 == 0] = "TT"
#   pt2[pt2 == 1] = "AT"
#   pt2[pt2 == 2] = "AA"
#   pt2[is.na(pt2)] = "00"
#   
#   qgt.out = cbind(pt1, pt2)
#   
#   write.table(qgt.out, "qgt.ped", 
#               row.names = F, col.names = F, quote = F, sep = " ")
#   #write .map file
#   map = af %>% 
#     mutate(variant = paste0(chr, ":", pos),
#            cm = 0) %>% 
#     select(chr, variant, cm, pos)
#   
#   write.table(map, "qgt.map", 
#               row.names = F, col.names = F, quote = F, sep = " ")
#   
# #write raw worker object
#   #write worker for plink .ped
#   pt1 = data.frame(FID = colnames(gt2), IID = colnames(gt2)) %>% 
#     mutate(father = 0, mother = 0, sex = 2, pheno = 0)
#   
#   pt2 = t(gt2)
#   pt2[pt2 == 0] = "TT"
#   pt2[pt2 == 1] = "AT"
#   pt2[pt2 == 2] = "AA"
#   pt2[is.na(pt2)] = "00"
#   
#   wgt.out = cbind(pt1, pt2)
#   
#   write.table(wgt.out, "wgt.ped", 
#               row.names = F, col.names = F, quote = F, sep = " ")
#   #write .map file
#   map = af %>% 
#     mutate(variant = paste0(chr, ":", pos),
#            cm = 0) %>% 
#     select(chr, variant, cm, pos)
#   
#   write.table(map, "wgt.map", 
#               row.names = F, col.names = F, quote = F, sep = " ")
#   
# 
# 
# #calc GRM
#   #remove missing > 50%
#   missing.count = colSums(is.na(qwgt))
#   
#   sum(missing.count > .5 * nrow(af))
#   
#   names.remove = names(missing.count)[missing.count > .5 * nrow(af)]
#   names.remove = c(names.remove, paste0(names.remove, "_w"))
#   
#   qwgt.filter = qwgt[,!colnames(qwgt) %in% names.remove]
#   
#   
#   #center by subtracting 2*f
#   qwgt.c = apply(qwgt.filter, 2, function(x){
#     x - (2*af$f)
#   })
#   
#   
#   #assume 0 where NA???
#     #try imputation??
#   qwgt.c[is.na(qwgt.c)] = 0
#   
#   #vanraiden scaling parameter
#   k = 2 * sum(af$f * (1-af$f))
#   
#   #G matrix
#   G.qw = (t(qwgt.c) %*% qwgt.c) / k
#   
#   G.qw = round(G.qw, 4)
#   
#   #heatmap(G.w)
#   
#   #fix symmetry
#   for(i in 1:dim(G.qw)[1]){
#     for(j in 1:i){
#       G.qw[j,i] = G.qw[i,j]
#     }
#   }
#   is.positive.definite(G.qw)
#   
#   write.table(G.qw, "queenworkerGRM.txt", sep = " ",
#               quote = F)
# 
#   
  
  
