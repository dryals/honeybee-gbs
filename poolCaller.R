library(tidyverse)
library(doParallel)
library(matrixcalc)
library(AGHmatrix)

#Dylan Ryals 25 AUG 2025

#last edit: 26 DEC 2025
  
  #re-doing with AlphaPeel allele frequencies

#TODO:
  #edit to accept a filename argument or hardcode a filename...
  #doulbe check REF and ALT alleles throughout pipelines

setwd("/scratch/negishi/dryals/gbs/24CBH/analysis")

sync = read.delim("chrPar/24CBH-ap.sync", sep = "", header = F)
#sync = read.delim("data/test4.sync", sep = "", header = F)

#sample names
sampnames = read.delim("24CBHpool.bamlist", header = F)
  sampnames = gsub(".*refmapping/([0-9A-Za-z]*)[-].*", "\\1", sampnames$V1)

  #count total depth, ignoring N and indels
  sync.g = sync[,-(1:3)]
  
  countdepth = function(x){
    # x = sync.g[,1]
    x.tab = (str_split(x, ":", simplify = T))
    x.tab = apply(x.tab, 2, as.numeric)
    #remove N and dels
    x.tab = x.tab[,1:4]
    return(rowSums(x.tab))
  }
  
  samples.depth = apply(sync.g, 2, countdepth)
  
  #count ref observations
  sync.r = sync[,3] %>% toupper()
  
  countref = function(x){
    #x = sync.g[,1]
    x.tab = (str_split(x, ":", simplify = T))
    x.tab = apply(x.tab, 2, as.numeric)
    #remove N and dels
    x.tab = x.tab[,1:4]
    #count ref
    colnames(x.tab) = c("A", "T", "C", "G")
    x.ref = vector(mode = "numeric", length = nrow(x.tab))
    for(i in 1:nrow(x.tab)){
      x.ref[i] = x.tab[i,sync.r[i]]
    }
    return(x.ref)
  }

  samples.ref = apply(sync.g, 2, countref) 
  
  sync.c = sync[,1]


#call queen genotypes
  # queens.geno = samples.depth
  #   queens.geno[,] = NA
# 
#     #MLE of allele freq at site, VERY SLOW and hardly different
#     MLE = function(l){
#           allf = data.frame(freq = seq(0, 1, 0.01),
#                             lik = NA)
#           for(i in 1:nrow(allf)){
# 
#             f = allf$freq[i]
#             obs = data.frame(d = samples.depth[l,],
#                              x = samples.ref[l,])
# 
#             allf$lik[i] = prod(apply(obs, 1, function(x){dbinom(x[2], x[1], f)}))
# 
#           }
#           return(allf$freq[which.max(allf$lik)])
#     }
#     frq0 = sapply(1:100, MLE)
#     
#     plot(frq[1:100], frq0)
# 
#     frq = rowSums(samples.ref) / rowSums(samples.depth)
    
  frq = read.delim()
    
    allchrs = unique(sync.c)
    
  #use max prob to call queen genotype
    #TODO: write NA if too few reads!!!
  callQueenPool = function(CHR){
    #just sites within CHR
    samples.depth.chr = samples.depth[sync.c==CHR,]
    samples.ref.chr = samples.ref[sync.c==CHR,]
    frq.chr = frq[sync.c==CHR]
    queens.geno.chr = samples.ref.chr
      queens.geno.chr[,] = NA
    #loop sites
    for(l in 1:nrow(queens.geno.chr)){
      #loop colonies 
      for(c in 1:ncol(queens.geno.chr)){
        #return NA if too few observations
        if(samples.depth.chr[l,c] < 6){
          queens.geno.chr[l,c] = NA
        } else{
          #ref allele freq at site
          f = frq.chr[l]
          #probability of each potential queen genotype
                #prob q is gt  *prob worker obs given q is gt and f
          pg2 = f**2 *          dbinom(samples.ref.chr[l,c], samples.depth.chr[l,c], (f + 1)/2)
          pg1 = 2*(f * (1-f)) * dbinom(samples.ref.chr[l,c], samples.depth.chr[l,c], (f + 0.5)/2)
          pg0 = (1-f)**2 *      dbinom(samples.ref.chr[l,c], samples.depth.chr[l,c],  f/2)
          
          queens.geno.chr[l,c] = c(2,1,0)[which.max(c(pg2, pg1, pg0))]
        }
      }
    }
    return(cbind(sync[sync.c==CHR,1:2],queens.geno.chr))
  }
  
  registerDoParallel(cores = 16)
  
  qgt.out = foreach(chr=allchrs, .combine=rbind) %dopar%
    callQueenPool(chr)
  
  colnames(qgt.out) = c("chr", "pos", sampnames)
  
  #how much missing?
    sum(is.na(qgt.out))
    qgt.m = qgt.out[,-(1:2)]
    hist(colSums(is.na(qgt.m)))
    
  #remove samps with >20% missing sites
    keep = colSums(is.na(qgt.m)) < (.2 * nrow(qgt.out))
    #does this code actualy work?
    finalqgt = cbind(qgt.out[,1:2],
                     qgt.m[,keep])
  
  #workers
    workers.geno = round(2 * samples.ref / samples.depth, 4)
    colnames(workers.geno) = paste0(sampnames, "_w")
    #remove same samples as with queens...
    workers.geno = workers.geno[,keep]
  
    dim(finalqgt)
    dim(workers.geno)
  
  #combine workers and queens, rename with _w...
  qwgt24 = cbind(finalqgt, workers.geno)
  
  
#combine with 2023 genotypes .... 
  
    #read in 2023 genotypes
  qwgt23 = read.delim("../../23CBH/analysis/23CBH_qw_ref.geno", sep = " ")
    colnames(qwgt23) = gsub("X", "", colnames(qwgt23))
    
    #verify worker genotypes are float
    hist(as.matrix(qwgt23[,grepl("_w", colnames(qwgt23))]))
    
    #ensure both years use REF or ALT
      #ensure sites are identical
    chrrename = read.delim("chrsrename.txt", header = F, sep = "")
      colnames(chrrename) = c("chr", "long")
      
    qwgt23 = qwgt23 %>% left_join(chrrename) %>% 
      select(chr = long, pos, starts_with("23"))
    
    qwgt23[qwgt23 == 5] = NA
    
    #combine into massive matrix
    
    qwgt.c = qwgt23 %>% left_join(qwgt24, by = c("chr", "pos"))
    
    dim(qwgt.c)
    
    2 * (sum(!grepl("_w", colnames(qwgt23))) + 
           sum(!grepl("_w", colnames(qwgt24))) - 4)
    
    #TODO: compare allele freqs between years
    
      #search for and resolve duplicated queens/colonies
        queennames = data.frame(full = colnames(qwgt.c)[-(1:2)])
          queennames = queennames %>% 
            filter(!grepl("_w", full)) %>% 
            mutate(short = gsub("^2[34]", "", full)) %>% 
            group_by(short) %>% 
            summarise(n = n())
          
          #065, 193, 221, II96 are duplicated ...
          #leave in for now, run GRM, and see if they're related!
  
#calculate GRM
#TODO: double check this method!
# try loading genotypes into PLINK and using KING est?
          
    qwgt.grm = qwgt.c[,-(1:2)]
    
    
    #remove missing > 50%
    missing.count = colSums(is.na(qwgt.grm))
    
    sum(missing.count > .5 * nrow(qwgt.grm))
    
    names.remove = names(missing.count)[missing.count > .5 * nrow(qwgt.grm)]
    #add some problematic samples
    names.remove = c(names.remove, "23CBH001", "23CBH246", "23CBH266")
    
    names.remove = c(names.remove, paste0(names.remove, "_w"))
    
    qwgt.grm.filter = qwgt.grm[,!colnames(qwgt.grm) %in% names.remove]
    
    #verify worker genotypes are float
    hist(as.matrix(qwgt.grm.filter[,grepl("_w", colnames(qwgt.grm.filter))]))
    hist(as.matrix(qwgt.grm.filter[,grepl("^23.*_w", colnames(qwgt.grm.filter))]))
    hist(as.matrix(qwgt.grm.filter[,grepl("^24.*_w", colnames(qwgt.grm.filter))]))
    
    #folloing vanRaden 2008
    
    #test queens only
    #qwgt.grm.filter = qwgt.grm.filter[,!grepl("_w", colnames(qwgt.grm.filter))]
    
    M = (2 - qwgt.grm.filter) -1 # +1 carries two ALT alleles
    
    M = t(M)
    
    #minor alleles in base popn
    P = read.delim("/scratch/negishi/dryals/gbs/23CBH/analysis/23CBH-updated.frq",
                   header = F)
      #this counts alternate allele
      colnames(P) = c("chr", "pos", "f")
      #match to sites 
      P = P %>% left_join(chrrename) %>% 
        select(chr = long, pos, f)
      P = qwgt.c %>% select(chr, pos) %>% 
        left_join(P)
      #in case we want MINOR allele freq ... ???
      #P$f[P$f > 0.5] = 1 - P$f[P$f > 0.5]
      
      P.c = 2*(P$f - 0.5)
      
      #for each individual, adjust genotypes by frequencies
      Z = apply(M, 1, function(x){
        as.numeric(x - P.c)
      })
      
      #transpose 
      Z = t(Z)
      
      hist(colMeans(Z, na.rm = T))
      #hist(colMeans(Z, na.rm = T))
      
      #impute missing values
      impute = function(x){
        x[is.na(x)] = mean(x, na.rm = T)
        return(x)
      }
      
      Z = apply(Z, 2, impute)
      
      G = Z %*% t(Z) /
        (2 * sum(P$f * (1-P$f)))

    # #vanraiden scaling parameter
    # k = 2 * sum(af * (1-af))
    # 
    # #G matrix
    # G.qw = (t(qwgt.grm.c) %*% qwgt.grm.c) / k
    # 
    # G.qw = round(G.qw, 4)
    # 
    # #heatmap(G.w)
    # 
    # #fix symmetry
    # for(i in 1:dim(G.qw)[1]){
    #   for(j in 1:i){
    #     G.qw[j,i] = G.qw[i,j]
    #   }
    # }
    # is.positive.definite(G.qw)

    write.table(G, "24CBHqwGRM.txt", sep = " ",
                quote = F)


#testing and visualizing
    hist(diag(G))
    hist(diag(G)[grepl("_w", colnames(G))])
    hist(diag(G)[!grepl("_w", colnames(G))])
    
    G.q = G[!grepl("_w", colnames(G)), !grepl("_w", colnames(G))]
    
    heatmap(G.q)
    
    hist(diag(G.q))
    #which(diag(G.q) > 1.8)
    
    colnames(G.q)[grepl(".x", colnames(G.q))]
    #065, 193, 221, II96
    #"23CBH065.x"  "23CBH221.x"  "23CBHII96.x"
    
    #well that's not a great sign ...
    G.q["23CBH065.x", "23CBH065.y"]
    G.q["23CBH221.x", "23CBH221.y"]
    #G.q["23CBHII96.x", "23CBHII96.y"]
    
    #load pedigree
    ped = read.csv("fullpedigree.csv")
    
    testgroup = ped$queen_id[ped$mother%in% c("23_065")]
      testgroup = testgroup[!grepl("p", testgroup)]
      testgroup = gsub("_", "CBH", testgroup)
      testgroup = testgroup[testgroup %in% colnames(G.q)]
      
      G.test = G.q[testgroup, testgroup]
      
      hist(G.test[upper.tri(G.test, diag = F)])
      mean(G.test[upper.tri(G.test, diag = F)])
      
      hist(G.q[upper.tri(G.q, diag = F)])
      mean(G.q[upper.tri(G.q, diag = F)])
      
      heatmap(G.test)
    
  
#write out...
  

    
  
  
  
  
  
    
    #OLD pileup code
    # #read and format pileup results
    #   pileup = read.delim("data/pileuptest4.out", header = F)
    #   
    #   #pull samples
    #   samples = pileup[,-(1:3)]
    #   locations = pileup[,(1:2)]
    #   
    #   #cat first two entries, dump third
    #   samples2 = samples
    #     for(i in 1:(ncol(samples)/3)){
    #   
    #       samples2[,i] = paste(samples2[,i], samples2[,(i+1)])
    #       samples2 = samples2[,-((i+1):(i+2))]
    #   
    #     }
    #   
    #   #total depth
    #   samples.depth = apply(samples2, 2, 
    #                         function(x){as.numeric(gsub("([0-9]*).*", "\\1", x))})
    #   
    #   #counts matching ref
    #   samples.ref = apply(samples2, 2, 
    #                       function(x){as.numeric(str_count(x, "[.,]"))})
    #   
    #   #ref allele frequencies (for worker grm)
    #   samples.frq = round(
    #       apply(samples.ref, 2, as.numeric) / apply(samples.depth, 2, as.numeric),
    #                       4)
    #     #annoying, but should leave as NA for now
    #     #samples.frq[is.nan(samples.frq)] = 0
    #   
    #   # #add location info
    #   # samples.depth = cbind(locations, samples.depth)
    #   # samples.ref = cbind(locations, samples.ref)

#old version, beethoven format


# #read
# pileup = read.delim("data/pileuptest3.out", header = F)
# 
# #pull samples
# samples = pileup[,-(1:3)]
# locations = pileup[,(1:2)]
#   
# #cat first two entries, dump third
# samples2 = samples
#   for(i in 1:(ncol(samples)/3)){
#     
#     samples2[,i] = paste(samples2[,i], samples2[,(i+1)])
#     samples2 = samples2[,-((i+1):(i+2))]
#     
#   }
# 
# #total depth
# samples.depth = apply(samples2, 2, function(x){gsub("([0-9]*).*", "\\1", x)})
# 
# #counts matching ref
# samples.ref = apply(samples2, 2, function(x){str_count(x, "[.,]")})
# 
# #ref allele frequencies (for worker grm)
# samples.frq = round(
#     apply(samples.ref, 2, as.numeric) / apply(samples.depth, 2, as.numeric), 
#                     4)
#   #annoying, but should leave as NA for now
#   #samples.frq[is.nan(samples.frq)] = 0
# 
# #add location info 
# samples.depth = cbind(locations, samples.depth)
# samples.ref = cbind(locations, samples.ref)
# 
# #write out
# write.table(file = "data/tmp.depth", samples.depth,
#             quote = F, col.names = F, row.names = F)
# 
# write.table(file = "data/tmp.ref", samples.ref,
#             quote = F, col.names = F, row.names = F)
# 
# write.table(file = "data/tmp.frq", samples.frq,
#             quote = F, col.names = F, row.names = F)
