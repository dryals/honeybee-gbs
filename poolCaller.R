library(tidyverse)

#Dylan Ryals 25 AUG 2025

#laste edit: 26 AUG 2025

#TODO:
  #edit to accept a filename argument or hardcode a filename...
  #use parallel processing and/or save intermediate files

setwd("/scratch/negishi/dryals/gbs/24CBH/analysis")

sync = read.delim("24CBH.sync", sep = "", header = F)

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


#call queen genotypes
  queens.geno = samples.depth
    queens.geno[,] = NA
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

    frq = rowSums(samples.ref) / rowSums(samples.depth)
    
    #use max prob to call queen genotype
      #probably faster as apply functions
      #could re-write to run chrs in parallel
    for(l in 1:nrow(queens.geno)){
      for(c in 1:ncol(queens.geno)){
        f = frq[l]
        #probability 2 ref alleles
              #prob q is gt  #prob worker obs given q is gt and f
        pg2 = f**2 *          dbinom(samples.ref[l,c], samples.depth[l,c], (f + 1)/2)
        pg1 = 2*(f * (1-f)) * dbinom(samples.ref[l,c], samples.depth[l,c], (f + 0.5)/2)
        pg0 = (1-f)**2 *      dbinom(samples.ref[l,c], samples.depth[l,c],  f/2)
        
        queens.geno[l,c] = c(2,1,0)[which.max(c(pg2, pg1, pg0))]
      }
    }
    
  workers.geno = round(samples.ref / samples.depth, 4)
  
  
#combine with 2023 genotypes .... 
  
    #read in 2023 genotypes
  
    #ensure both years use REF or ALT; ensure sites are identical 
  
    #search for and resolve duplicated individuals
  
    #combine into massive matrix
  
#calculate GRM
  
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
