#install.packages('vcfR')
library(tidyverse)
library(vcfR)
library(purrr)

#read files
  #read in sample names
  samples = read.delim("data/header.txt", header = F) %>% t() %>% 
    as.data.frame() %>% 
    rename(sample_id = 1) %>% 
    filter(grepl("23CBH", sample_id)) %>% 
    mutate(queen_id = gsub("_[0-9]*", "", sample_id),
           colony_id = gsub("23CBH", "", queen_id)) 
  
  #select samples to remove
  s5 = read.delim("data/s5summary.txt", sep = "")
  
  s5$bad = (is.nan(s5$nsites) | s5$nsites == 0 |
              s5$reads_consens < 2e4)
  s5$sample_id = rownames(s5)
  s5 = s5 %>% select(sample_id, reads_consens, nsites, bad)
  #add sample info
  s5 = s5 %>% 
    mutate(queen_id = gsub("_[0-9]*", "", sample_id),
           colony_id = gsub("23CBH", "", queen_id)) 
  
  s5.sum = s5 %>% group_by(colony_id) %>% 
    summarise(nbad = sum(bad))
  
  # hist(s5$reads_consens)
  # hist(s5.sum$nbad)
  # 
  # sum(s5$reads_consens[!s5$bad], na.rm = T)
  # 
  # #output for removal
  # s5.remove = as.character(s5$sample_id[s5$bad])
  # write.table(s5.remove, "data/toremove.txt", quote = F,
  #             row.names = F, col.names = F)
  # 
  # #sample threshold for 90% complete data?
  # sum(!s5$bad) * 0.90
  
  #read in VCF
  vcf.raw = read.vcfR("data/23CBH-filter.vcf")
  
  #convert to dosage
    gt = extract.gt(vcf.raw)
    
    gt2 = gt
      gt2[gt == "0/0"] = "0"
      gt2[gt2 == "1/1"] = "2"
      gt2[gt2 == "1/0" | gt2 == "0/1"] = "1"
      gt2 = as.data.frame(gt2)
      
    rm(gt, vcf.raw)
  
  #read in pedigree
  pheno = read.csv("data/CBH_raw_phenotypes.csv") %>% 
    mutate(colony_id = str_pad(colony_id, 3, "left", "0"))
  samples = samples %>% left_join(pheno %>% select(colony_id, breeder))
  
  # #create representative sample
  # set.seed(123)
  # balance = samples %>% group_by(breeder) %>% 
  #   summarise(colony_id = unique(colony_id)) %>% 
  #   slice_sample(n = 10) %>% 
  #   left_join(samples)
  # 
  # #write out for analysis
  # write.table(balance$sample_id, file = "data/balanceSet.txt",
  #             quote = F, row.names = F, col.names = F)
  # 
  # #calculate allele freq for representative sample
  
  #read in allele frequencies ...
  af = read.delim("data/23CBH.frq", header = F) 
    colnames(af) = c("chr", "pos", "f")
    #af = af %>% filter(chr == 11)
  
    
#TODO: fix missing breeders and 16workers from some colonies

#colony diversity  
#####
    
  #TODO: double check sd and dist calculations...
      
  #define distance function
      cd = function(x){
        #cd of each site
        cdm = rep(NA, nrow(x))
        for(l in 1:nrow(x)){
          cdm[l] = mean(dist(t(x[l,])), na.rm = T)
        }
        return(mean(cdm, na.rm = T))
      }
  #create object
    coldiv = data.frame(queen_id = unique(samples$queen_id), sd.gt = NA, 
                        dist.gt = NA)
  #loop through colonies
    for(i in 1:nrow(coldiv)){
      #gather all relevant genotypes
      c = coldiv$queen_id[i]
      cgt = gt2[,grepl(c, names(gt2))]
      #sd of genotype
      coldiv$sd.gt[i] = mean(
                              apply(cgt, 1, function(x){sd(x, na.rm = T)}),
                              na.rm = T)
      #genoytpe distance
      coldiv$dist.gt[i] = cd(cgt)
    }    
    
    hist(coldiv$sd.gt)
    plot(coldiv$sd.gt, coldiv$dist.gt)
      
    #write_csv(coldiv, file = "data/coldiv.csv")
    
    
    
  #quickly check against phenos ... 
    pheno$queen_id = paste0("23CBH", pheno$colony_id)
    pheno = pheno %>% left_join(coldiv)
    sum(!is.na(pheno$dist.gt))
    
    pheno.f = pheno %>% 
      filter(!is.na(dist.gt), !is.na(pattern)) %>% 
      mutate(apiary_id = as.character(apiary_id))
    
    mod = lm(honey ~ apiary_id + established + breeder,  data = pheno.f)
    
    test = data.frame(resid = mod$residuals, dist = pheno.f$dist.gt)
    
    #randomize for fun
    #test$resid = rnorm(nrow(test), 0, 1)
    
    plot(test$dist, test$resid)
    
    #bins
    se = function(x) sd(x) / sqrt(length(x))
    test = test %>% 
      mutate(bins = cut(dist, 
                        breaks = seq(min(dist) - 0.01, max(dist), length.out = 12))) 
    
    plot(test$bins, test$resid)
    
    test.sum = test %>% group_by(bins) %>% 
    summarise(mean = mean(resid),
              sd = sd(resid),
              n = n()) %>% 
      filter(n > 10)
    
    ggplot(test.sum, aes(x = bins, y = sd)) + geom_point()
    


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
#write object 
  #write queen for plink .ped
    pt1 = data.frame(IID = names(qgt)) %>% 
      left_join(samples %>% select(IID = queen_id, FID = breeder), multiple = 'first') %>% 
      mutate(FID = gsub("-", "", FID)) %>% 
      select(FID, IID) %>% 
      mutate(father = 0, mother = 0, sex = 2, pheno = 0)
    
    pt2 = t(qgt)
      pt2[pt2 == 0] = "TT"
      pt2[pt2 == 1] = "AT"
      pt2[pt2 == 2] = "AA"
      pt2[is.na(pt2)] = "00"
      
    
    qgt.out = cbind(pt1, pt2)
    
    write.table(qgt.out, "data/qgt.ped", 
                row.names = F, col.names = F, quote = F, sep = " ")
    #TODO: write .map file
    
  #TODO: write worker
    #??? maybe beagle format?

#calculate the A matrix
  #this should be checked against the pedigree!
  # 




#Theoretical Sim
#####


#TODO: remove sites with narrow probabilities 

#using joint likelihood of het and hom proportions of workers

#simulation: error of real genotypes
lsim = function(p, qg){
  if(qg == 2) rqg = c(1,1)
  if(qg == 1) rqg = c(0,1)
  if(qg == 0) rqg = c(0,0)
  
  
  # rqg = c(1,0)
  # p = 0.8
  reps = 500
  nworker = 8
  #queen genotype likelihood
  qgl = data.frame("w2" = c(p, p/2, 0), 
                   "w1" = c(1-p, 1/2, p), 
                   "w0" = c(0, (1-p)/2, 1-p))
  rownames(qgl) = c("q2", "q1", "q0")
  
  #simulate workers
  call = rep(NA, reps)
  for(i in 1:reps){
    swg = sample(rqg, nworker, replace = T) + 
      sample(c(1,0), nworker, prob = c(p, 1-p), replace = T)
    
    #w2, w1, w0
    swg.counts = c(sum(swg == 2), sum(swg == 1), sum(swg == 0))
    
    #multinomial
      #lik of q2
      lq2 = dmultinom(swg.counts, nworker, as.numeric(qgl[1,]))
      #likelihood of q1
      lq1 = dmultinom(swg.counts, nworker, as.numeric(qgl[2,]))
      #lik of q0
      lq0 = dmultinom(swg.counts, nworker, as.numeric(qgl[3,]))
      
      c(lq2, lq1, lq0)

    call[i] = c("q2", "q1", "q0")[which.max(c(lq2, lq1, lq0))]
  }
  
  if(qg == 0) return(sum(call != "q0") / length(call))
  if(qg == 1) return(sum(call != "q1") / length(call))
  if(qg == 2) return(sum(call != "q2") / length(call))
} 
  
#run for a range of p and each genotype
s = seq(0, 1, 0.1)
plotdb = data.frame(
  p = s,
  q0 = s|> map_vec(lsim, qg = 0),
  q1 = s|> map_vec(lsim, qg = 1),
  q2 = s|> map_vec(lsim, qg = 2)
)

#pivot and plot
plotdb.long = plotdb %>% 
  pivot_longer(cols = starts_with("q"), names_to = "realGT")

ggplot(plotdb.long) + 
  geom_line(aes(x = p, y= value, color = realGT), size = 2) + 
  theme_bw() + 
  labs(y = "error rate", x = "pop'n allele freq", color = "real queen\ngenotype")







#simulation2: error of called genotypes given p
  #i think this only holds if there are equal probablilities of 
  #q0,q1,q2, which is not a good assumption
  #I could change the simulation to also choose a queen gt based on p
  #but is that a good assumption?

lsim2 = function(p, qg){
  if(qg == 2) rqg = c(1,1)
  if(qg == 1) rqg = c(0,1)
  if(qg == 0) rqg = c(0,0)
  
  
  # rqg = c(1,0)
  # p = 0.8
  reps = 2000
  nworker = 8
  #queen genotype likelihood
  qgl = data.frame("w2" = c(p, p/2, 0), 
                   "w1" = c(1-p, 1/2, p), 
                   "w0" = c(0, (1-p)/2, 1-p))
  rownames(qgl) = c("q2", "q1", "q0")
  
  #simulate workers
  call = rep(NA, reps)
  for(i in 1:reps){
    swg = sample(rqg, nworker, replace = T) + 
      sample(c(1,0), nworker, prob = c(p, 1-p), replace = T)
    
    #w2, w1, w0
    swg.counts = c(sum(swg == 2), sum(swg == 1), sum(swg == 0))
    
    #multinomial
    #lik of q2
    lq2 = dmultinom(swg.counts, nworker, as.numeric(qgl[1,]))
    #likelihood of q1
    lq1 = dmultinom(swg.counts, nworker, as.numeric(qgl[2,]))
    #lik of q0
    lq0 = dmultinom(swg.counts, nworker, as.numeric(qgl[3,]))
    
    #c(lq2, lq1, lq0)
    
    call[i] = c("q2", "q1", "q0")[which.max(c(lq2, lq1, lq0))]
  }
  
  out = data.frame(q2 = sum(call == "q2"),
                   q1 = sum(call == "q1"),
                   q0 = sum(call =="q0"))
  return(out)
} 

#run for a range of p and each genotype
s = seq(0, 1, 0.1)

plotdb2 = rbind(s |> map_dfr(lsim2, qg = 2),
                s |> map_dfr(lsim2, qg = 1),
                s |> map_dfr(lsim2, qg = 0))
  plotdb2$p = rep(s, 3)
  plotdb2$rqg = c(rep(2, length(s)), rep(1, length(s)), rep(0, length(s)))

#pivot and plot
plotdb2.long = plotdb2 %>% 
  pivot_longer(cols = starts_with("q"), names_to = "calledGT") %>% 
  mutate(rqg = as.character(rqg))
  plotdb2.long$calledGT = gsub("q", "", plotdb2.long$calledGT)
  
  plotdb2.long$error = (plotdb2.long$rqg != plotdb2.long$calledGT)

  sim2 = plotdb2.long %>% group_by(p, calledGT, error) %>% 
    summarise(count = sum(value)) %>% 
    pivot_wider(names_from = error, values_from = count) %>% 
    mutate(error_rate = `TRUE` / (`FALSE` + `TRUE`))
  
  ggplot(sim2, aes(x = p, y = error_rate, color = calledGT)) +
    geom_line(size = 2) + 
    theme_bw()
  


#likelihood cutoff will work, but will bias towards homozygous sites and 
  #intermediate allele freqs

#phasing will help, but complicated with GBS data
  #pedigreee is good, which will help with phasing 
  #shapeit2 for GBS
  #


#using total allele count 
  #... is less accurate...
# lsim2 = function(p, qg){
#   if(qg == 1) rqg = c(0,1)
#   if(qg == 0) rqg = c(0,0)
#   if(qg == 2) rqg = c(1,1)
#   
#   reps = 1000
#   nworker = 8
#   #expeted worker freq
#   ewf = data.frame("q0" = p/2, 
#                    "q1" = (p+0.5)/2, 
#                    "q2" = (p+1)/2)
#   
#   #simulate workers
#   call = rep(NA, reps)
#   for(i in 1:reps){
#     #simulate worker geno
#     swg = sample(rqg, nworker, replace = T) + 
#       sample(c(1,0), nworker, prob = c(p, 1-p), replace = T)
#     #total alleles
#     sswg = sum(swg)
#     
#     #total likelihood of q0
#     lq0 = dbinom(sswg, 2*nworker, ewf[1, 'q0'])
#     
#     #total likelihood of q1
#     lq1 = dbinom(sswg, 2*nworker, ewf[1, 'q1'])
#     
#     #total likelihood of q2
#     lq2 = dbinom(sswg, 2*nworker, ewf[1, 'q2'])
#     
#     call[i] = c("q0", "q1", "q2")[which.max(c(lq0, lq1, lq2))]
#   }
#   
#   if(qg == 0) return(sum(call != "q0") / length(call))
#   if(qg == 1) return(sum(call != "q1") / length(call))
#   if(qg == 2) return(sum(call != "q2") / length(call))
# } 
# 
# s = seq(0, 1, 0.1)
# plotdb2 = data.frame(
#   p = s,
#   q0 = s|> map_vec(lsim2, qg = 0),
#   q1 = s|> map_vec(lsim2, qg = 1),
#   q2 = s|> map_vec(lsim2, qg = 2)
# )
# 
# ggplot(data = plotdb2) + 
#   geom_line(aes(x = p, y = q0), color = 'orange') + 
#   geom_line(aes(x = p, y = q1), color = 'blue') + 
#   geom_line(aes(x = p, y = q2), color = 'red') +
#   theme_bw()