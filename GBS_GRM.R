#install.packages('vcfR')
library(tidyverse)
library(vcfR)
library(purrr)

#read files
  #read in VCF
  
  #read in sample names
  
  #read in pedigree


#create object for queen and worker group gt

#loop through sites
  #calculate allele freq
  #loop through colonies
    #pull workers
    #estimate queen gt and error
      #multinomial max likelihood
    #calculate average worker gt

#write object 
  #???

#calculate A matrix
  #???




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



  
#read in test data
bag13 = 1
  
#calculate queen GRM
  #calculate queen genotype
    #split into sister groups
    #


#calculate worker group GRM

