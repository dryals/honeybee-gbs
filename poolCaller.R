library(dplyr)
library(readxl)

#read
pileup = read.delim("data/pileuptest.out", header = F)

#pull samples
samples = pileup[,-(1:3)]
locations = pileup[,(1:2)]
  
#cat first two entries, dump third
samples2 = samples
  for(i in 1:(ncol(samples)/3)){
    
    samples2[,i] = paste(samples2[,i], samples2[,(i+1)])
    samples2 = samples2[,-((i+1):(i+2))]
    
  }

samples.depth = apply(samples2, 2, function(x){gsub("([0-9]*).*", "\\1", x)})
samples.ref = apply(samples2, 2, function(x){str_count(x, "[.,]")})

samples.depth = cbind(locations, samples.depth)
samples.ref = cbind(locations, samples.ref)

write.table(file = "data/tmp.depth", samples.depth,
            quote = F, col.names = F, row.names = F)

write.table(file = "data/tmp.ref", samples.ref,
            quote = F, col.names = F, row.names = F)
