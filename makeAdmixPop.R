library(dplyr)

reffam = read.table("data/refData.txt", header = T, sep = "\t")
  id = read.delim("/scratch/negishi/dryals/gbs/analysis/admix.fam", sep = " ", header = F) %>% 
    select(1)
  colnames(id) = "id"

sup = id %>% left_join(reffam %>% select(id = SRR, lineage))
sup$lineage[is.na(sup$lineage)] = "-"

write.table(sup$lineage, file = "/scratch/negishi/dryals/gbs/analysis/admix.pop",
            sep = "\t",
            quote = F, row.names = F, col.names = F)
