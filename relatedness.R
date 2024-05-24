library(tidyverse)
library(Matrix)
library(qgraph)
library(reshape2)


#plink rel
# rel = read.delim("data/test.relatedness",
#                  header = T)
#   colnames(rel) = c("id1", "id2", "r")
# 
# 
# relid = rel$id2[rel$id1 == rel$id1[1]]
# 
# 
# relmat = matrix(NA, nrow = length(relid), ncol = length(relid))
# 
# colnames(relmat) = rownames(relmat) = relid
# 
# relmat[lower.tri(relmat,diag=T)] <- rel$r
# 
# relmat = forceSymmetric(relmat,uplo="L") %>% as.matrix


#ngs rel
reladmix = read.delim("data/bag13.rel", header = T)

reladmix$r = reladmix$k2 + (reladmix$k1 / 2)

rel = reladmix %>% select(id1 = ind1, id2 = ind2, r)

relid = read.delim("data/admix.fam",
                   header = F, sep = " ") %>% select(id = V1) %>%
  filter(!grepl("SRR", id))


relmat = matrix(NA, nrow = nrow(relid), ncol = nrow(relid))
colnames(relmat) = rownames(relmat) = relid$id
relmat[lower.tri(relmat,diag=F)] <- rel$r
diag(relmat) = 1
relmat = forceSymmetric(relmat,uplo="L") %>% as.matrix



#check a's and b's (replicates)
relmat["23-T5w01a", "23-T5w01b"]

meta = data.frame(id = relid$id)
meta$col = gsub("23-(.*)[wd].*", "\\1", meta$id)


#helpful labels
new.names = matrix( data =
                      c("II17", "Ii1a",
                        "II18", "Ii1b",
                        "II19", "Tx1",
                        "II20", "Ix1",
                        "II22", "Ti1",
                        "II28", "Ix2",
                        "II29", "Ti2",
                        "II32", "Ii2",
                        "II40", "Ix3",
                        "II42", "Ti3",
                        "II43", "Tx3",
                        "II44", "Ii3"),
                    nrow = 12, byrow = T
)

meta$new = meta$col
for (i in 1:nrow(new.names)){
  meta$new[grepl(new.names[i,1], meta$id)] = new.names[i,2]
}

meta$cross = paste0(meta$col, ":", meta$new)



qgraph(relmat, layout='spring', vsize=5,
       cut = 0, repulsion = 1, diag = F,
       minimum = 0.01, maximum = 0.75,
       labels = meta$col,
       #label.cex = 1.2,
       groups = meta$cross, palette = "pastel",
       #edge.labels = T,
       legend = T)

#plot against pedigree relatedness
write.table(meta %>% select(id, col), "data/ids.txt", col.names = F, row.names = F, quote = F)

GRM = relmat

PRM = read.delim("/home/dylan/Documents/bees/harpurlab/misc/ped/divExPed.tab",
                 header = T, sep = "")
colnames(PRM) = rownames(PRM)


PRM2 = PRM[colnames(GRM), colnames(GRM)]

d = GRM - PRM2

d.t = lower.tri(d, diag = F) %>% as.vector()

colnames(misses) = "value"

ggplot(data = misses, aes(x = value)) + geom_histogram()


ggplot(misses, aes(x = ))

offset = rowSums(abs(PRM2 - GRM))

off2 = cbind(meta, offset)

off2 %>% group_by(col) %>% summarise(md = mean(offset))








#plot admixture
#####

#read in reference data
reffam = read.delim("data/refData.txt")
refadmix = cbind( read.delim("data/admix.4.Q", header = F, sep =""),
                  read.delim("data/admix.fam", sep = "", header = F) %>% select(oldid = V1)) %>%
  left_join(reffam %>% select(lineage, oldid = SRR))

refadmix.melt = refadmix %>% melt(id.vars = c("oldid", "lineage"),
                                  measure.vars = c("V1", "V2", "V3", "V4"),
                                  variable.name = "fam")



#references to find lineage assignments
ggplot(data = refadmix.melt) +
  geom_bar(aes(x = oldid, y = value, fill = fam), 
           stat='identity', width = 1) +
  facet_grid(cols = vars(lineage), scales = "free_x") +
  scale_fill_brewer(palette = "PRGn") +
  labs(x = "sample", y = "Proportion of Genome", fill = "Lineage",
       title = "maf 0.05") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())


#just the samples
admix = refadmix %>% filter(!grepl("SRR", oldid)) %>%
  rename(id = oldid) %>%
  rename(M = V1, C = V2, O = V3, A = V4) %>%
  mutate(col = gsub("23-(.*)[wd].*", "\\1", id))

admix.melt = admix %>% melt(id.vars = c("id", "col"),
                            measure.vars = c("A", "M", "C", "O"),
                            variable.name = "fam")

ggplot(data = admix.melt) +
  geom_bar(aes(x = id, y = value, fill = fam), 
           stat='identity', width = 1) +
  facet_grid(cols = vars(col), scales = "free_x") +
  scale_fill_brewer(palette = "PRGn") +
  labs(x = "sample", y = "Proportion of Genome", fill = "Lineage") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank())


