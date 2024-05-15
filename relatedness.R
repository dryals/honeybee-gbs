library(dplyr)
library(Matrix)
library(qgraph)



rel = read.delim("data/test.relatedness",
                 header = T)
  colnames(rel) = c("id1", "id2", "r")


relid = rel$id2[rel$id1 == rel$id1[1]]


relmat = matrix(NA, nrow = length(relid), ncol = length(relid))

colnames(relmat) = rownames(relmat) = relid

relmat[lower.tri(relmat,diag=T)] <- rel$r

relmat = forceSymmetric(relmat,uplo="L") %>% as.matrix

meta = data.frame(id = relid)
meta$col = gsub("23-(.*)[wd].*", "\\1", meta$id)

qgraph(relmat, layout='spring', vsize=6,
       cut = 0, repulsion = 1, diag = F,
       minimum = 0.1, maximum = 0.75,
       labels = meta$id,
       #label.cex = 1.2,
       groups = meta$col, palette = "pastel",
       #edge.labels = T,
       legend = F)

View(diag(relmat))
