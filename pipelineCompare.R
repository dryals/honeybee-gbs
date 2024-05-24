library(tidyverse)

dr = read.delim("data/bag13_reads_all_.txt", sep = "")
ml = read.delim("data/Keyfile_bee_bag13_ReadsPerSample.log", sep = "")

both.raw = dr %>% select(sample, raw_ipyrad = reads_raw) %>%
  left_join(ml %>% select(sample=FullSampleName, raw_gatk = goodBarcodedReads)) %>%
  pivot_longer(cols = starts_with("raw"))

ggplot(both.raw, aes(x = sample, y = value, fill = name)) + 
  geom_col(position = "dodge") + theme_bw()

ggplot(both.raw, aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(facets = vars(name), nrow = 2) +
  theme_bw()


both.mapped = dr %>% select(sample, mapped_ipyrad = refseq_mapped_reads) %>%
  left_join(ml %>% select(sample=FullSampleName, mapped_gatk = goodReadsMatchedToDataBase)) %>%
  pivot_longer(cols = starts_with("mapped"))

ggplot(both.mapped, aes(x = value)) + 
  geom_histogram() + 
  facet_wrap(facets = vars(name), nrow = 2) +
  labs(y = "samples", x = "number of reads") +
  theme_bw()
  

#remove samples with few clusters
dr %>% arrange(desc(loci_in_assembly)) %>%
  mutate(order = 1:nrow(dr)) %>%
  ggplot(aes(x = order, y = loci_in_assembly)) + geom_col()


dr %>% arrange(loci_in_assembly) %>% slice(1:20) %>% select(sample, loci_in_assembly)

remove = dr %>% filter(loci_in_assembly < 10000) %>% select(sample)
write.table(remove$sample, "data/bag13-lowqual.txt", col.names = F, row.names = F, quote = F)

#remove drones

drones = dr$sample[grepl("d", dr$sample)]
write.table(drones, "data/bag13-drones.txt", col.names = F, row.names = F, quote = F)

#remove both
remove2 = unique(c(drones, remove$sample))
write.table(remove2, "data/bag13-remove.txt", col.names = F, row.names = F, quote = F)

