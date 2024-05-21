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
  