library(tidyverse)


depths = read.delim("data/test-filter.depth", header = F, sep = "")
ids = read.delim("data/test-filter.names", header = F, sep = "")
colnames(depths) = c("chr", "pos", ids$V1)

#snps per chr
depths %>% filter(!grepl("N", chr)) %>% 
  mutate(chr = factor(chr, levels = as.character(1:16))) %>%
  group_by(chr) %>% summarise(n = n()) %>%
  ggplot(aes(x = chr, y = n)) + geom_col() +
  theme_bw()

#typical sample depth by chromosome
depths.long = depths %>% pivot_longer(cols = starts_with("23"))

depths.chrs = depths.long %>% filter(!grepl("N", chr)) %>%
  mutate(chr = factor(chr, levels = as.character(1:16))) 
  
  ggplot(depths.chrs, aes(x = value)) + geom_histogram() +
    facet_wrap(facet = vars(chr), scales = "free")+
    theme_bw()
  
nsamp = length(ids$V1)

#depth for the best and worst samples
sample.depth = depths.chrs %>% group_by(name) %>%
  summarise(mean_depth = mean(value))

sum(sample.depth$mean_depth > 14) / nrow(sample.depth)

mean(sample.depth$mean_depth)

sample.depth %>%
  arrange(desc(mean_depth)) %>% slice(c(1:3, (nsamp -2) : nsamp)) %>%
  ggplot(aes(x = name, y = mean_depth)) + geom_col() +
  theme_bw()
