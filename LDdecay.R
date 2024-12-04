library(tidyverse)

ld = read.delim("/scratch/bell/dryals/gbs/analysis/indy2.ld", header = T, sep = "")
ld$dist = (ld$BP_B - ld$BP_A) * (ld$CHR_A == ld$CHR_B)

ld.small = ld %>% filter(dist < 300000)


ld.bin = ld.small %>% mutate(bin = cut_interval(dist, length = 100)) %>%
  group_by(bin) %>% summarise(r2 = mean(R2)) %>%
  mutate(dist_kb = seq(0, max(ld.small$dist), 100) / 1000)

ggplot(ld.bin, aes(x = dist_kb, y = r2)) +
  geom_point(alpha = 0.5) + 
  geom_smooth(method = 'loess', span = 0.1) +
  theme_bw()


ld.bin$bin[1]
