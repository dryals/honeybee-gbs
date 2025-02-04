#select samples to remove

s5 = read.delim("data/s5summary.txt", sep = "")

hist(s5$reads_consens)
hist(s5$nsites)

plot(s5$nsites, s5$reads_consens)

sum(s5$reads_consens < 3e4, na.rm = T)
sum(s5$nsites < 5e6, na.rm = T)
