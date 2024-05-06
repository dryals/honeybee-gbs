library(dplyr)

#read file
bc = read.delim("barcodes/Keyfile_bee_bag13.txt")

#grab first instance of each sample name
bc.unique = bc %>% group_by(FullSampleName) %>% 
  slice(1) %>% ungroup()

#write an output for each pool
write.table(bc.unique %>% filter(GBS_Plate == 1) %>% select(FullSampleName, Barcode),
            file = "barcodes/P1barcodes.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(bc.unique %>% filter(GBS_Plate == 2) %>% select(FullSampleName, Barcode),
            file = "barcodes/P2barcodes.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")
