library(dplyr)
library(readxl)

#read file
bc = read.delim("barcodes/Keyfile_bee_bag13.txt")

# #grab first instance of each sample name
# bc.unique = bc %>% group_by(FullSampleName) %>% 
#   slice(1) %>% ungroup()
#
# #write an output for each pool
# write.table(bc.unique %>% filter(GBS_Plate == 1) %>% select(FullSampleName, Barcode),
#             file = "barcodes/P1barcodes.txt",
#             row.names = F, col.names = F, quote = F, sep = "\t")
# 
# write.table(bc.unique %>% filter(GBS_Plate == 2) %>% select(FullSampleName, Barcode),
#             file = "barcodes/P2barcodes.txt",
#             row.names = F, col.names = F, quote = F, sep = "\t")

#second version of barcodes
write.table(bc %>% filter(GBS_Plate == 1) %>% select(FullSampleName, Barcode),
            file = "barcodes/P1barcodes-full.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

#second version of barcodes
write.table(bc %>% filter(GBS_Plate == 2) %>% select(FullSampleName, Barcode),
            file = "barcodes/P2barcodes-full.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

#pull names
bc$name = gsub("^23[-](.*)[d|w].*", "\\1", bc$FullSampleName)
unique(bc$name)




### CBH barcodes
  filename = "barcodes/23CBH_dilution plates_barcoded.xlsx"
  
  bc = data.frame(well = NA, barcodes = NA, sample = NA, plate = NA)
  
  sheets = excel_sheets(filename)
  for(s in sheets){
    
    bc.part = rbind(read_excel(filename, sheet = s, range = "B1:D49"),
                    read_excel(filename, sheet = s, range = "H1:J49"))
    bc.part$plate = s
    
    colnames(bc.part) = c("well", "barcodes", "sample", "plate")
    
    bc = rbind(bc, bc.part)
    
  }

  bc = bc %>% filter(!is.na(sample)) %>% 
    rename(colony = sample) 
  bc$sample = NA

  #create unique sample names
  for(c in unique(bc$colony)){
    bc$sample[bc$colony == c] = 
      paste0(bc$colony[bc$colony == c], "_", 1:sum(bc$colony == c))
  }
  
#write out
  for (PLATE in unique(bc$plate)){
    p = gsub("Plate ", "", PLATE)
    filename = paste0("barcodes/23CBH_formatted/23CBH_", p, ".txt")
    write.table(bc %>% filter(plate == PLATE),
                filename,
                row.names = F, quote = F, col.names = F, sep = '\t')
  }




