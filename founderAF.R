#create founder allele frequencies using pedigree and alphaPeel

#load libraries
  setwd("/home/dylan/Documents/bees/harpurlab/project/honeybee-gbs")
  #setwd("~/ryals/honeybee-gbs")
  library(tidyverse)
  #library(Matrix)
  library(readxl)
  library(vcfR)



#create KING: see cbh_analysissh

### Load and sort data
  
  #load KING
    king = read.delim("data/23CBH.king", header = F)
    ids = read.delim("data/23CBH.king.id")
      ids$id = paste0(ids$X.FID, "_", ids$IID)
      ids$pedid = gsub("CBH", "_", ids$X.FID)
      
    colnames(king) = rownames(king) = ids$id
    king = as.matrix(king)
    king[is.infinite(king)] = NA
    
    #remove outliers!
    RM = rowMeans(king, na.rm = T)
    hist(RM)
    RMnames = names(RM)[RM < -10]
    dim(king)
    king = king[!colnames(king) %in% RMnames, !colnames(king) %in% RMnames]
    dim(king)
    
  
  #load pedigree
    ped = read_excel("../gensel/data2025/pedigree_edit_15DEC25.xlsx")
    ped.f = ped %>% filter(queen_id %in% ids$pedid)
    
  #load pedigree relatedness matrix
    load("/home/dylan/Documents/bees/harpurlab/project/gensel/data2025/qwA_15DEC25.Rdat")
    
    # setwd("/home/dylan/Documents/bees/harpurlab/misc/ped/AINV-honeybees/v20/")
    # load("CBH_Adirect_28NOV25.Rdata")
    # 
    # ident = read.table("pedigree_ident.txt")
    # names = ident %>% select(1, 6) %>% rename(ident = 1, dam = 2)
    # names$id = names$ident
    # names$id[grepl("^30", names$id)] = 
    #   paste0(names$dam[grepl("^30", names$id)], "_w")
    # 
    # colnames(A) = rownames(A) = names$id
    # 
    # setwd("/home/dylan/Documents/bees/harpurlab/project/honeybee-gbs")
    
  #filter ped matrix
    workers = paste0(colnames(A)[colnames(A) %in% ids$pedid],"_w")
    A.f = A[workers, workers]
    
### detect pedigree outliers from KING
    
    #shortking
    u = unique(ids %>% filter(!id %in% RMnames) %>% select(X.FID) %>% unique )
    u = u$X.FID
    shortking = matrix(NA, nrow = length(u), ncol = length(u))
    colnames(shortking) = rownames(shortking) = u
    
    for(i in 1:length(u)){
      for(j in 1:i){
        
        #rel within col
        if(i == j){
          colking = king[grepl(u[i], colnames(king)), grepl(u[i], colnames(king))]
          shortking[i,j] = mean(colking[(upper.tri(colking, diag = F))], na.rm = T)
        } else{
          
          #rel between cols
          allgenrel = king[grepl(u[i], colnames(king)), grepl(u[j], colnames(king))]
          shortking[i,j] = mean(unlist(as.vector(allgenrel)), na.rm = T)
          
        }

      }
    }
    
    
#compare!
    
    trans = data.frame(genid = u)
    trans$pedid = gsub("CBH", "_", trans$genid)
    trans$pedid = paste0(trans$pedid, "_w")
    
    trans2 = trans %>% filter(genid %in% colnames(shortking),
                              pedid %in% colnames(A.f)) %>% 
      mutate(queen_id = gsub("_w", "", pedid)) %>% 
      left_join(ped %>% select(queen_id, mother))
    
   shortking.full = forceSymmetric(shortking, uplo = "L") %>% as.matrix

   A.out = A.f[trans2$ped, trans2$pedid]
   #A.out = A.out[upper.tri(A.out, diag = F)]
   k.out = shortking.full[trans2$gen, trans2$genid]
   k2.out = k.out
    k2.out[k2.out<0] = 0
   #k.out = k.out[upper.tri(k.out, diag = F)]
   # 
   # mat.comp = (A.out - (2 * k2.out)) / A.out
   # mat.comp = (k2.out - (0.5 * A.out))
   # 
   # hist(colMeans(abs(mat.comp)))
   # 
   # goodcols = trans2[colMeans(abs(mat.comp)) < 0.0075,]
   # 
   # table(goodcols$mother)
   # 
   
   
   #  
   #  
   #  badcols = trans2[colMeans(mat.comp) > 0.015,]
   #  badcols = badcols %>% 
   #    left_join(ped.f %>% 
   #                mutate(pedid = paste0(queen_id, "_w")) %>% 
   #                select(pedid, mother))
   #  table(badcols$mother)
    
    
    
#sort out errors
   
   
  brds = trans2 %>% group_by(mother) %>% summarise(n = n()) %>% 
    filter(n>2)
  
  suspect = vector("character")
  for(B in brds$mother){
    
    #B="II87"
    
    #A.test = A.out[trans2$mother == B, trans2$mother == B]
    k.test = k.out[trans2$mother == B, trans2$mother == B]
    # k2.test = k.test
    #   k2.test[k2.test<0] = 0
    
    #comp.test = abs(A.test - (k2.test * 2))
    #RM = base::rowMeans(comp.test, na.rm = T) %>% scale()
    RM = rowMeans(k.test, na.rm = T) %>%  scale()
    RM = RM[,1]
    
    hist(RM, main = B)
    
    suspect = append(suspect, names(RM)[RM < -2])
    
  }
  
  trans2 %>% filter(genid %in% suspect)
  
  #better placement?
  brdplace = matrix(nrow = length(suspect), ncol = nrow(brds))
    rownames(brdplace) = suspect
    colnames(brdplace) = brds$mother
    
  for(B in brds$mother){
    
    #B = "BQ02"
    k.test = k.out[suspect, trans2$mother == B & !trans2$genid %in% suspect]
    RS = rowMeans(k.test, na.rm = T)
    
    brdplace[,B] = RS

  }
    
    brdplace
    hist(brdplace %>% as.vector() %>% unlist)
    brdplace[brdplace > -0.08]
  
    
  #remove suspect
  A.good = A.out[!trans2$genid %in% suspect, !trans2$genid %in% suspect]
  k.good = k.out[!trans2$genid %in% suspect, !trans2$genid %in% suspect]
  out = data.frame(ped = A.good[upper.tri(A.good, diag = F)],
                   king = k.good[upper.tri(k.good, diag = F)])

  # #all
  # out = data.frame(ped = A.out[upper.tri(A.out, diag = F)],
  #                  king = k.out[upper.tri(k.out, diag = F)])

  # #focus breeder
  # B = "I69"
  # A.good = A.out[trans2$mother == B, trans2$mother == B]
  # k.good = k2.out[trans2$mother == B, trans2$mother == B]
  # out = data.frame(ped = A.good[upper.tri(A.good, diag = F)],
  #                  king = k.good[upper.tri(A.good, diag = F)])

   ggplot(out, aes(x = ped, y = king)) +
     geom_jitter(width = 0.001, alpha = 0.5) +
     geom_abline(aes(intercept = 0, slope = 1/2)) +
     geom_smooth(method = 'lm', se = F) + 
     lims(x = c(0,0.075))
   
   mean(out$king[out$ped == 0])
   mean(out$king[out$ped == 0.0625])
   mean(out$king[out$ped == 0.125])
   

#write out!
   
   goodworkers = ids %>% 
     filter(!id %in% RMnames,
            !X.FID %in% suspect)
   
   write.table(goodworkers %>% select(X.FID, id), 
               file = "data/goodworkers.txt",
               col.names = F, row.names = F)
   
#prepare pedigree
   #ensure ids sare same
   ped.ap = ped %>% select(-tech) %>% filter(yob < 2024)
   #0 for unknown
   ped.ap$mother[grepl("stock", ped.ap$mother)] = 0
   ped.ap$mating[grepl("stock", ped.ap$mating)] = 0
   ped.ap$mating[ped.ap$mating == "-"] = 0
   ped.ap$mother[ped.ap$mother == "-"] = 0
   #change to genomic IDs
   ped.ap$queen_id = gsub("23_", "23CBH", ped.ap$queen_id)
   ped.ap$mother = gsub("23_", "23CBH", ped.ap$mother)
   #add genids for the breeders
   genbrd = ids$X.FID %>% 
     unique() %>% 
     gsub("23CBH", "",.)
     genbrd = genbrd[grepl("[A-Za-z]", genbrd)]
     genbrd = genbrd[!grepl("[a]", genbrd)]
     
   ped.ap$queen_id[ped.ap$queen_id %in% genbrd] = 
     paste0("23CBH", ped.ap$queen_id[ped.ap$queen_id %in% genbrd])
   
   ped.ap$mother[ped.ap$mother %in% genbrd] = 
     paste0("23CBH", ped.ap$mother[ped.ap$mother %in% genbrd])
   
   ped.ap$mating[ped.ap$mating %in% genbrd] = 
     paste0("23CBH", ped.ap$mating[ped.ap$mating %in% genbrd])
   
   #fix
   ped.ap$mother[ped.ap$mother == "23CBH0"] = 0
   ped.ap$mating[ped.ap$mating == "23CBH0"] = 0
   
   #check
   sum(ped.ap$queen_id %in% ids$X.FID)
   length(unique(ids$X.FID))
   #add all genotyped workers
   new.ap = ids %>% 
     select(queen_id = id, mother = X.FID) %>% 
     mutate(yob = 2023, mating = 0)
   ped.ap = rbind(ped.ap, new.ap)
   #check
   nrow(ped.ap)
   length(unique(ped.ap$queen_id))
   ped.ap %>% group_by(queen_id) %>% 
     summarise(n = n()) %>% filter(n >1)
   
   #translate to physical pedigree!!!
   ped.phys = ped.ap %>% select(queen_id, mother, mating)
   ped.phys$father = 0
   for(i in 1:nrow(ped.phys)){
      if(sum(ped.phys$queen_id == ped.phys$mother[i]) != 0){
       ped.phys$father[i] =
         ped.phys$mating[which(ped.phys$queen_id == ped.phys$mother[i])]
      }
   }
   
  #write out
   write.table(ped.phys %>% select(queen_id, father, mother),
               file = "data/APped.txt", row.names = F, col.names = F)
   


   
   
#write genotype file
   #read in allele frequencies ...
   setwd("/scratch/negishi/dryals/gbs/23CBH/analysis")
  
   #read in genotypes
   vcf.raw = read.vcfR("23CBH-ap.vcf")
   
   
   #convert to dosage
   gt = extract.gt(vcf.raw)
   
   gt2 = gt
   gt2[gt == "0/0"] = "0"
   gt2[gt2 == "1/1"] = "2"
   gt2[gt2 == "1/0" | gt2 == "0/1"] = "1"
   gt2 = as.data.frame(gt2)
   
   rm(gt, vcf.raw)
   
  #write it
   gt.out = cbind(colnames(gt2), t(gt2))
   gt.out[is.na(gt.out)] = 9
   gt.out[1:10,1:10]
   
   write.table(gt.out,
               file = "ap/ap.geno",
               col.names = F,
               row.names = F,
               quote = F)
   
   
