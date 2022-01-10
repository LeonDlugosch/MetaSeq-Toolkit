##### Gene Catalogue taxonomy Compiler
### This Skript will integrate Kaiju taxonomy data from the Kaiju Ref-Seq Genomes and Progenomes database. 

rm(list = ls(all=T))
##### Loading Libraries and resources
library(stringr)
library(dplyr)
source("C:/Users/icbmadmin/Dropbox/R/Functions.R")

##### Setting Working Dir
setwd("C:/MG_Analysis_v2/Data/Kaiju_Taxonomy/AOM")
##### Loading ProGenome and RefSeq Annotations

IN = list()
IN[[1]] = as.data.frame(read.delim("AOM_ProGenomes_names.txt", header = F, stringsAsFactors = F))
IN[[2]] = as.data.frame(read.delim("AOM_RefSeq_names.txt", header = F, stringsAsFactors = F))
TAX = list()

for (i in 1:2){
  print(paste("Reading dataframe", i, "..."))
  ### Preparing dataframes and Splitting Taxonomy in Taxonomic levels
  IN[[i]] = as.data.frame(IN[[i]][,-1])
  colnames(IN[[i]]) = c("SeqID", "KaijuID", "Taxonomy")
  split = str_split(IN[[i]]$Taxonomy, ";")
  
  if (i == 1){ 
    TAX[[i]] = data.frame(SeqID = IN[[i]]$SeqID,
                          Pro.Kingdom = trimws(unlist(lapply(split, "[[", 1)), which = "left"),
                          Pro.Phylum = trimws(unlist(lapply(split, "[[", 2)), which = "left"),
                          Pro.Class = trimws(unlist(lapply(split, "[[", 3)), which = "left"),
                          Pro.Order = trimws(unlist(lapply(split, "[[", 4)), which = "left"),
                          Pro.Family = trimws(unlist(lapply(split, "[[", 5)), which = "left"),
                          Pro.Genus = trimws(unlist(lapply(split, "[[", 6)), which = "left"),
                          Pro.Species = trimws(unlist(lapply(split, "[[", 7)), which = "left"))
  }
  if (i == 2){ 
    TAX[[i]] = data.frame(SeqID = IN[[i]]$SeqID,
                          Ref.Kingdom = trimws(unlist(lapply(split, "[[", 1)), which = "left"),
                          Ref.Phylum = trimws(unlist(lapply(split, "[[", 2)), which = "left"),
                          Ref.Class = trimws(unlist(lapply(split, "[[", 3)), which = "left"),
                          Ref.Order = trimws(unlist(lapply(split, "[[", 4)), which = "left"),
                          Ref.Family = trimws(unlist(lapply(split, "[[", 5)), which = "left"),
                          Ref.Genus = trimws(unlist(lapply(split, "[[", 6)), which = "left"),
                          Ref.Species = trimws(unlist(lapply(split, "[[", 7)), which = "left"))
  }
  #TAX[[i]] = as.data.frame(TAX[[i]][,-3]) 
  TAX[[i]][TAX[[i]] == "NA"] = NA
  TAX[[i]]$SeqID = ExtractField(as.character(TAX[[i]]$SeqID), 1, ";")
  
  
  if (i == 2){
    print("Merging dataframes...")
    combined.TAX = full_join(TAX[[2]], TAX[[1]], by = "SeqID", .keep_all = T)
    ProG = combined.TAX[,c(1,9:15)]
    RefS = combined.TAX[,c(1:8)]
    SeqID = combined.TAX$SeqID 
    #rm(IN, TAX, split, combined.TAX)
    
    #RefS = head(RefS, n = 500)
    #ProG = head(ProG, n = 500)
    
    na.ProG = apply(ProG, 1, function(x) sum(is.na(x)))
    na.RefS = apply(RefS, 1, function(x) sum(is.na(x)))
    
    RefS = paste(RefS$Ref.Kingdom, RefS$Ref.Phylum, RefS$Ref.Class, RefS$Ref.Order, RefS$Ref.Family, RefS$Ref.Genus, RefS$Ref.Species, sep = ";")
    ProG = paste(ProG$Pro.Kingdom, ProG$Pro.Phylum, ProG$Pro.Class, ProG$Pro.Order, ProG$Pro.Family, ProG$Pro.Genus, ProG$Pro.Species, sep = ";")
    int.Tax = ifelse(na.RefS <= na.ProG, RefS, ProG)
    DB = ifelse(na.RefS <= na.ProG, "RefSeq", "ProGenomes")
    
    split = str_split(int.Tax, ";")
    int.Tax = data.frame(SeqID = SeqID,
                         DB = DB,
                         Kingdom = unlist(lapply(split, "[[", 1)),
                         Phylum = unlist(lapply(split, "[[", 2)),
                         Class = unlist(lapply(split, "[[", 3)),
                         Order = unlist(lapply(split, "[[", 4)),
                         Family = unlist(lapply(split, "[[", 5)),
                         Genus = unlist(lapply(split, "[[", 6)),
                         Species = unlist(lapply(split, "[[", 7)))
    int.Tax[int.Tax == "NA"] = NA
    int.Tax$TaxID = as.numeric(as.factor(as.character(paste(int.Tax$Kingdom, int.Tax$Phylum, int.Tax$Class, int.Tax$Order, int.Tax$Family, int.Tax$Genus, int.Tax$Species, sep = ";"))))
  }
}

TaxIDs = levels(as.factor(int.Tax$TaxID))
for (i in 1:length(TaxIDs)){
  if(i == 1){
    
    for (k in 1:ncol(int.Tax)){
      int.Tax[,k] = as.character(int.Tax[,k])  
    }
    new.tax1 = as.data.frame(NULL)
    new.tax2 = as.data.frame(NULL)
    new.tax3 = as.data.frame(NULL)
    new.tax4 = as.data.frame(NULL)
    #names(new.tax) = names(int.Tax)[-10]
    SAR.search = glob2rx("SAR*")
    c = 0
  }
  
  ### Getting Taxonomy of Kaiju ID 
  #print("0")
  ID.df = int.Tax[which(int.Tax$TaxID == TaxIDs[i]),1:9]
  ID = ID.df[1,]
  
  ### 1. Skipping ID if Kingdom is not knwon (= no information available) OR no taxonomy level is missing 
  #print("1")
  if(is.na(ID$Kingdom[1]) | sum(is.na(ID[1,])) == 0){
    c = c+nrow(ID.df)
    print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); ", round(c/nrow(int.Tax)*100, 3),"% done", sep =""))
    if(c/nrow(int.Tax) <= 0.25){new.tax1 = rbind(new.tax1, ID.df)
    }
    if(0.25 < c/nrow(int.Tax) & c/nrow(int.Tax) <= 0.50){new.tax2 = rbind(new.tax2, ID.df)
    }
    if(0.50 < c/nrow(int.Tax) & c/nrow(int.Tax) <= 0.75){new.tax3 = rbind(new.tax3, ID.df)
    }
    if(0.75 < c/nrow(int.Tax)){new.tax4 = rbind(new.tax4, ID.df)
    }
    next}
  ### 2. Substitution of NAs in unranked Eukaryota
  #print("2")
  if(ID$Kingdom[1] == "Eukaryota"){
    for (f in 3:9){if(is.na(ID[1,f])){f.c = f}}
    for (r in 9:3){if(is.na(ID[1,r])){r.c = r}}
    if(r.c != 9){
      if ((r.c != f.c) | (r.c == f.c & !is.na(ID[1,r.c+1]))){
        ID[1,f.c:r.c] = paste(ID[1,f.c-1],"(unranked)", sep = " ")
      }
    }
  }
  ### 3. Substitution of NAs in unranked Viruses
  #print("3")
  if(ID$Kingdom[1] == "Viruses"){
    ID[,c(4:5)] = "Viruses (unranked)"
  }
  ### 4. Substitution on NA in species if genus is known
  #print("4")
  if(!is.na(ID$Genus[1]) & is.na(ID$Species[1])){
    if(ID$Kingdom[1] != "Viruses"){
      ID$Species = paste(ID$Genus[1], "sp.", sep = " ")
    }else{
      ID$Species = paste("unid.", ID$Genus[1], sep = " ")
    }  
  }
  
  ### 5. Substituting NAs in taxonomy with last known taxonomic rank
  #print("5")
  if(ID$Kingdom[1] != "Viruses"){
    for (j in 3:9){
      if(j == 3){m = NULL}
      if(j == 9 & !is.na(ID[1,j])){break}
      if(is.na(ID[1,j])){
        
        if (is.null(m)){
          ID[,j] = paste("unid.", ID[1,j-1], sep = " ")
          m = ID[1,j-1]
        }else{
          ID[,j] = paste("unid.", m, sep = " ")  
        }
      }else{m = NULL}
    }
  }else{
    for (j in 5:9){
      if(j == 3){m = NULL}
      if(j == 9 & !is.na(ID[1,j])){break}
      if(is.na(ID[1,j])){
        if (is.null(m)){
          ID[1,j] = paste("unid.", ExtractField(ID[1,j-1], 1, sep = " "), sep = " ")
          m = ExtractField(ID[1,j-1], 1, sep = " ")
        }else{
          ID[1,j] = paste("unid.", m, sep = " ")  
        }
      }else{m = NULL}
    } 
  }
  ### 6. Cyanobacteria Class substitution
  #print("6")
  if(ID$Phylum[1] == "Cyanobacteria"){
    ID$Class[1] = "Cyanophyceae"
  }
  
  ### 7. Substitutes NA's in Order-Genus in SAR86, SAR92, SAR202, SAR324 and SAR116 sequneces to "SAR86" or "SAR116" (and so on...) respectively  
  #print("7")
  if(j == 9){
    SAR = ExtractField(as.character(ID[1,9]), 1, " ")
    if (grepl(SAR.search, SAR)){
      ID[,6:8] = SAR
    }
    ID[1,8] = ifelse(grepl("HOT", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl("HTCC", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl("MED", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl("HIMB", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl("SCGC ", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]])-1, sep = " "), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), sep = " ")),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl("KM", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" GC", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" HF", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" EB", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" DG", ID[1,9], ignore.case = F),
                     paste(gsub("unid. ", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" GW", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" RBG", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" SCN", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" SM", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" Ant", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" MedDCM", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" GOM", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" MS0", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    
    ID[1,8] = ifelse(grepl(" QH", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    ID[1,8] = ifelse(grepl(" QS", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
    ID[1,8] = ifelse(grepl(" SW", ID[1,9], ignore.case = F),
                     paste(gsub("unid.", "", ID[1,8]), ExtractField(ID[1,9], n = length(strsplit(ID[1,9], " ")[[1]]), " "), sep =  " "),
                     ID[1,8])
  }
  
  ### 8. Replacing old information with new Taxonomy
  for (k in 3:9){
    ID.df[,k] = ID[1,k]
  }
  
  c = c+nrow(ID.df)
  print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); ", round(c/nrow(int.Tax)*100, 3),"% done", sep =""))
  if(c/nrow(int.Tax) <= 0.25){new.tax1 = rbind(new.tax1, ID.df)
  #print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); Writing in df1", sep =""))
  }
  if(0.25 < c/nrow(int.Tax) & c/nrow(int.Tax) <= 0.50){new.tax2 = rbind(new.tax2, ID.df)
  #print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); Writing in df2", sep =""))
  }
  if(0.50 < c/nrow(int.Tax) & c/nrow(int.Tax) <= 0.75){new.tax3 = rbind(new.tax3, ID.df)
  #print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); Writing in df3", sep =""))
  }
  if(0.75 < c/nrow(int.Tax)){new.tax4 = rbind(new.tax4, ID.df)
  #print(paste(i, " of ", length(TaxIDs), " (KaijuID: ", TaxIDs[i],"; ", nrow(ID.df), " sequences); Writing in df4", sep =""))
  }
  
  if(i == length(TaxIDs)){
    new.tax = rbind(new.tax1, new.tax2, new.tax3, new.tax4)
    rm(new.tax1, new.tax2, new.tax3, new.tax4)}
  
}
write.csv2(file = "AOM_Integratend_ProGenome_RefSeq_taxonomy_v3.csv", new.tax, row.names = F)


bla = new.tax[grepl("SCGC ", new.tax$Species, ignore.case = F),]
