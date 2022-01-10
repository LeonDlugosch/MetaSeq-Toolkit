### 1 .Download KO Brite from KEGG (https://www.kegg.jp/brite/ko00001) 
### 2. Open with text editor (Notepad++), remove first 3 and last 5 lines (the stuff that looks like html and does not belong to KO reference)
### 3. Save as .txt

### Libraries and custom functions
library(stringr)

ExtractField = function(x, n, sep = "_"){
  return(unlist(lapply(strsplit(x, sep), "[", n)))
}
DeleteField = function(x, n, sep = "_"){
  list = strsplit(x, sep)
  l = unlist(rapply(list, f = length, how = "list"))
  for(i in 1:length(list)){
    if(i == 1){new.strings = NULL}
    seq = 1:l[i]
    seq = seq[-n]
    for(j in 1:length(seq)){
      if(j == 1){
        string = list[[i]][seq[j]]
      }else{
        string = paste(string, list[[i]][seq[j]], sep = sep)
      }
    }
    new.strings = c(new.strings, string)
  }
  return(new.strings)
}
PasteVector = function(v = NULL, sep = "_"){
  if(is.null(v)){
    stop("No input!")
  }
  for(i in 1:length(v)){
    if(i == 1){
      r = v[1]
    }else{
      r = paste(r, v[i], sep = sep)
    }
  }
  return(r)
}

### converting to tabular format
KOs = as.data.frame(readLines("KO_BRITE.txt"))
names(KOs) = "V1"
KOs$V1 = str_squish(KOs$V1)

for(i in 1:nrow(KOs)){
  if(i == 1){
    KO.out = as.data.frame(matrix(nrow = 0, ncol = 6))
    names(KO.out) = c("KNr", "A", "B", "C", "Gene", "Function")
    c = 1
    }
  class = ExtractField(KOs$V1[i], 1," ")
  if(substr(class, 1,1) == "A"){A = DeleteField(KOs$V1[i], 1, " ")}
  if(substr(class, 1,1) == "B"){B = DeleteField(KOs$V1[i], 1:2, " ")}
  if(substr(class, 1,1) == "C"){C = DeleteField(KOs$V1[i], 1:2, " ")}
  if(substr(class, 1,1) == "D"){
    KO.out[c,1] = ExtractField(KOs$V1[i], 2," ")
    KO.out[c,2] = A
    KO.out[c,3] = B
    KO.out[c,4] = C
    
    gene = ExtractField(KOs$V1[i], 3:length(strsplit(KOs$V1[i], " ")[[1]])," ")
    
    fc = 0
    for (k in 1:length(gene)){
      if(substr(gene[k], nchar(gene[k]), nchar(gene[k])) == ","){
          fc = fc+1
      }else{
        break
        } 
    }
    KO.out[c,5] = gsub(";", "", PasteVector(ExtractField(KOs$V1[i], 3:(3+fc), " "), sep = " "))
    
    KO.out[c,6] = DeleteField(KOs$V1[i], 1:(3+fc), " ")
    c = c+1
  }
}

### dereplicating KOs and adjusting for presence in multiple pathways
KNrs = levels(as.factor(KO.out$KNr))
for(i in 1:length(KNrs)){
  if(i == 1){
    KOs.adj = as.data.frame(matrix(nrow = 0, ncol = 6))
    names(KOs.adj) = c("KNr", "A", "B", "C", "Gene", "Function")
    c = 1
  }
  tmp = KO.out[which(KO.out$KNr == KNrs[i]),]
  tmp.vec = as.data.frame(tmp[1,])
  
  if(length(table(tmp$C)) == 1){
    KOs.adj[c,] = tmp.vec
    c = c+1
    next
  }
  if(length(table(tmp$B)) == 1 & length(table(tmp$C) > 1)){
    tmp.vec[1,4] = paste("general", tmp.vec[1,3], sep = " ")
    KOs.adj[c,] = tmp.vec
    c = c+1
    next
  }
  if(length(table(tmp$A)) == 1 & length(table(tmp$B) > 1)){
    tmp.vec[1,c(3,4)] = paste("general", tmp.vec[1,2], sep = " ")
    KOs.adj[c,] = tmp.vec
    c = c+1
    next
  }
  if(length(table(tmp$A)) > 1){
    tmp.vec[1,c(2:4)] = "multiple functions" 
    KOs.adj[c,] = tmp.vec
    c = c+1
    next
  }
  if(i == length(KNrs)){
    KOs.adj$A = tolower(KOs.adj$A)
    KOs.adj$B = tolower(KOs.adj$B)
    KOs.adj$C = tolower(KOs.adj$C)
  }
}

write.table(x = KOs.adj, file = "KO_refTable.txt", sep = "\t", quote = F, row.names = F)
