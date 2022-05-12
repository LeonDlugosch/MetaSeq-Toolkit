rm(list = ls(all=T))

ExtractField = function(x, n, sep = "_"){
  return(unlist(lapply(strsplit(x, sep), "[", n)))
}

library(readr, quietly = T)
out_paths = c(Sys.getenv("faaOut"), Sys.getenv("fnaOut"))
in_paths = c(Sys.getenv("faaDir"), Sys.getenv("fnaDir"))
file = dir()[2]
mincov = 3
minNuc = 210
minAA = minNuc/3 

i = 2
for(i in 1:2){
  if(i == 1){filetype = ".faa"}else{filetype = ".fna"}
  wd = getwd()
  #path = paste(wd, in_paths[i], sep = "/")
  minlen = ifelse(filetype == ".fna", minNuc, minAA)
  name = substr(file, start = 1, stop = nchar(file)-4)
  
  print(paste("Filtering ", ifelse(i == 1, "Proteins", "Genes"), " from ",  substr(file, start = 1, stop = nchar(file)-4), "..."), sep = "")
  
  file.in = paste(path, name, filetype, sep = "")
  file.out = paste(wd, "/", substr(out_paths[i], start = 3, stop = nchar(out_paths[i])), "/", name, "_filtered", filetype, sep = "")
  fasta = read_lines(file)
  seq = as.character(fasta[seq(from = 2, to = length(fasta), by = 2)])
  len = nchar(seq)
  id = as.character(fasta[-seq(from = 2, to = length(fasta), by = 2)])
  fasta = data.frame(SeqID = id, Seq = seq)
  cov = as.numeric(ExtractField(fasta$SeqID, 6, "_"))
  fasta = fasta[which(cov >= mincov & len >= minlen),]
  fasta$SeqID = paste(name, "_", ExtractField(fasta$SeqID, 2, "_"), "_", ExtractField(ExtractField(fasta$SeqID, 6, "_"), 1, " "), sep = "")
  file.out = paste(substr(out_paths[i], start = 3, stop = nchar(out_paths[i])), "/", name, "_filtered", filetype, sep = "")
  print(file.out)
  write.table(file = file.out,  fasta,  sep = "\t",  row.names = F,  col.names = F, quote = F)
}
