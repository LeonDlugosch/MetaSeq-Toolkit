paths = c(Sys.getenv("faaOut"), Sys.getenv("fnaOut"))
file = Sys.getenv("file")
mincov = Sys.getenv("minCov")
minAA = Sys.getenv("minNuc")
minNuc = minAA*3 


for(i in 1:2){
  setwd(paths[i])
  filetype = substr(file, start = nchar(file)-3, stop = nchar(file))
  minlen = ifelse(filetype = ".fna", minNuc, minAA)
  name = substr(file, start = 1, stop = nchar(file)-4)
  fasta = read.table(paste(name, filetype, sep = ""), header = F)
  seq = as.character(fasta[seq(2, nrow(fasta), 2),])
  len = nchar(seq)
  fasta = as.data.frame(fasta[-seq(2, nrow(fasta), 2),])
  fasta$Seq = seq
  rm(seq)
  names(fasta) = c("SeqID", "Seq")
  cov = as.numeric(unlist(lapply(strsplit(fasta$SeqID, "_"), "[", 6)))
  fasta = fasta[which(cov >= mincov & minlen >= minlen),]
  fasta$SeqID = paste(name, "_", as.numeric(unlist(lapply(strsplit(fasta$SeqID, "_"), "[", 2))), "_", as.numeric(unlist(lapply(strsplit(fasta$SeqID, "_"), "[", 7))), sep = "") 
  write.table(file = paste(name, "_gf", filetype, sep = ""),  fasta,  sep = "\t",  row.names = F,  col.names = F, quote = F)
}