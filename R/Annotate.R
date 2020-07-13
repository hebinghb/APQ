Annotate <- function(APQ_data) {
  ids <- row.names(APQ_data)
  ids <- unique(ids)
  uids <- ids[-grep(";",ids)]
  mids <- ids[grep(";",ids)]
  l <- floor(length(uids)/100)
  Gene.names <-data.frame()
  for (i in 1:l) {
    si <- (i-1)*100 + 1
    ei <- i*100
    temp <- uniprot_mapping(uids[si:ei])
    Gene.names <- rbind(Gene.names,temp)
  }
  si <- l*100+1
  ei <- length(uids)
  temp <- uniprot_mapping(uids[si:ei])
  Gene.names <- rbind(Gene.names,temp)
  for (i in mids) {
    temp.ids <- unlist(strsplit(i, ";"))
    temp <- uniprot_mapping(temp.ids)
    temp.gene.names <- paste(temp$Gene.names...primary.., collapse=";") 
    temp.gene.description <- paste(temp$Protein.names, collapse=";")
    temp.data <- data.frame(Entry=i, Gene.names...primary..=temp.gene.names,Protein.names=temp.gene.description)
    Gene.names <- rbind(Gene.names,temp.data)
  }
  row.names(Gene.names) <- Gene.names$Entry
  colnames(Gene.names) <- c("Protein","Gene","Description")
  Gene.names <- Gene.names[row.names(APQ_data),]
  new.APQ.data <-cbind(Gene.names,APQ_data)
  return(new.APQ.data)
}
uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab&columns=id,genes(PREFERRED),protein names'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat <- subset(dat, Entry %in% ids)
  dat
}
