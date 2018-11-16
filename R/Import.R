Import<-function(filename, filetype, ...) UseMethod("Import")
Import.default <- function(filename, ...) {
	if (!grepl("\\.csv$",filename) && !grepl("\\.tsv$",filename) ) {
		stop("Please use a valid output filetype from Spectronaut. e.g.:  \"csv\", \"tsv\"")
	}
  if (grepl("\\.csv$",filename)) {
		# Spectronaut csv
		raw_data <- read.csv(file=filename)
		ms_data<-raw_data[c("R.FileName","PG.ProteinAccessions","EG.StrippedSequence","F.PeakArea")]
		return(ms_data)
  }
  if (grepl("\\.tsv$",filename)) {
		# Spectronaut tsv
		raw_data <- read.csv(file=filename,sep = "\t")
		ms_data<-raw_data[c("R.FileName","PG.ProteinAccessions","EG.StrippedSequence","F.PeakArea")]
		return(ms_data)
  } 
  if(grepl("peptides\\.txt$",filename)){
  
  } 
}
Import.DDA <- function(filename, ...) {
	if (!grepl("peptides\\.txt$",filename) ) {
		stop("Please use peptides.txt from MaxQuant")
	}
		# MaxQuant peptides
		raw_data <- read.csv(file=filename,sep = "\t")
		colname <- colnames(raw_data)
		# Remove Contamination
		Proteins <- raw_data$Proteins
		temp_data <- raw_data[-(grep("CON_",Proteins)),]
		Proteins <- temp_data$Leading
		clean_data <- temp_data[-(grep("REV_",Proteins)),]
		ms_data<-clean_data[c("Proteins","Gene.names","Protein.names",colname[grep("Intensity.",colname)])]
		return(ms_data)
  
}