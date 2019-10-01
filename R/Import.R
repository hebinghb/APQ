Import<-function(filename, filetype, ...) UseMethod("Import")

Import.default <- function(filename, filetype, ...) {
	if (!filetype %in% c("spectronaut","skyline")) {
		stop("Please select a valid filetype. Options:  \"spectronaut\", \"skyline\"")
	}
	if (!grepl("\\.csv$",filename) && !grepl("\\.tsv$",filename) ) {
		stop("Please use a valid file format. e.g.:  \"csv\", \"tsv\"")
	}
	if(filetype == "spectronaut"){ 
	  cat("Filetype: spectronaut\n")
    if (grepl("\\.csv$",filename)) {
    cat("File format: csv\n")
		# Spectronaut csv
		cat("Filename: ",filename,"\n")
		cat("Loading data...\n")
		raw_data <- read.csv(file=filename)
		cat("Done\n")
		cat("Transforming data...\n")	
		ms_data <- raw_data[c("R.FileName","PG.ProteinAccessions","EG.StrippedSequence","F.PeakArea")]
		ms_data <- subset(ms_data,ms_data$F.PeakArea != "NaN")
		colnames(ms_data) <- c("SampleName","ProteinName","PeptideSequence","PeakArea")
		cat("Done\n")
		return(ms_data)
   }
    if (grepl("\\.tsv$",filename)) {
    cat("File format: tsv\n")
		# Spectronaut tsv
		cat("Filename: ",filename,"\n")
		cat("Loading data...")
		raw_data <- read.csv(file=filename,sep = "\t")
		cat("Done\n")
		cat("Transforming data...\n")			
		ms_data <- raw_data[c("R.FileName","PG.ProteinAccessions","EG.StrippedSequence","F.PeakArea")]
		ms_data <- subset(ms_data,ms_data$F.PeakArea != "NaN")
		colnames(ms_data) <- c("SampleName","ProteinName","PeptideSequence","PeakArea")			 
		cat("Done\n")
		return(ms_data)
   } 
  } else if(filetype == "skyline") {
    cat("Filetype: skyline\n")
    if (grepl("\\.csv$",filename)) {
    cat("File format: csv\n")
		# Skyline csv
		cat("Filename: ",filename,"\n")
		cat("Loading data...\n")
		raw_data <- read.csv(file=filename)
		cat("Done\n")
		cat("Transforming data...\n")			
		ms_data <- raw_data[c("Replicate.Name","Protein.Name","Peptide.Sequence","Fragment.Ion","Area")]
		ms_data <- subset(ms_data,ms_data$Area != "#N/A")
		ms_data <- Reformat(ms_data)
		ms_data <- ms_data[c("Replicate.Name","Protein.Name","Peptide.Sequence","Area")]
		colnames(ms_data) <- c("SampleName","ProteinName","PeptideSequence","PeakArea")
		ms_data$PeakArea <- as.numeric(as.character(ms_data$PeakArea))
		cat("Done\n")
		return(ms_data)
    }
    if (grepl("\\.tsv$",filename)) {
    cat("File format: tsv\n")
		# Skyline tsv
		cat("Filename: ",filename,"\n")
		cat("Loading data...\n")		
		raw_data <- read.csv(file=filename,sep = "\t")
		cat("Done\n")
		cat("Transforming data...\n")			
		ms_data <- raw_data[c("Replicate.Name","Protein.Name","Peptide.Sequence","Fragment.Ion","Area")]
		ms_data <- subset(ms_data,ms_data$Area != "#N/A")
		ms_data <- Reformat(ms_data)
		ms_data <- ms_data[c("Replicate.Name","Protein.Name","Peptide.Sequence","Area")]
		colnames(ms_data) <- c("SampleName","ProteinName","PeptideSequence","PeakArea")
		ms_data$PeakArea <- as.numeric(as.character(ms_data$PeakArea))
		cat("Done\n")
		return(ms_data)
    }      
  }
}

Import.DDA <- function(filename, ...) {
	if (!grepl("peptides\\.txt$",filename) ) {
		stop("Please use peptides.txt from MaxQuant")
	}
		# MaxQuant peptides
		cat("Filetype: maxquant\n")
		cat("Filename: ",filename,"\n")
		cat("Loading data...\n")
		raw_data <- read.csv(file=filename,sep = "\t")
		cat("Done\n")
		cat("Transforming data...\n")	
		colname <- colnames(raw_data)
		# Remove Contamination
		Proteins <- raw_data$Proteins
		if(length(grep("CON_",Proteins))!=0){
		   temp_data <- raw_data[-(grep("CON_",Proteins)),]
		}else{
		   temp_data <- raw_data
		}
		Proteins <- temp_data$Leading.razor.protein
		if(length(grep("REV_",Proteins))!=0){
		   clean_data <- temp_data[-(grep("REV_",Proteins)),]
		}else{
		   clean_data <- temp_data
		}
		ms_data <- clean_data[c("Proteins",colname[grep("Intensity.",colname)])]
		temp_data <- clean_data[c(colname[grep("Intensity.",colname)])]
		ms_data <- ms_data[apply(temp_data, 1, function(x) !all(x == 0)),]
		cat("Done\n")
		return(ms_data)
  
}

Reformat <- function(ms_data){
   temp_data1 <- data.frame(ID=paste0(ms_data$Replicate.Name,ms_data$Peptide.Sequence,ms_data$Area,ms_data$Fragment.Ion),Protein.Name=ms_data$Protein.Name)
   temp_data2 <- unique(data.frame(ID=paste0(ms_data$Replicate.Name,ms_data$Peptide.Sequence,ms_data$Area,ms_data$Fragment.Ion),Replicate.Name=ms_data$Replicate.Name,Peptide.Sequence=ms_data$Peptide.Sequence,Area=ms_data$Area))
   temp_data1 <- as.data.frame.table(by(data=temp_data1,INDICES = temp_data1$ID,FUN = function(x){paste0(x$Protein.Name,collapse = ';')}))
   colnames(temp_data1) <- c("ID","Protein.Name")
   temp_data <- merge(temp_data1,temp_data2)
   temp_data <- temp_data[c("Replicate.Name","Protein.Name","Peptide.Sequence","Area")]
   return(temp_data)
   
}
