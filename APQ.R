APQ<-function(ms_data, ...) UseMethod("APQ")
APQ.default <- function(ms_data, datatype, ...) {
     	if (!datatype %in% c("DIA","DDA")) {
		   stop("Please select a valid datatype. e.g.:  \"DIA\", \"DDA\"")
	   }
	   if(datatype == "DIA"){
     total_ms2 <- as.list(tapply(ms_data$F.PeakArea,ms_data$R.FileName,sum))
     samples <- unique(ms_data$R.FileName)
     final_APQ <- list();
     for(sample_name in samples){
        sample_data <- subset(ms_data, ms_data$R.FileName==sample_name) 
        protein_ms2 <- as.matrix(tapply(sample_data$F.PeakArea,sample_data$PG.ProteinAccessions,sum))
        new_ms2 <- Redistribution(protein_ms2)
        protein_APQ <- 1000*new_ms2/as.numeric(total_ms2[sample_name])
        final_APQ <- cbind(final_APQ,protein_APQ)
     }
     colnames(final_APQ) <- samples
     return(final_APQ)
     }
     if(datatype == "DDA"){
     total_ms1 <- apply(ms_data[-(1:3)],2,sum)
     samples <- colnames(ms_data[-(1:3)])
     final_APQ <- list();
     for(sample_name in samples){
        protein_ms1 <- na.omit(as.matrix(tapply(ms_data[[sample_name]],ms_data$Proteins,sum)))
        new_ms1 <- Redistribution(protein_ms1)
        protein_APQ <- 1000*new_ms1/as.numeric(total_ms1[sample_name])
        final_APQ <- cbind(final_APQ,protein_APQ)
     }
     colname <- gsub("Intensity.","",samples)
     colnames(final_APQ) <- colname
     
     return(final_APQ)
     }
}

Redistribution <- function(protein_ms){
      proteins <- rownames(protein_ms)
      share2uniq_ms <- protein_ms
      for(protein in proteins[grep(";", proteins)]){ 
          total_uniq_ms <- 0 
          for(temp_protein in unlist(strsplit(protein,";"))){
             if(temp_protein %in% proteins){ 
               total_uniq_ms <- total_uniq_ms + protein_ms[temp_protein,1]
             }
          }
          
          if(total_uniq_ms > 0){
             for(temp_protein in unlist(strsplit(protein,";"))){
                if(temp_protein %in% proteins){ 
                  share2uniq_ms[temp_protein,1] <- share2uniq_ms[temp_protein,1] + (protein_ms[temp_protein,1]/total_uniq_ms)*protein_ms[protein,1]
                  }
             }
             share2uniq_ms[protein,1] <- 0
          }
       }
      return(share2uniq_ms)
}
