APQ<-function(ms_data, ...) UseMethod("APQ")
APQ.default <- function(ms_data, datatype, ...) {
     	if (!datatype %in% c("DIA","DDA")) {
		   stop("Please select a valid datatype. e.g.:  \"DIA\", \"DDA\"")
	   }
	   if(datatype == "DIA"){
	       total_ms2 <- as.list(tapply(ms_data$PeakArea,ms_data$SampleName,sum))
         samples <- unique(ms_data$SampleName)
         final_APQ <- list()
         for(sample_name in samples){
            sample_data <- subset(ms_data, ms_data$SampleName==sample_name) 
            protein_ms2 <- as.matrix(tapply(sample_data$PeakArea,sample_data$ProteinName,sum))
            new_ms2 <- Redistribution(protein_ms2)
            protein_APQ <- 1000*new_ms2/as.numeric(total_ms2[sample_name])
            final_APQ <- cbind(final_APQ,protein_APQ)
         }
         colnames(final_APQ) <- samples
         final_APQ[final_APQ == 0] <- NA
         if(!all(!is.na(final_APQ))){
            final_APQ <- final_APQ[apply(final_APQ, 1, function(x) !all(is.na(x))),]
         }
         return(final_APQ)
     }
     if(datatype == "DDA"){
       ms_data$Proteins<-as.character(ms_data$Proteins)
       total_ms1 <- apply(ms_data[-1],2,sum)
       samples <- colnames(ms_data[-1])
       final_APQ <- list()
       for(sample_name in samples){
          protein_ms1 <- as.matrix(tapply(ms_data[[sample_name]],ms_data$Proteins,sum))
          new_ms1 <- Redistribution(protein_ms1)
          protein_APQ <- 1000*new_ms1/as.numeric(total_ms1[sample_name])
          final_APQ <- cbind(final_APQ,protein_APQ)
       }
       colname <- gsub("Intensity.","",samples)
       colnames(final_APQ) <- colname
       final_APQ[final_APQ == 0] <- NA
       final_APQ <- final_APQ[apply(final_APQ, 1, function(x) !all(is.na(x))),]   
       final_APQ <- Clean(final_APQ)
       final_APQ[final_APQ == 0] <- NA
       if(!all(!is.na(final_APQ))){
          final_APQ <- final_APQ[apply(final_APQ, 1, function(x) !all(is.na(x))),]       
       }
       return(final_APQ)
     }
}

Redistribution <- function(protein_ms){
      proteins <- rownames(protein_ms)
      share2uniq_ms <- protein_ms
      proteinfrequency <- vector()
      for(protein in proteins[grep(";", proteins)]){ 
          total_uniq_ms <- 0 
          for(temp_protein in unlist(strsplit(protein,";"))){
             if(temp_protein %in% proteins && !is.na(protein_ms[temp_protein,1])){ 
               total_uniq_ms <- total_uniq_ms + protein_ms[temp_protein,1]
             }
             

          }
          
          if(total_uniq_ms > 0 & !is.na(protein_ms[protein,1])){
             for(temp_protein in unlist(strsplit(protein,";"))){
                if(temp_protein %in% proteins && !is.na(protein_ms[temp_protein,1])){ 
                  share2uniq_ms[temp_protein,1] <- share2uniq_ms[temp_protein,1] + (protein_ms[temp_protein,1]/total_uniq_ms)*protein_ms[protein,1]
                  }
             }
             share2uniq_ms[protein,] <- NA
          }
       }
      return(share2uniq_ms)
} 

Clean <- function(APQ){
      proteins <- rownames(APQ)
      clean_apq <- APQ
      clean_apq[is.na(clean_apq)] <- 0
      #Count frequency of protein in protein groups
      proteinfrequency <- vector()
      for(protein in proteins[grep(";", proteins)]){
          for(temp_protein in unlist(strsplit(protein,";"))){
             if(is.na(proteinfrequency[temp_protein])){
                if(temp_protein %in% proteins){
                  proteinfrequency[temp_protein] <- 1000
                }else{
                  proteinfrequency[temp_protein] <- 1
                }
             }else{
                proteinfrequency[temp_protein] <- proteinfrequency[temp_protein]+1
             } 
           }        
      }
      #Clean low confidence protein groups
      for(protein in proteins[grep(";", proteins)]){
          temp_proteins <- unlist(strsplit(protein,";"))
          if(max(proteinfrequency[temp_proteins]) != min(proteinfrequency[temp_proteins])){
             max_proteins <- c()
             for(temp_protein in temp_proteins){
                 if(proteinfrequency[temp_protein] == max(proteinfrequency[temp_proteins])){
                    max_proteins <- c(max_proteins,temp_protein)
                 }
             }
             new_protein <- paste0(sort(max_proteins),collapse=";")
             temp_rowname <- rownames(clean_apq)
             if(!new_protein %in% temp_rowname){
               temp_rowname[which(temp_rowname == protein)] <- new_protein
               rownames(clean_apq) <- temp_rowname
             }else{
             clean_apq[new_protein,] <- unlist(clean_apq[new_protein,]) + unlist(clean_apq[protein,]) 
             clean_apq[protein,] <- NA
             }
          }else{
               new_protein <- paste0(sort(temp_proteins),collapse=";")
               temp_rowname <- rownames(clean_apq)
               if(!new_protein %in% temp_rowname){  
                 temp_rowname[which(temp_rowname == protein)] <- new_protein
                 rownames(clean_apq) <- temp_rowname
               }else if (new_protein!=protein){
                 clean_apq[new_protein,] <- unlist(clean_apq[new_protein,]) + unlist(clean_apq[protein,]) 
                 clean_apq[protein,] <- NA
	       }
          }
      } 
      return(clean_apq)
}
