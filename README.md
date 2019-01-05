# APQ

#Label-Free Absolute Protein Quantification with Data-Independent Acquisition

#This is an R package for label-free absolute protein quantification (APQ) using data-independent acquisition (DIA).

#This package is developed based on TPA method and an algorithm that redistribute the MS signals from shared peptides to individual proteins or isoforms.


########################################################

#Installation on R:

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

library(devtools)

devtools::install_github("hebinghb/APQ")


#######################################################

#Demo code for DIA-APQ analysis:

library(APQ)

data<-Import("20181028_202636_HLS9_36_WILD_Report.csv")

#20181028_202636_HLS9_36_WILD_Report.csv is the output file of Spectronaut. The output file should include following columns:"R.FileName","PG.ProteinAccessions","EG.StrippedSequence","F.PeakArea". CSV and TSV format are supported in current version.

quantity<-APQ(data,"DIA") #DIA indicates you are using DIA data.

write.table(quantity,file="DIA_demo.txt")


#######################################################

#This package also supports APQ analysis using data-dependent acquisition (DDA).

#Demo code for DDA-APQ analysis:

library(APQ)

data<-Import.DDA("peptides.txt") 

#peptides.txt is the output file of MaxQuant. It locates in /combined/txt/ under your MaxQuant output directory. 

quantity<-APQ(data,"DDA") #DDA indicates you are using DDA data.

write.table(quantity,file="DDA_demo.txt")
