# APQ

#Label-Free Absolute Protein Quantification with Data-Independent Acquisition

#This is an R package for label-free absolute protein quantification (APQ) using data-independent acquisition (DIA).

#This package is developed based on TPA method and an algorithm that redistribute the MS signals from shared peptides to individual proteins or isoforms.


########################################################

#Installation of this package on R:

#Run this command if package "devtools" wasn't installed

install.packages(c("devtools", "roxygen2", "testthat", "knitr"))

#Load package "devtools"

library(devtools)

#Install package "APQ"

devtools::install_github("hebinghb/APQ")


#######################################################

#Demo code for DIA-APQ analysis:

#Load package "APQ"
library(APQ)

#Load DIA data

##For MS2 data from Spectronaut

data<-Import(filename="20181028_202636_HLS9_36_WILD_Report.csv",filetype="spectronaut")

##20181028_202636_HLS9_36_WILD_Report.csv should include following columns: "R.FileName", "PG.ProteinAccessions", "EG.StrippedSequence", "F.PeakArea". CSV and TSV format are supported in current version.

##For MS2 data from Skyline

data<-Import(filename="Transition Results.csv",filetype="skyline")

##"Transition Results.csv"  should include following columns: "Replicate.Name", "Protein.Name", "Peptide.Sequence", "Fragment.Ion", "Area". CSV and TSV format are supported in current version.

#APQ analysis

quantity<-APQ(data,"DIA") #DIA indicates you are using DIA data.

#Write APQ result to a file

write.table(quantity,file="DIA_demo.txt")


#######################################################

#This package also supports APQ analysis using data-dependent acquisition (DDA).

#Demo code for DDA-APQ analysis:

#Load package "APQ"

library(APQ)

#Load DDA data

data<-Import.DDA("peptides.txt")

##peptides.txt is the output file of MaxQuant. It locates in /combined/txt/ under your MaxQuant output directory. 

#APQ analysis

quantity<-APQ(data,"DDA") #DDA indicates you are using DDA data.

#Write APQ result to a file

write.table(quantity,file="DDA_demo.txt")
