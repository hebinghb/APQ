# Label-Free Absolute Protein Quantification with Data-Independent Acquisition

This is an R package for label-free absolute protein quantification (APQ) using data-independent acquisition (DIA).

This package was developed based on TPA method and an algorithm that redistribute the MS signals from shared peptides to individual proteins or isoforms. If there is any question or error, please email me: hebinghb@gmail.com

### Citation 
Bing He, Jian Shi, Xinwen Wang, Hui Jiang, Hao-Jie Zhu."Label-free absolute protein quantification with data-independent acquisition."*J Proteomics*. 2019 May 30;200:51-59. doi: [10.1016/j.jprot.2019.03.005](https://doi.org/10.1016/j.jprot.2019.03.005).

##
### Installation of "APQ" package on R:
#### Run this command if package "devtools" wasn't installed
    install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
#### Load package "devtools"
    library(devtools)
#### Install package "APQ"
    devtools::install_github("hebinghb/APQ")

##
### Demo code for DIA-APQ (DIA-TPA) analysis:
#### Load package "APQ"
    library(APQ)
#### Load DIA data
##### For MS2 data from Spectronaut
    data<-Import(filename="20181028_202636_HLS9_36_WILD_Report.csv",filetype="spectronaut")

    #Note: 20181028_202636_HLS9_36_WILD_Report.csv should include following columns: "R.FileName", "PG.ProteinAccessions", "EG.StrippedSequence", "F.PeakArea". CSV and TSV format are supported in current version.
##### For MS2 data from Skyline
    data<-Import(filename="Transition Results.csv",filetype="skyline")

    #Note: "Transition Results.csv"  should include following columns: "Replicate.Name", "Protein.Name", "Peptide.Sequence", "Fragment.Ion", "Area". CSV and TSV format are supported in current version.
#### APQ analysis
    quantity<-APQ(data,"DIA") #DIA indicates you are using DIA data.
#### Write APQ result to a file
    write.csv(quantity,file="DIA_demo.csv")
#### Protein id annotation
    Annotated.quantity<-Annotate(quantity) #Annotate protein id with gene symbol and description. Note: This function needs internet connection.
#### Write Annotated APQ result to a file
    write.csv(Annotated.quantity,file="DIA_demo.csv",row.names=FALSE)



###
##
>**Note:** This package also supports APQ analysis using data-dependent acquisition (DDA).
### Demo code for DDA-APQ (DDA-TPA) analysis:
#### Load package "APQ"
    library(APQ)
#### Load DDA data
    data<-Import.DDA("peptides.txt")
    #Note: peptides.txt is the output file of MaxQuant. It locates in /combined/txt/ under your MaxQuant output directory. 
#### APQ analysis
    quantity<-APQ(data,"DDA") #DDA indicates you are using DDA data.
#### Write APQ result to a file
    write.csv(quantity,file="DDA_demo.csv")
