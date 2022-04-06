#how to install package once
devtools::install_github("stjohn3/R_packages",subdir="BceenetPCAPackage")

#how to load package in the future
library(BceenetPCAPackage)

?BceenetPCAPackage::run.pca.function

##Example
library(librarian)
shelf(seqinr)

setwd("~/Desktop/BCEEnet_ShinyApp/")
fasta.data.glau<-read.alignment(file = "./MakePhylogeny/Glaucomys_ 2192835840/MUSCLE_Fasta",format="fasta")

run.pca.function(data=fasta.data.glau, title.input="PUT YOUR TITLE HERE" )
