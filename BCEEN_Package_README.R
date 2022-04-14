#how to install package 
devtools::install_github("stjohn3/R_packages",subdir="BceenetPCAPackage")

#how to load package 
library(BceenetPCAPackage)

#Look at the help page
?BceenetPCAPackage::run.pca.function

##Example
library(librarian)
shelf(base,
      datasets,
      graphics,
      grDevices,
      methods,
      stats,
      seqinr,
      tidyverse,
      dplyr,
      magrittr,
      factoextra,
      stringr,
      utils,
      ggrepel)

setwd("~/Desktop/BCEEnet_ShinyApp/")
fasta.data.glau<-read.alignment(file = "./MakePhylogeny/Glaucomys_ 2192835840/Full_Fasta",format="fasta")

run.pca.function(data=fasta.data.glau, title.input="PUT YOUR TITLE HERE" )
