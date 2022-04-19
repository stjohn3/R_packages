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


run.pca.function("./MakePhylogeny/Glaucomys_ 2192835840/MUSCLE_Fasta", "Glaucomys" ,title.input = "Glaucomys")
run.pca.function("./MakePhylogeny/Nlepida_ 112046140/MUSCLE_file", "lepida" ,title.input = "N. lepida")
run.pca.function("./MakePhylogeny/Mcalifornicus_ 1073931755/MUSCLE_file", "Microtus" ,title.input = "californicus")
run.pca.function("./MakePhylogeny/Mcalifornicus_1428083543/MUSCLE_file", "Microtus" ,title.input = "californicus")
run.pca.function("./MakePhylogeny/Mcalifornicus_ 148727633/MUSCLE_file", "Microtus" ,title.input = "californicus")
