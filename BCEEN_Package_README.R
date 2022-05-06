#how to install package 
devtools::install_github("stjohn3/R_packages",subdir="BceenetPCAPackage", force=TRUE)

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
      ggrepel,
      fuzzyjoin,
      RColorBrewer,
      scales)

#set the working directory where you want the package to output the .csv
setwd("~/Desktop/BCEEnet_ShinyApp/")

run.pca.function("~/Desktop/BCEEnet_ShinyApp//MakePhylogeny/Glaucomys_ 2192835840/MUSCLE_Fasta", "Glaucomys" ,title.input = "Glaucomys")
run.pca.function("./MakePhylogeny/Nlepida_ 112046140/MUSCLE_file", "lepida" ,title.input = "N. lepida")
run.pca.function("./MakePhylogeny/Mcalifornicus_ 1073931755/MUSCLE_file", "Microtus" ,title.input = "californicus 107")
run.pca.function("./MakePhylogeny/Mcalifornicus_1428083543/MUSCLE_file", "Microtus" ,title.input = "californicus 142")
run.pca.function("./MakePhylogeny/Mcalifornicus_ 148727633/MUSCLE_file", "Microtus" ,title.input = "californicus 148")



