#how to install package 
devtools::install_github("stjohn3/R_packages",subdir="BceenetPCAPackage", force=TRUE)

3#how to load package 
library(BceenetPCAPackage)

#Look at the help page
?BceenetPCAPackage::PC.Labels()

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
      stringi,
      utils,
      ggrepel,
      fuzzyjoin,
      RColorBrewer,
      scales,
      sjmisc,
      ggpubr)

#Run 
setwd("/Users/mickey7210/Desktop/BCEEnet_ShinyApp/Check_Matching_IDs/")

Files<-c(
   "./AlignedFastaFiles/Taricha_75857949.fasta",
   "./AlignedFastaFiles/Elgaria.aligned.fasta",
   "./AlignedFastaFiles/Baeolophus_50844830.fasta",
   "./AlignedFastaFiles/Cyanositta_CA_Cicero_amendedheader.fasta",
   "./AlignedFastaFiles/Charina.aligned.fasta",
   "./AlignedFastaFiles/Diadophis.aligned.fasta",
   "./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta",
   "./AlignedFastaFiles/Microtus_californicus_1073931755_aligned.fasta",
   "./AlignedFastaFiles/Microtus_californicus_148727633.fasta",
   "./AlignedFastaFiles/Glaucomys_2192835840_aligned.fasta",
   "./AlignedFastaFiles/Thomomys.aligned.fasta",
   "./AlignedFastaFiles/Contia.aligned.fasta",
   "./AlignedFastaFiles/Batrachoseps_675617998.fasta",
   "./AlignedFastaFiles/Sorex_1764667868.fasta",
   "./AlignedFastaFiles/Ensatina_eschscholtzii_339521622_aligned_original.fasta",
   "./AlignedFastaFiles/Batrachoseps_nigriventris_1867190792.fasta",
   "./AlignedFastaFiles/Aneides_lugubris_33641515.fasta",
   "./AlignedFastaFiles/Batrachoseps_407098425.fasta",
   "./AlignedFastaFiles/Neotoma_112046140.fasta"
)
Species<-c(
            "Taricha",
            "Elgaria",
            "Baeolophus",
            "Cyanositta",
            "Charina",
            "Diadophis",
            "Artemisiospiza",
            "Microtus_1073",
            "Microtus_1487",
            "Glaucomys",
            "Thomomys",
            "Contia",
            "Batrachoseps_675",
            "Sorex",
            "Ensatina",
            "Batrachoseps_186",
            "Aneides",
            "Batrachoseps_407",
           "Neotoma")


for(i in 1:length(Files)){
   fasta.to.pca(Files[i], Species[i])
}

#### Run one file, Example: ####
fasta.to.pca("./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta","Artemisiospiza" )

## Example to check pipeline ##
### Pipe line####

#Step 1
get.potential.voucher.numbers("./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta")%>%
   match.vernet.to.fasta()->matched.ID.list
  
#Step 2 
subset.fasta.file("./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta",matched.ID.list)%>%
   make.pca.data.frame()->PCA.dataframe

#Step 3   
run.pca.analysis(PCA.dataframe,matched.ID.list,subset.fasta.file("./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta",matched.ID.list))->PCA.results

#step 3.5
PC.Labels(PCA.dataframe,matched.ID.list,subset.fasta.file("./AlignedFastaFiles/Artemisiospiza_belli_397327424.fasta",matched.ID.list))->graph.labels

#Step 5   
plot.PCA.Results(PCA.results, "Testing Pipeline Artemisiospiza", PC.lables.input=graph.labels)
