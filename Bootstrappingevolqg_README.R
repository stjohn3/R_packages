#how to install package once
devtools::install_github("stjohn3/R_packages",subdir="Bootstrappingevolqg")

#how to load package in the future
library(librarian)
shelf(Boostrappingevolqg)

#help pages 
#this is the main function to use
?Boostrappingevolqg::meanandci()


#View these pages for the specifics
?Boostrappingevolqg::boot.meanmatrixstats()
?Boostrappingevolqg::boot.randomskewer()
?Boostrappingevolqg::boot.pca.similarity()
?Boostrappingevolqg::boot.matrix.distance()