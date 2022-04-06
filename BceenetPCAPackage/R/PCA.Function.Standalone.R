#' Make PCA
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param data This parameter accepts a FASTA file containing the samples of interest. Use the seqinr package and the 
#'"read.alignment" function to read data into R. Make sure your FASTA files have been aligned!
#'@param title.input This parameter accepts a character string (make sure to put your string in quotes!). This string will be used
#'as the title of the PCA and as part of the .csv output name
#' @return This function returns a PCA describing variation in the FASTA file input, and a .csv output file that roughly
#' groups individuals into clusters so you can match points on the PCA to an individual. 
#' @export
run.pca.function<-function(data, title.input){
  fasta.data<-data
  
  ####get number of individuals in dataset####
  num.inds<-as.numeric(fasta.data$nb)
  
  ####get number of loci in the data####
  num.loci<-str_length(fasta.data$seq[[2]])
  
  #head(fasta.data)
  
  #####for loop to take the strings of ATCG's for each individual and: #####
  ####1) break them into single characters, ######
  ####2) convert them to numbers######
  ####and 3) calculate proportion of individuals with each nucleotide ######
  
  store.df<<-matrix(ncol = num.loci, nrow = 0)%>%as.data.frame()
  
  for(i in 1:num.inds){
    ind.sample<-str_split(unlist(fasta.data$seq[[i]]), pattern="")%>%unlist()
    s2n(ind.sample)%>%return()
    
    store.df<<-rbind(store.df,s2n(ind.sample))
  }
  
  count.variants<-function(data){
    zero<-0
    one<-0
    two<-0
    three<-0
    missing<-0
    for(i in 1:length(data)){
      
      if(is.na(data[i])==TRUE){missing<-missing+1}
      else if(data[i]==0){zero<-zero+1}
      else if(data[i]==1){one<-one+1}
      else if(data[i]==2){two<-two+1}
      else if(data[i]==3){three<-three+1}
      
    }
    return(paste(zero, one, two, three, missing))
  }
  
  
  lapply(store.df, count.variants)%>%as.matrix()->count.of.SNPS
  rownames(count.of.SNPS)->positions
  
  count.of.SNPS.seperated<-data.frame(position=as.character(positions),
                                      count.A=as.numeric(0),
                                      count.T=as.numeric(0),
                                      count.C=as.numeric(0),
                                      count.G=as.numeric(0))
  
  for(i in 1:nrow(count.of.SNPS)){
    str_split(unlist(count.of.SNPS[i]), pattern=" ")%>%
      unlist()%>%as.data.frame()->counts.seperated
    
    count.of.SNPS.seperated$count.A[i]<-as.numeric(counts.seperated[1,1]) #zero is A, which corresponds to row 1
    count.of.SNPS.seperated$count.T[i]<-as.numeric(counts.seperated[4,1]) #Three is T, which corresonds to row 4
    count.of.SNPS.seperated$count.C[i]<-as.numeric(counts.seperated[2,1]) #one is C, which corresponds to row 2
    count.of.SNPS.seperated$count.G[i]<-as.numeric(counts.seperated[3,1]) #two is G, which corresponds to row 3
  }
  
  
  count.of.SNPS.seperated%<>%
    mutate(individuals=num.inds,
           prop.A=(count.A/num.inds), #zero is A
           prop.T=(count.T/num.inds), #Three is T
           prop.C=(count.C/num.inds), #one is C
           prop.G=(count.G/num.inds)) #two is G
  #return(count.of.SNPS.seperated)
  
  ### filter loci####
  count.of.SNPS.seperated%>%
    dplyr::select("position", "prop.A", "prop.T", "prop.C", "prop.G")%>%
    pivot_longer(cols=prop.A:prop.G,
                 names_to = "Nucleotide",
                 values_to = "Proportion")%>%
    group_by(position)%>%
    dplyr::summarise(similarity=max(Proportion))%>%
    filter(similarity<=.9)->list.of.positions.to.keep.for.PCA
  
  ####PCA####
  #select columns
  store.df%>%
    dplyr::select(list.of.positions.to.keep.for.PCA$position)%>%
    dplyr::select(starts_with("X0"))->pca.graph.data
  
  
  replace.na.custom<-function(data){
    replace_na(data, as.numeric(0.05))
  }
  
  lapply(pca.graph.data, replace.na.custom)->pca.graph.data
  
  #Oh yes this was my problem! When we use princomp Your columns cannot exceed
  #your samples, ignore for now. 
  #princomp(pca.graph.data[,1:10])
  
  locus.pca<-prcomp(as.data.frame(pca.graph.data))
  
  row.names(locus.pca$x)<-fasta.data$nam%>%unlist()

  #######################
  ##OLD plotting of PCA##
  #######################
  # I want to keep this here incase we want to go pack ot this way of visualizing
  #Final.PCA.plot<-fviz_pca_ind(locus.pca, 
  #                             geom=c("point", "text"),
  #                             pointsize = 3, 
  #                             repel = TRUE,
  #                             geom.var=c("point", "text"))
  
  PCA.ggplot.data<-locus.pca$x[,c(1:2)]
  
  data.frame(PC1.group=round(PCA.ggplot.data[,1], digits=0),
             PC2.group=round(PCA.ggplot.data[,2], digits=0),
             Specimen.ID=row.names(PCA.ggplot.data))%>%
    group_by(PC1.group, PC2.group)%>%
    dplyr::summarise(number.inds=n(),
      list.IDs=str_c(Specimen.ID, collapse = ","))%>%
    write.csv(., file.path("./", paste0(toString(title.input),"_table.csv")))
    
    #mutate(PC1.group=round("PC1", digits=0),
    #       PC2.group=round("PC2", digits=0),
    #       specimen.id=row.names(PCA.ggplot.data))
    
  
  Final.PCA.plot<-ggplot(PCA.ggplot.data, aes(x=PC1, y=PC2, label=row.names(PCA.ggplot.data)))+
    geom_point(size=4,position = position_jitter(width=.5, height=.5), alpha=.6)+
    geom_label_repel()
    
  Final.PCA.plot+
    ggtitle(toString(title.input))+
    theme_classic(18)%>%
    return()
  
}
