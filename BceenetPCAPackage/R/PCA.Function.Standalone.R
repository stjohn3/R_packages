#' make file to assign groups
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#'@param file.path Place the file path to your aligned fastafile here... PUT IT IN QUOTES
#'@param Specimen.name What specimen are we working with? choose one of the following: "Aneides lugubris,Artemisiospiza belli belli,Batrachoseps nigriventris,Batrachoseps incognitus,Batrachoseps luciae,Batrachoseps gavilanensis,Batrachoseps attenuatus,Batrachoseps sp.,Charina bottae bottae,Charina bottae,Cyanocitta stelleri carbonacea,Diadophis punctatus vandenburghii,Elgaria multicarinata multicarinata,Elgaria multicarinata,Microtus californicus californicus,Sorex ornatus californicus,Sorex ornatus salarius,Taricha torosa,Thomomys bottae navus,Batrachoseps major,Diadophis punctatus,Microtus californicus sanctidiegi,Neotoma bryanti intermedia,Sorex ornatus ornatus,Thomomys bottae bottae,Microtus californicus kernensis,Microtus californicus aestuarinus,Sorex ornatus,Sorex ornatus sinuosus,Thomomys bottae pascalis,Contia tenuis,Contia longicaudae,Cyanocitta stelleri frontalis,Elgaria coerulea,Ensatina eschscholtzii oregonensis,Microtus californicus eximius,Sorex vagrans vagrans,Sorex ornatus sinuosus x Sorex vagrans vagrans,Taricha rivularis,Thomomys bottae laticeps,Artemisiospiza belli canescens,Elgaria multicarinata webbii,Microtus californicus mohavensis,Neotoma lepida lepida,Thomomys bottae perpallidus,Thomomys bottae riparius,Thomomys bottae albatus,Microtus californicus vallicola,Elgaria panamintina,Baeolophus inornatus inornatus,Glaucomys sabrinus flaviventris,Glaucomys sabrinus fuliginosus,Glaucomys sabrinus lascivus,Taricha granulosa,Baeolophus ridgwayi,Ensatina eschscholtzii platensis,Batrachoseps bramei,Batrachoseps gregarius,Batrachoseps relictus,Batrachoseps simatus,Batrachoseps stebbinsi,Batrachoseps regius,Batrachoseps altasierrae,Batrachoseps robustus,Batrachoseps diabolicus,Batrachoseps kawia,Ensatina eschscholtzii ssp.,Glaucomys sabrinus LASCIVUS,Neotoma bryanti x Neotoma lepida,Neotoma macrotis streatori,Batrachoseps minor,Neotoma fuscipes bullatior,Charina umbratica,Diadophis punctatus modestus,Glaucomys sabrinus californicus,Thomomys bottae nigricans
#'@return Returns a dataframe matching IDs to locations and assigns them colors
#'@export
assign.individuals.to.ecoregions<-function(file.path.fasta, Specimen.name=character()){
  
  #get fasta headers to extract voucher numbers
  read.fasta(file = file.path(file.path.fasta))->temp.fasta
  getAnnot(temp.fasta)%>%unlist()%>%as.data.frame()->annotation.list
  names(annotation.list)<-"annotation"
  
  #print(annotation.list)
  matched.data.frame<<-NULL
  
  #print("start for loop")
  for(i in 1:length(annotation.list$annotation)){
    #print(i)
    #print(annotation.list$annotation[i])
    #Fasta
    scan(text = annotation.list$annotation[i], what = "", quiet=TRUE)%>%as.data.frame()->to.match
    
    to.match%<>%
      str_split(., "MVZ")%>%as.data.frame()
    
    names(to.match)<-"potential.catalognumber"
    
    #print("from scan")
    #print((to.match))
    
    to.match%<>%
      filter(!stringr::str_detect(to.match$potential.catalognumber, ">"),
             !stringr::str_detect(to.match$potential.catalognumber, "-"))%>%
      #stringr::str_detect(to.match$potential.catalognumber, "[0-9]+")==TRUE,
      #!stringr::str_detect(to.match$potential.catalognumber, ">"),
      #stringr::str_detect(potential.catalognumber, "([A-Z]+)[0-9]+")==TRUE)%>%
      mutate(potential.catalognumber=str_extract(potential.catalognumber, "[0-9]+"))%>%
      na.omit()
    
    if(nrow(to.match)==0){
      fill.missing<-data.frame(potential.catalognumber="Missing")
      matched.data.frame<<-rbind(matched.data.frame, fill.missing)
    }
    else{matched.data.frame<<-rbind(matched.data.frame, to.match)}
    
  }
  
  
  working.label.df<-ecoregions%>%
    filter(grepl(Specimen.name, scientific))%>%
    dplyr::select(catalognum,decimallat,decimallon )%>%
    dplyr::mutate(catalog.number.only=str_extract(catalognum, "(?<=:)[0-9]+"))%>%
    dplyr::mutate(catalog.number.only=ifelse(is.na(catalog.number.only)==TRUE, catalognum, catalog.number.only))%>%
    unique()
  
  stringdist_left_join(matched.data.frame,working.label.df,  
                       by=c("potential.catalognumber"="catalog.number.only"),
                       method = c("lv"),
                       max_dist = 0)%>%
    dplyr::select(potential.catalognumber,catalognum,decimallat,decimallon)->out
  return(out)
}

#' Make PCA
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#'@param file.path Place the file path to your aligned fasta file here... PUT IT IN QUOTES
#'@param Specimen.name What specimen are we working with? choose one of the following: "Aneides lugubris,Artemisiospiza belli belli,Batrachoseps nigriventris,Batrachoseps incognitus,Batrachoseps luciae,Batrachoseps gavilanensis,Batrachoseps attenuatus,Batrachoseps sp.,Charina bottae bottae,Charina bottae,Cyanocitta stelleri carbonacea,Diadophis punctatus vandenburghii,Elgaria multicarinata multicarinata,Elgaria multicarinata,Microtus californicus californicus,Sorex ornatus californicus,Sorex ornatus salarius,Taricha torosa,Thomomys bottae navus,Batrachoseps major,Diadophis punctatus,Microtus californicus sanctidiegi,Neotoma bryanti intermedia,Sorex ornatus ornatus,Thomomys bottae bottae,Microtus californicus kernensis,Microtus californicus aestuarinus,Sorex ornatus,Sorex ornatus sinuosus,Thomomys bottae pascalis,Contia tenuis,Contia longicaudae,Cyanocitta stelleri frontalis,Elgaria coerulea,Ensatina eschscholtzii oregonensis,Microtus californicus eximius,Sorex vagrans vagrans,Sorex ornatus sinuosus x Sorex vagrans vagrans,Taricha rivularis,Thomomys bottae laticeps,Artemisiospiza belli canescens,Elgaria multicarinata webbii,Microtus californicus mohavensis,Neotoma lepida lepida,Thomomys bottae perpallidus,Thomomys bottae riparius,Thomomys bottae albatus,Microtus californicus vallicola,Elgaria panamintina,Baeolophus inornatus inornatus,Glaucomys sabrinus flaviventris,Glaucomys sabrinus fuliginosus,Glaucomys sabrinus lascivus,Taricha granulosa,Baeolophus ridgwayi,Ensatina eschscholtzii platensis,Batrachoseps bramei,Batrachoseps gregarius,Batrachoseps relictus,Batrachoseps simatus,Batrachoseps stebbinsi,Batrachoseps regius,Batrachoseps altasierrae,Batrachoseps robustus,Batrachoseps diabolicus,Batrachoseps kawia,Ensatina eschscholtzii ssp.,Glaucomys sabrinus LASCIVUS,Neotoma bryanti x Neotoma lepida,Neotoma macrotis streatori,Batrachoseps minor,Neotoma fuscipes bullatior,Charina umbratica,Diadophis punctatus modestus,Glaucomys sabrinus californicus,Thomomys bottae nigricans
#'@param title.input Place the title of your graph and file here-- PUT IT IN QUOTES
#'@return This function returns a PCA describing variation in the FASTA file input, and a .csv output file that roughly
#'groups individuals into clusters so you can match points on the PCA to an individual. 
#'@export
run.pca.function<-function(file.path.fasta, Specimen.name ,title.input){
  lat.long.dataframe<-assign.individuals.to.ecoregions(file.path.fasta, Specimen.name)
  
  fasta.data<-read.alignment(file = file.path(file.path.fasta),format="fasta")
  ####get number of individuals in dataset####
  num.inds<-as.numeric(fasta.data$nb)
  
  #print(num.inds)
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
    filter(similarity<=.95)->list.of.positions.to.keep.for.PCA
  
  ####PCA####
  #select columns
  store.df%>%
    dplyr::select(list.of.positions.to.keep.for.PCA$position)%>%
    dplyr::select(starts_with("X"))->pca.graph.data
  
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
  
  PCA.ggplot.data<-cbind(PCA.ggplot.data, lat.long.dataframe)
  
  data.frame(PC1.group=round(PCA.ggplot.data[,1], digits=0),
             PC2.group=round(PCA.ggplot.data[,2], digits=0),
             Specimen.ID=PCA.ggplot.data$potential.catalognumber)%>%
    group_by(PC1.group, PC2.group)%>%
    dplyr::summarise(number.inds=n(),
                     list.IDs=str_c(Specimen.ID, collapse = ","))%>%
    write.csv(., file.path("./", paste0(toString(title.input),"_table.csv")))
  
  ### Color Calculations ####
  max.lat<-42
  max.long<-124
  
  
  PCA.ggplot.data%<>%
    mutate(xcol=decimallat/ max.lat,
           ycol=(abs(decimallon)/ max.long))%>%
    mutate(xcol=rescale(xcol, to=c(0,1)),
           ycol=rescale(ycol, to=c(0,1)))
  
  PCA.ggplot.data$xcol<-replace_na(PCA.ggplot.data$xcol, 0)
  PCA.ggplot.data$ycol<-replace_na(PCA.ggplot.data$ycol, 0)
  
  PCA.ggplot.data%<>%
    mutate(color.assignment=ifelse(is.na(decimallon)=="TRUE", "#808080", rgb(1, xcol, ycol)))
  
  
  col <- PCA.ggplot.data$color.assignment%>%as.data.frame()%>%unique()
  names(col)<-"color"
  
  
  Final.PCA.plot<-ggplot(data=PCA.ggplot.data, aes(x=PC1, y=PC2, colour=color.assignment))+
    geom_point(data=PCA.ggplot.data, aes(x=PC1, y=PC2,colour=color.assignment),
               size=4,position = position_jitter(width=.5, height=.5), alpha=.6)+
    scale_colour_manual(values = col$color)
  
  
  
  #Final.PCA.plot<-ggplot(as.data.frame(PCA.ggplot.data), aes(x=PC1, y=PC2, label=row.names(PCA.ggplot.data)))+
  #  geom_point(size=4,position = position_jitter(width=.5, height=.5), alpha=.6)+
  #  geom_label_repel()
  
  Final.PCA.plot+
    ggtitle(toString(title.input))+
    theme_classic(18)+
    theme(legend.position = "none")%>%
    return()
  
}

