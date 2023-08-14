#' Input fasta file, get annotation strip out potential IDs
#'
#' This function takes a fasta file and outputs a dataframe containing potential ID's to be matched with the vertnet data
#'
#'@param file.path.fasta full path to fasta files.
#'@return Returns a dataframe of potential ID's to match with vernet data
#'@export
get.potential.voucher.numbers<-function(file.path.fasta){
  #read in Fasta file
  read.fasta(file = file.path(file.path.fasta))->temp.fasta

  #Get the headers and make into dataframe
  getAnnot(temp.fasta) %>%
    unlist() %>%
    as.data.frame() -> annotation.list

  #Change column name
  names(annotation.list) <- "annotation"

  #Create new column
  annotation.list$potential.catalognumber<-NULL

  #for loop to grab potential voucher number
  for(i in 1:nrow(annotation.list)){

    # the annotations are sometimes seperated by spaces and somtiems by under scores. this section and the following if else statment strips the voucher
    # numbers depending on the seperator
    scan(text = annotation.list$annotation[i], what = "", sep = c(""), quiet=TRUE)%>%
      as.data.frame()%>%nrow()->seperator.type

    if(seperator.type>1){
      scan(text = annotation.list$annotation[i], what = "", sep = c(""), quiet=TRUE)%>%
        as.data.frame()%>%
        dplyr::rename("potential.catalognumber"=".")%>%
        filter(!stringr::str_detect(potential.catalognumber, ">"),
               !stringr::str_detect(potential.catalognumber, "-"),
               !stringr::str_detect(potential.catalognumber, "\\("),
               stringr::str_detect(potential.catalognumber, "\\d"))%>%
        dplyr::slice(1)%>%
        as.character()->potential.catalognumber

      annotation.list$potential.catalognumber[i]<-potential.catalognumber
    }
    else{scan(text = annotation.list$annotation[i], what = "", sep = c("_"), quiet=TRUE)%>%
        as.data.frame()%>%
        dplyr::rename("potential.catalognumber"=".")%>%
        filter(!stringr::str_detect(potential.catalognumber, ">"),
               !stringr::str_detect(potential.catalognumber, "-"),
               !stringr::str_detect(potential.catalognumber, "\\("),
               stringr::str_detect(potential.catalognumber, "MVZ"))%>%
        dplyr::slice(1)%>%
        as.character()->potential.catalognumber

      annotation.list$potential.catalognumber[i]<-potential.catalognumber
    }

  }

  #make new column to remove any prefix from the voucher number and eliminate any straggler number possibilities
  annotation.list%<>%
    dplyr::mutate(catalog.number.only=str_extract(potential.catalognumber, "[0-9]+"))%>%
    filter(str_length(catalog.number.only)>3)

  return(annotation.list)

}

#' Input dataframe of potential IDs and matches them with vertnet data
#'
#'  This function uses the ecoregions dataset to match up IDs between the vertnet data and the fasta file and to  make columns for genus, species, and subspecies for each sample
#' It then outputs a dataframe containing this info.
#'
#'
#'@param input.file input should be the output from the get.potential.voucher.numbers() function or a dataframe that contains a list of potential IDs to be matched with vertnet data
#'@return returns a list of matched IDs that are found in both the fasta files and the vertnet dataframes
#'@export
match.vernet.to.fasta<-function(input.file){

  input.file->annotation.list

  #grab only relevant columns from ecoregion dataframe
  lat.long<-ecoregions%>%
    dplyr::select(scientific, catalognum, decimallat,decimallon,New_label)%>%
    dplyr::mutate(catalog.number.only=str_extract(catalognum, "(?<=:)[0-9]+"))%>%
    dplyr::mutate(catalog.number.only=ifelse(is.na(catalog.number.only)==TRUE, catalognum, catalog.number.only))%>%
    unique()


  #make comparison lists to use in for loop
  does.this.list<-lat.long$catalognum
  contain.these.values<-annotation.list$catalog.number.only

  #initialize data frame for loop
  matching.id.data.frame<<-data.frame(scientific=character(),
                                      catalognum=character(),
                                      decimallat=numeric(),
                                      decimallon=numeric(),
                                      New_label=character(),
                                      catalog.number.only=numeric(),
                                      annotation=character(),
                                      potential.catalognumber=character(),
                                      catalog.number.only=numeric())


  #for loop asking whether the lat long list from vertnet contains the stripped catalog/voucher numbers from the fasta files.
  for (i in 1:length(does.this.list)) {
    for (a in 1:length(contain.these.values)) {
      str_contains(does.this.list[i], contain.these.values[a]) -> link

      if (link == TRUE) {
        cbind(lat.long[i, ], annotation.list[a, ]) -> int
        matching.id.data.frame <<- rbind(matching.id.data.frame, int)
      }
    }
  }

  #rename columns and reorder columns
  names(matching.id.data.frame)<-c("scientific","catalognum","decimallat","decimallon","ecoregion_label","vernet.catalog.number.only",
                                   "annotation","potential.catalognumber","genbank.catalog.number.only")

  #Make unique columns for genus species and subspecies

  matching.id.data.frame$Genus<-NA
  matching.id.data.frame$Species<-NA
  matching.id.data.frame$subspecies<-NA
  matching.id.data.frame$matching<-NA

  #for loop to grab genus species and subspecies, and to double check that the matched IDs are with the correct genus.
  for (i in 1:nrow(matching.id.data.frame)) {
    ### Make columns for genus species and subspecies
    scan(text = matching.id.data.frame$scientific[i], what = "", quiet = TRUE)[1] -> genus
    scan(text = matching.id.data.frame$scientific[i], what = "", quiet = TRUE)[2] -> species
    scan(text = matching.id.data.frame$scientific[i], what = "", quiet = TRUE)[3] -> subspecies

    ## print(genus)
    ## print(species)
    ## print(subspecies)

    matching.id.data.frame$Genus[i] <- genus
    matching.id.data.frame$Species[i] <- species
    matching.id.data.frame$subspecies[i] <- subspecies


    if(str_contains(matching.id.data.frame$annotation[i],matching.id.data.frame$Genus[i])==TRUE){
      matching.id.data.frame$matching[i] <- 1
    }
    else if(matching.id.data.frame$Genus[i]=="Artemisiospiza" & str_contains(matching.id.data.frame$annotation[i],"Amphispiza")){
      matching.id.data.frame$matching[i] <- 1
    }
    else if(matching.id.data.frame$Genus[i]=="Cyanocitta" & str_contains(matching.id.data.frame$annotation[i],"Cyanosita")){
      matching.id.data.frame$matching[i] <- 1
    }
    else{matching.id.data.frame$matching[i] <- 0
    }

  }

  matching.id.data.frame%>%
    dplyr::filter(matching==1)%>%
    dplyr::select(scientific:ecoregion_label, annotation)->Matched.vertnet.fasta.samples

  return(Matched.vertnet.fasta.samples)
}

#' Subsets fasta file to only include individuals with matching vertnet information
#'
#' This function takes the output from the match.vertnet.to.fasta() function OR a dataframe containing matched Fasta ID's and vernet ID's and subsets the fasta file
#' to only include individuals with matching vertnet info.
#'
#'@param file.path.fasta full path to fasta file.
#'@param list.matching.annotations list of matching fasta and vernet IDs to keep
#'@return returns a list containing sequence information for only individuals with vertnet data
#'@export
fasta.file.subset<-function(file.path.fasta, list.matching.annotations){
  read.fasta(file = file.path(file.path.fasta))->temp.fasta
  list.matching.annotations->keep.list
  my_fasta_sub <- temp.fasta[str_contains(keep.list$annotation,names(temp.fasta))==TRUE]
  return(my_fasta_sub)
}

#' Takes the subsetted sequence data, count variants per site, filter out sites with >95% similarity, make data frame where each column represents a position
#'
#' This function takes the output from fasta.file.subset() or a fasta file (converted to a list using the sequinr package) and produces a dataframe that can be used for a pca.
#'
#'@param subsetted.fasta subsetted list of sequence data
#'@return returns a dataframe where each column is a sequence position, and each row represents an individual. This dataframe is suitable for running a pca analysis
#'@export
make.pca.data.frame<-function(subsetted.fasta){

  fasta.data<-as.alignment(nb = length(subsetted.fasta), nam = names(subsetted.fasta),
                           seq = getSequence(subsetted.fasta), com = NA)

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
    rownames(store.df)[i]<<-fasta.data$nam[[i]]
  }

  colnames(store.df)<- paste0("position", seq(1,900))

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

  print(tail(count.of.SNPS.seperated))

  count.of.SNPS.seperated%<>%
    mutate(individuals=num.inds,
           prop.A=(count.A/num.inds), #zero is A
           prop.T=(count.T/num.inds), #Three is T
           prop.C=(count.C/num.inds), #one is C
           prop.G=(count.G/num.inds)) #two is G

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
    dplyr::select(list.of.positions.to.keep.for.PCA$position)->pca.graph.data
  #dplyr::select(starts_with("X"))->pca.graph.data

  sapply(pca.graph.data, replace_na, value=as.numeric(.01))%>%as.data.frame()->pca.graph.data
  rownames(pca.graph.data)<-fasta.data$nam

  #print(pca.graph.data)
  return(pca.graph.data)
}

#' Runs pca analysis on sequence dataframe
#'
#' This function takes the dataframe from the make.pca.data.frame() function and runs a pca analysis
#'
#'@param fasta.to.pca.data input dataframe from make.pca.data.frame() function
#'@param matched.vertnet.and.fasta input dataframe from match.vernet.to.fasta() function
#'@param subsetted.fasta input dataframe from fasta.file.subset() function
#'@return returns data frame with PC1 and PC2 information matched up per individual from fasta file
#'@export
run.pca.analysis<-function(fasta.to.pca.data, matched.vertnet.and.fasta,subsetted.fasta){

  #Step 1
  fasta.to.pca.data->pca.graph.data
  subsetted.fasta->subset.names

  #Step 2
  locus.pca<-prcomp(as.data.frame(pca.graph.data))

  #Step 3
  row.names(locus.pca$x)<-names(subset.names)

  #Step 4
  PCA.ggplot.data<-locus.pca$x[,c(1:2)]

  PCA.ggplot.data%<>%
    as.data.frame()%>%
    rownames_to_column("FASTA.ID")

  #Step 5
  matched.vertnet.and.fasta%>%
    mutate(FASTA.ID=str_extract(matched.vertnet.and.fasta$annotation, pattern = "(?<=>)[^/ /]+"))->matched.vernet.and.fasta.forjoining

  PCA.ggplot.data<-left_join(PCA.ggplot.data,matched.vernet.and.fasta.forjoining, by=c("FASTA.ID"))%>%
    dplyr::select(annotation, scientific:ecoregion_label, PC1, PC2)

  return(PCA.ggplot.data)

}

#' Grabs PC1 and PC2 percent variance explained for the graph
#'
#' This function takes the dataframe from the make.pca.data.frame() function and runs a pca analysis and outputs a matrix with the labels for the graph
#'
#'@param fasta.to.pca.data input dataframe from make.pca.data.frame() function
#'@param matched.vertnet.and.fasta input dataframe from match.vernet.to.fasta() function
#'@param subsetted.fasta input dataframe from fasta.file.subset() function
#'@return returns data frame with PC1 and PC2 information matched up per individual from fasta file
#'@export
PC.Labels<-function(fasta.to.pca.data, matched.vertnet.and.fasta,subsetted.fasta){
  fasta.to.pca.data->pca.graph.data
  subsetted.fasta->subset.names

  locus.pca<-prcomp(as.data.frame(pca.graph.data))
  summary(locus.pca)["importance"]%>%as.data.frame()->summary.out

  summary.out["Proportion of Variance", 1:2]%>%as.matrix()->Graph.labels

  return(Graph.labels)

}

#' make a table to export that describes the lat long range and PC1&2 values for each ecoregion
#'
#' This function takes the output from run.pca.analysis() function and summarises the information
#'
#'@param data.frame.graph.pca output from run.pca.analysis() function
#'@param title.input informative name for table (e.g. cyprinodon_variegatus will be output as cyprinodon_variegatus_table.csv)
#'@return outputs a .csv file containing summary info. NOTE: Table will be saved to the location of the working directory.
#'@export
make.table.of.ecoregiongroups<-function(data.frame.graph.pca, title.input){

  data.frame.graph.pca->PCA.ggplot.data

  PCA.ggplot.data%>%
    group_by(ecoregion_label)%>%
    dplyr::summarise(number.inds=n(),
                     latitude=paste(min(decimallat), "-", max(decimallat)),
                     longitude=paste(min(decimallon), "-", max(decimallon)),
                     PC1=mean(PC1),
                     PC2=mean(PC2),
                     list.IDs=str_c(catalognum, collapse = ","))%>%
    write.csv(., file.path("./", paste0(toString(title.input),"_table.csv")))

}

#' Plot PCA results
#'
#' This function takes the pca data frame from run.pca.analysis() function, and produces a PCA plot
#'
#'@param data.frame.graph.pca input dataframe from the output of run.pca.analysis() function
#'@param title.input input informative title that will be used in the plot and for saving the figure.
#'@return Saves out a pdf of the pca to the working directory
#'@export
plot.PCA.Results <- function(data.frame.graph.pca, title.input, PC.lables.input) {

  PCA.ggplot.data <- left_join(data.frame.graph.pca, unique(colors),
                               by = c("decimallat", "decimallon"))%>%arrange(annotation)
  PC.lables.use<-PC.lables.input

  pc1.label<-paste0("PC1 ",(PC.lables.use[1]*100)%>%round(., digits=2), "%")
  pc2.label<-paste0("PC2 ",(PC.lables.use[2]*100)%>%round(., digits=2), "%")

  Final.PCA.plot <- ggplot(data = PCA.ggplot.data, aes(x = PC1, y = PC2, colour = color)) +
    geom_point(
      data = PCA.ggplot.data,
      aes(
        x = PC1, y = PC2,
        colour = color,
        shape = ecoregion_label
      ),
      size = 4,
      stroke=2,
      #position = position_jitter(width = .5, height = .5),
      position = position_jitterdodge(jitter.width = .9, jitter.height = .9, dodge.width=.9),
      alpha = 1
    ) +
    scale_colour_identity() +
    scale_shape_manual(values = c("Northern California Coast Ranges and Coast"=16,
                                  "Central California Coast Ranges and Coast"=17,
                                  "Klamath Mountains"=18,
                                  "Southern California Coast"=19,
                                  "Central Valley"=3,
                                  "Southern California Mountains and Valleys"=4,
                                  "Sierra Nevada"=5,
                                  "Southern Cascades"=6,
                                  "Mojave_Sonoran"=8,
                                  "Modoc Plateau"=1,
                                  "Basin"=11,
                                  "Colorado Desert"=0
    ))+
    labs(x=pc1.label, y=pc2.label)+
    theme_classic(18) +
    theme(legend.position = "none")


  Final.PCA.plot +
    ggtitle(toString(title.input))->save.graph
  ggexport(save.graph, filename=toString(paste(getwd(), paste0(title.input,"_PCAgraph.pdf"), sep="/")))

  Final.PCA.plot +
    ggtitle(toString(title.input)) %>%
    return()
}

#' Runs the entire fasta to pca pipeline
#'
#' This function turns a fasta input into a pca based on the matches found in teh vertnet dataframe
#'
#'@param file.path.fasta full path to fasta files.
#'@param save.filename Input informative title to be used for tables and graphs
#'@return returns .csv table with summary information and a PCA graph. Both files will be saved to wherever the working directory is set.
#'@export
fasta.to.pca<-function(file.path.fasta, save.filename){
  #Variables that are reused
  file.location<-file.path.fasta
  title.table<-save.filename
  title.graph<-save.filename
  #step1
  #print("step 1: get.potential.voucher.numbers")
  get.potential.voucher.numbers(file.path.fasta =file.location)->potential.matches
  #step2
  #print("step 2: match.vernet.to.fasta")
  match.vernet.to.fasta(input.file = potential.matches)->final.matched.vertnet.and.fasta
  #step3
  #print("step 3: fasta.file.subset")
  fasta.file.subset(file.path.fasta = file.location,
                    list.matching.annotations = final.matched.vertnet.and.fasta)->subsetted.fasta.file
  #step 4
  #print("step 4: make.pca.data.frame")
  make.pca.data.frame(subsetted.fasta = subsetted.fasta.file)->subsetted.fasta.to.pca.df
  #step 5
  #print("step 5: run.pca.analysis")
  run.pca.analysis(fasta.to.pca.data = subsetted.fasta.to.pca.df,
                   matched.vertnet.and.fasta = final.matched.vertnet.and.fasta,
                   subsetted.fasta = subsetted.fasta.file)->data.frame.graph.pca

  #Step 5b
  #get labels for percent variance graph
  PC.Labels(fasta.to.pca.data = subsetted.fasta.to.pca.df,
                   matched.vertnet.and.fasta = final.matched.vertnet.and.fasta,
                   subsetted.fasta = subsetted.fasta.file)->labels.out
  #step 6A
  make.table.of.ecoregiongroups(data.frame.graph.pca, title.input = as.character(title.table))
  #step 6B
  #print("step 6: Plot pca")
  plot.PCA.Results(data.frame.graph.pca,
                   title.input = as.character(title.graph),
                   PC.lables.input=labels.out)%>%return()
}

