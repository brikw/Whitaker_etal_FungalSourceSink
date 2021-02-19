# functions to summarize blast results

#summarize blast results -- took too long...
summarize_hits_PCRsim <- function(blasthits){
  
  # categorize into Plant vs Fungi hits
  blasthits %>%
    mutate(SubjectID_type = ifelse(grepl("k__Fungi", blasthits$SubjectID), "Fungi","Plant")) -> blasthits
  
  # pick the top hit for each read
  blasthits %>%
    group_by(QueryID) %>%
    summarize(max.e = max(E)) -> read.maxe
  read.maxe <- data.frame(read.maxe)
  tophit.list <- list()
  i<-0
  for(i in 1:dim(read.maxe)[1]){
    curr.read <- read.maxe[i,"QueryID"]
    curr.maxe <- read.maxe[i,"max.e"]
    blasthits %>%
      filter(QueryID == curr.read) %>%
      filter(E == curr.maxe) -> curr.maxe.df
    tophit.list[[i]] <- unique(curr.maxe.df$SubjectID_type)
    print(paste(i, "of", dim(read.maxe)[1]))
  }
  names(tophit.list) <- read.maxe$QueryID
  tmp <- data.frame(SubjectID_type = unlist(tophit.list))
  tophit.df <- data.frame(QueryID = row.names(tmp), tmp, row.names = NULL)
  
  return(tophit.df)
  
}

# summarize blast hits
summ_hits_quickly <- function(blasthits){
  
  # add type
  blasthits %>%
    mutate(SubjectID_type = ifelse(grepl("k__Fungi", blasthits$SubjectID), "Fungi","Plant")) -> blasthits.cat
  
  # pick the top hit for each read
  blasthits.cat %>%
    group_by(QueryID) %>%
    summarize(max.e = max(E)) -> read.maxe
  read.maxe <- data.frame(read.maxe)
  
  # subset based on max.e for each read and summarize by QueryID (ASV)
  blasthits.cat %>%
    left_join(read.maxe) %>%
    filter(E == max.e)  %>%
    group_by(QueryID) %>%
    summarize(n = length(SubjectID),
              mean.Perc.Ident = mean(Perc.Ident),
              sd.Perc.Ident = sd(Perc.Ident),
              uniq.Subject = paste(unique(SubjectID), collapse = "----"),
              uniq.Subject.Type = paste(unique(SubjectID_type), collapse = "----")) -> topblasthit.asv
  
  return(topblasthit.asv)
  
}

# quantify the proportion of panicum reads in each sample
calc_propPlantReads <- function(tophits, asvlist){
  
  # make a list of Plant ASVs
  #unique(tophits$uniq.Subject.Type) # there are no ambigous Plant/Fungi ASVs
  tophits %>%
    filter(uniq.Subject.Type == "Plant") -> plant.asv.df
  
  # select plant ASVs from the ASV table
  asv.tab <- asvlist$asv.tab
  asv.plant.tab <- asv.tab[,colnames(asv.tab) %in% plant.asv.df$QueryID]
  plant.reads <- rowSums(asv.plant.tab)
  all.reads <- rowSums(asv.tab)
  
  # calc proportion
  df <- data.frame(sample.name.match = row.names(asv.tab), plant.reads, all.reads)
  df %>%
    mutate(prop.plant.reads = plant.reads/all.reads) -> df
  
  return(df)
  
}

# function to remove plant ASVs from the ASV table and add this as a new table to the asvlist object
remove_plantASVs <- function(tophits, asvlist){
  
  # make a list of Plant ASVs
  type.levels <- unique(tophits$uniq.Subject.Type) # there are no ambigous Plant/Fungi ASVs
  both <- c("Plant----Fungi","Fungi----Plant")
  if(sum(both %in% type.levels) !=0){
    print("Warning! There is at least 1 ASV that is ambigiously classified")
    result <- NA
  }else{
    tophits %>%
      filter(uniq.Subject.Type != "Plant") -> noplant.asv.df
    
    # select non-Plant ASVs from the ASV table and add it as a new table to the asvlist object
    asv.tab <- asvlist$asv.tab
    asv.noplant.tab <- asv.tab[,colnames(asv.tab) %in% noplant.asv.df$QueryID]
    result <- list(asv.tab = asvlist$asv.tab, 
                   asv.tab.noplant = asv.noplant.tab,
                   asv.indx = asvlist$asv.indx)
    
  }
  
  return(result)
  
}
