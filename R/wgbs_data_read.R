ReadCooMet <- function(Config){

  # Reads the methylation call files and gets the motif targets' coordinates
  # and their methylation frequencie.
  
  # Methylation calls come from pb3seq output
  if (tolower(Config$InputFormat)=="pb3seq"){
	  CooMet <- get_methylation_coordinates_pb3Seq(Config, Config$DirDat)
    return(CooMet)
  }
  # Methylation calls come from bedmethyl outputs
  if (tolower(Config$InputFormat)=="bedmethyl"){
    CooMet <- get_methylation_coordinates_bedMethyl(Config)
    return(CooMet)
  }
  stop("Signal format not supported")
}

get_methylation_coordinates_bedMethyl <- function(config){
  
  files <- list.files(config$DirDat)
  
  if (length(files) == 0) stop("WGBS file not in directory")
  
  if (is.na(config$SignalFile)) {
    # Get methylation rates in the directory
    if (length(files) > 1){
      
      CooMet <- list()
      
      # As there are more than one file, we perform an outer join between each of them and calculate the average
      # methylation rate for each common methylation coordinate.
      CooMetFor <- list()
      CooMetRev <- list()
      CooMets <- list()
      for (i in 1:length(files)) {
        CooMets[[i]] <- get_methylation_coordinates_bedMethyl_(config, config$DirDat, files[i])
      }
      
      # Calculate methylation average for each coordinate.
      for (ci in 1:config$NumChrl){
        
        dfs = list()
        i = 1
        for (cm in CooMets) {
          dfs[[i]] <- cm[[1]][[ci]]
          i = i + 1
        }
        CooMetFor[[ci]] <- as.data.frame(rbindlist(dfs)[,.(ColMet=mean(ColMet)),ColCoo])
        
        dfs <- list()
        i = 1
        for (cm in CooMets) {
          dfs[[i]] <- cm[[2]][[ci]]
          i = i + 1
        }
        CooMetRev[[ci]] <- as.data.frame(rbindlist(dfs)[,.(ColMet=mean(ColMet)),ColCoo])
        
      }
      CooMet[[1]] <- CooMetFor
      CooMet[[2]] <- CooMetRev
      
    } else {
      fl <- files[1]
      
      CooMet <- get_methylation_coordinates_bedMethyl_(config, config$DirDat, fl)
    }
  }
  else {
    CooMet <- get_methylation_coordinates_bedMethyl_(config, config$DirDat, config$SignalFile)
  }
  
  
  return(CooMet)
  
}

get_methylation_coordinates_bedMethyl_ <- function(config, directory, file_name){
  
  # Reads the methylation call files and gets the motif targets' coordinates
  # and their methylation frequencies.
  CooMet <- list()
  CooMetFor <- list()
  CooMetRev <- list()
  
  FulNamFilTxt_for <- file.path(directory, file_name)
  print(paste("Reading file", FulNamFilTxt_for))
  
  # Can read .bed.gz or csv.gz files directly.
  df <- fread(FulNamFilTxt_for, colClasses=c('character', 'numeric', 'NULL','NULL','NULL', 'character','NULL','NULL','NULL','NULL','numeric'))
  
  div100 <- function(x) x/100
  df <- data.frame(df[[1]], df[[2]], df[[3]], df[[4]]/100 )
  colnames(df) <- c("ChrName","ColCoo","Sense","ColMet")
  df_for <- subset(df, Sense=="+")
  
  for(i1 in 1:config$NumAutosomes){
    # print(paste("Autosome",i1,"for"))
    CooMetFor[[i1]] <- subset(df_for, ChrName==paste0("chr",i1)) %>% select(ColCoo, ColMet)
  }
  
  for(i2 in config$Allosomes){
    # print(paste("Allosome",i2,"for"))
    i1 <- i1 + 1
    CooMetFor[[i1]] <- subset(df_for, ChrName==paste0("chr",i2)) %>% select(ColCoo, ColMet) 
  }
  
  rm(df_for)
  
  df_rev <- subset(df, Sense=="-")
  for(i1 in 1:config$NumAutosomes){
    # print(paste("Autosome",i1,"rev"))
    CooMetRev[[i1]] <- subset(df_rev, ChrName==paste0("chr",i1)) %>% select(ColCoo, ColMet)
  }
  
  for(i2 in config$Allosomes){
    i1 <- i1 + 1
    # print(paste("Allosome",i2,"rev"))
    CooMetRev[[i1]] <- subset(df_rev, ChrName==paste0("chr",i2)) %>% select(ColCoo, ColMet) 
  }
  
  CooMet <- list()
  CooMet[[1]] <- CooMetFor
  CooMet[[2]] <- CooMetRev
  return(CooMet)
}

get_methylation_coordinates_pb3Seq <- function(config, directory){
  
  CooMet <- list()
  CooMetFor <- list()
  CooMetRev <- list()
  
  for(i1 in 1:config$NumAutosomes){
    
    # Select current autosome's methylation call file.
    NamFilTxt_for <- paste("chr",i1,".fa.txt_forw.txt_",config$MotifTarget,sep="")
    # print(NamFilTxt_for)
    FulNamFilTxt_for <- file.path(directory,NamFilTxt_for)
    # Get the targets' coordinates' and methylation frequencies' columns from the file.
    DatCooMetFor <- read.table(FulNamFilTxt_for, colClasses=c('NULL', 'numeric', 'NULL','NULL','NULL', 'NULL','numeric','NULL'))
    colnames(DatCooMetFor) <- c("ColCoo","ColMet")
    CooMetFor[[i1]] <- DatCooMetFor
    # Select current autosome's methylation call file.
    NamFilTxt_rev <- paste("chr",i1,".fa.txt_rev.txt_",config$MotifTarget,sep="")
    # print(NamFilTxt_rev)
    FulNamFilTxt_rev <- file.path(config$DirDat,NamFilTxt_rev)
    # Get the targets' coordinates' and methylation frequencies' columns from the file.
    DatCooMetRev <- read.table(FulNamFilTxt_rev, colClasses=c('NULL', 'numeric', 'NULL','NULL','NULL', 'NULL','numeric','NULL'))
    colnames(DatCooMetRev) <- c("ColCoo","ColMet")
    # sub1 <- function(x) x-1
    # Substract one to coordinates for posterior DicWord.
    # DatCooMetRev <- data.frame( lapply(df[1], sub1), DatCooMetRev[2:2] )
    #print(paste("DEV Adding to",i1))
    CooMetRev[[i1]] <- DatCooMetRev
  }
  
  for(i2 in config$Allosomes){
    i1 <- i1 + 1
    # Select current allosome's methylation call file.
    NamFilTxt_for <- paste("chr",i2,".fa.txt_forw.txt_",config$MotifTarget,sep="")
    # print(NamFilTxt_for)
    FulNamFilTxt_for <- file.path(config$DirDat,NamFilTxt_for)
    # Get the targets' coordinates' and methylation frequencies' columns from the file.
    DatCooMetFor <- read.table(FulNamFilTxt_for, colClasses=c('NULL', 'numeric', 'NULL','NULL','NULL', 'NULL','numeric','NULL'))
    colnames(DatCooMetFor) <- c("ColCoo","ColMet")
    CooMetFor[[i1]] <- DatCooMetFor  
    # Select current allosome's methylation call file.
    NamFilTxt_rev <- paste("chr",i2,".fa.txt_rev.txt_",config$MotifTarget,sep="")
    # print(NamFilTxt_rev)
    FulNamFilTxt_rev <- file.path(config$DirDat,NamFilTxt_for)
    # Get the targets' coordinates' and methylation frequencies' columns from the file.
    DatCooMetRev <- read.table(FulNamFilTxt_rev, colClasses=c('NULL', 'numeric', 'NULL','NULL','NULL', 'NULL','numeric','NULL'))
    colnames(DatCooMetRev) <- c("ColCoo","ColMet")
    #print(paste("DEV Adding to",i1))
    CooMetRev[[i1]] <- DatCooMetRev
  }
  
  CooMet <- list()
  CooMet[[1]] <- CooMetFor
  CooMet[[2]] <- CooMetRev
  return(CooMet)
}




