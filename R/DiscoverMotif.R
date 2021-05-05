
ListFDR = DiscoverMotif <-function(
  FullInputDataDir,
  PathReferenceSequence,
  PathOutput,
  DatasetName,
  Species="homosapiens",
  Growing.Mode="Central",
  Clust.Metric="Cosine",
  Num.CPU=3,
  Prone.Threshold=0.85,
  Resistant.Threshold=0.5,
  Min.MotifLength=5,
  Max.MotifLength=21,
  SeqRepetition.Frequency=3,
  Target.Displacement=0,
  Coord.Displacement=1,
  DrawLogo=FALSE,
  AssignGenes=TRUE,
  ThrFDR = NA,
  MotifTarget = "CG",
  MetHist = FALSE,
  Scan.Type = 'sr',
  Dist.Difference = 'ks',
  Input.Format = "bedmethyl",                                      
  SignalFile = NA,
  lambda = 2,
  beta = 0.00001,
  fdr = 0.05,
  nbins = 1000,
  significance_level = 0.00001,
  cutoff_value = 0.75 ){
  
  # Fibroblast bedmethyl 
  # FullInputDataFile="/mnt/beegfs/german/DMMD_methylation_datasets/fibroblast"
  # IPSC PB3Seq
  # FullInputDataFile="/mnt/beegfs/german/10_Methylation_Call"


  # # TODO: this is an auxiliar path.
  # Fibroblast BedMethyl
  # PathReferenceSequence="/mnt/beegfs/german/REFGENOME/hg38"
  # PB3Seq
  # PathReferenceSequence="/mnt/beegfs/german/REFGENOME/hg19/fasta"
  
  PathOutput=paste(PathOutput, paste("d", Target.Displacement, sep = ""), sep="/")
  # # Draw motifs' logos.
  if ( is.na( ThrFDR ) ) ThrFDR = -Inf
  
  
  ###Parsing
  print("Entered DiscoverMotif")
  
  #Log
  # line <- "Starting logs file"
  LogFilNam <- paste("DMMD2019Disp",Target.Displacement,".logs",sep="") 
  # write(line, file=LogFilNam, append=FALSE)
  # line <- toString(Sys.Date())
  # write(line,file=LogFilNam,append=TRUE)
  
  #Log
  # line <- paste("DISPLACEMENT NUM " , Target.Displacement)
  # write(line,file=LogFilNam,append=TRUE)
  
  #Species
  SPECIES <- c("homosapiens","musmusculus")
  SpeciesVal <- pmatch(tolower(Species), SPECIES)
  if (is.na(SpeciesVal)) stop("invalid species")
  if (SpeciesVal == -1)  stop("ambiguous species")
  
  #Growing Mode
  if (!is.na(pmatch(tolower(Growing.Mode), "central")))  Growing.Mode <- "Central"
  GROWING.MODE <- c("central","latRight","latLeft")
  Growing.ModeVal <- pmatch(tolower(Growing.Mode), GROWING.MODE) 
  if (is.na(Growing.ModeVal)) stop("invalid growing mode")
  if (Growing.ModeVal == -1)  stop("ambiguous growing mode")
  
  #Clustering Metric
  if (!is.na(pmatch(tolower(Clust.Metric), "cosine")))  Clust.Metric <- "Cosine"
  CLUST.METRIC <- c("cosine","pearson")
  Clust.MetricVal <- pmatch(tolower(Clust.Metric), CLUST.METRIC)
  if (is.na(Clust.MetricVal)) stop("invalid clustering metric")
  if (Clust.MetricVal == -1)  stop("ambiguous clustering metric")
  
  #Scan type
  Scan.Type <- tolower(Scan.Type)
  SCAN.TYPE <- c("ss","sr", "tt")
  Scan.TypeVal <- pmatch(Scan.Type, SCAN.TYPE)
  if (is.na(Scan.TypeVal)) stop("invalid scan type")
  if (Scan.TypeVal == -1)  stop("ambiguous scan type")
  
  #Stadistical method to prove the difference between
  #resistant and prone binding distributions in each motif.
  Dist.Difference <- tolower(Dist.Difference)
  DIST.DIFFERENCE <- c("mw","ks")
  Scan.dd <- pmatch(Dist.Difference, DIST.DIFFERENCE)
  if (is.na(Scan.dd)) stop("Invalid distribution difference stadistics method.")
  if (Scan.dd == -1) stop("Ambiguous  distribution difference stadistics method.")
  
  #DatasetName
  if (!is.character(DatasetName)) stop("invalid dataset name. It must be a character")
  
  #NumCPU
  if (!is.numeric(Num.CPU)) stop("invalid number of CPUs. It must be numeric")
  NumCores=parallel::detectCores()
  if (Num.CPU<1 || Num.CPU>NumCores) stop("invalid number of CPUs")
  
  #Methylation Prone Threshold
  if (!is.numeric(Prone.Threshold)) stop("invalid threshold. It must be numeric")
  if (Prone.Threshold>0 & Prone.Threshold<0.5) stop("invalid threshold. It must be between the range [0.5 1]")
  
  #Methylation Resistant Threshold
  if (!is.numeric(Resistant.Threshold)) stop("invalid threshold. It must be numeric")
  if (Resistant.Threshold>0.5 & Resistant.Threshold<1) stop("invalid threshold. It must be between the range [0 0.5]")
  
  # #W_min
  # if (!is.numeric(Min.MotifLength)) stop("invalid minimum length. It must be numeric")
  # if (Min.MotifLength<4 || Min.MotifLength>12) stop("invalid minimum length. It must be between the range [4 12]")
  # 
  # #W_max
  # if (!is.numeric(Max.MotifLength)) stop("invalid maximum length. It must be numeric")
  # if (Max.MotifLength<14 || Max.MotifLength>100) stop("invalid maximum length. It must be between the range [14 100]")
  
  #Frequency
  if(!is.numeric(SeqRepetition.Frequency)) stop("invalid frequency. It must be numeric")
  if (SeqRepetition.Frequency<3 || SeqRepetition.Frequency>6) stop("invalid frequency. It must be between the range [3 6]")
  
  #X
  if (!is.numeric(Target.Displacement)) stop("Invalid value. Target displacement must be numeric")
  
  #DrawLogo and AssignGenes
  if(class(DrawLogo) != "logical")  stop("DrawLogo must be logical")
  if(class(AssignGenes) != "logical") stop("AssignGenes must be logical")
  
  # PATH MANAGEMENT
  ##################################
  # PathInputData
  
  # PathOutput
  if (is.na(PathOutput)){
    # Default path output inside path input 
    PathInputData=dirname(FullInputDataDir)
    PathOutput=file.path(PathInputData,"DMMD_GOUT")
  }
  
  if (dir.exists(PathOutput)) stop(paste("The path to the output directory:",PathOutput," already exists",sep=""))
  DMMD_OutputDir=dir.create(PathOutput,showWarnings = FALSE)
  if (DMMD_OutputDir==FALSE) stop(paste("The path to the output directory:",PathOutput," cannot be created",sep="")) 
  
  PathInputDataFileDir=file.path(PathOutput,DatasetName)
  dir.create(PathInputDataFileDir)
  
  PathFiguresDir=file.path(PathInputDataFileDir,"Figures")
  dir.create(PathFiguresDir)
  
  PathLogosDir=file.path(PathFiguresDir,"Logos")
  dir.create(PathLogosDir)
  
  PathScanDir=file.path(PathFiguresDir,"ScanPlot")
  dir.create(PathScanDir)
  
  PathHtmlDir=file.path(PathInputDataFileDir,"DocHtml")
  dir.create(PathHtmlDir)
  
  PathIndexHtml=file.path(PathInputDataFileDir,"Index.html")
  file.create(PathIndexHtml)
  
  PathForProneHtml=file.path(PathHtmlDir,paste(MotifTarget,"_",DatasetName,"_disp",Target.Displacement,"_p_prone.html",sep=""))
  file.create(PathForProneHtml)
  PathRevProneHtml=file.path(PathHtmlDir,paste(MotifTarget,"_",DatasetName,"_disp",Target.Displacement,"_n_prone.html",sep=""))
  file.create(PathRevProneHtml)
  PathForResisHtml=file.path(PathHtmlDir,paste(MotifTarget,"_",DatasetName,"_disp",Target.Displacement,"_p_resistant.html",sep=""))
  file.create(PathForResisHtml)
  PathRevResisHtml=file.path(PathHtmlDir,paste(MotifTarget,"_",DatasetName,"_disp",Target.Displacement,"_n_resistant.html",sep=""))
  file.create(PathRevResisHtml)
  
  ##### Create Config
  print("Create config")
  NumConfigElem=48
  Config = list()
  
  
  Config$Species=tolower(Species)
  # Added:
  # TODO: either set up debugging or delete this line
  Config$debugging="No"
  Config$LogFile=LogFilNam
  Config$MotifTarget=MotifTarget
  Config$PathOutput=PathOutput
  Config$PathInputDataFileDir=PathInputDataFileDir
  Config$PathFiguresDir=PathFiguresDir
  Config$Logos=PathLogosDir
  Config$Scan=PathScanDir
  Config$HtmlDir=PathHtmlDir
  Config$IndexHtml=PathIndexHtml
  Config$ForProneHtml=PathForProneHtml
  Config$RevProneHtml=PathRevProneHtml
  Config$ForResisHtml=PathForResisHtml
  Config$RevResisHtml=PathRevResisHtml
  
  Config$DirDat = FullInputDataDir
  Config$DirFas = PathReferenceSequence
  
  StrSpecies=switch(tolower(Species),
                    homosapiens=list(RefGenome="hg38",
                                     NumAutosomes=22,
                                     XChromosome="X",
                                     YChromosome="Y",
                                     MChromosome="M",
                                     NumChromosomes=25,
                                     Species="H"),
                    musmusculus=list(RefGenome="mm37",
                                     NumAutosomes=20,
                                     Species="M"))
  
  Config$RefGenome=StrSpecies$RefGenome
  Config$NumAutosomes=StrSpecies$NumAutosomes
  Config$XChr=StrSpecies$XChromosome
  Config$YChr=StrSpecies$YChromosome
  Config$MChr=StrSpecies$MChromosome
  Config$NumChr=StrSpecies$NumChromosomes
  Config$Species=StrSpecies$Species
  Config$Allosomes=c(Config$XChr,Config$YChr,Config$MChr)
  
  Config$DatasetName=DatasetName
  
  Config$nCPU=Num.CPU
  
  Config.GrowingMode=Growing.Mode
  GrowingMode=switch(tolower(Config.GrowingMode), 
                     central = "C", 
                     latRight = "R", 
                     latLeft="L")
  
  Config$GrowingMode=GrowingMode
  
  Config.MetricMode=Clust.Metric
  MetricMode=switch(tolower(Config.MetricMode), 
                    pearson = "P", 
                    cosine = "C")
  
  Config$MetricMode=MetricMode
  
  Config$MethProne = Prone.Threshold
  Config$MethResis = Resistant.Threshold
  Config$w_min = Min.MotifLength
  Config$w_max = Max.MotifLength
  Config$significance_level = significance_level
  Config$lambda = lambda
  Config$fdr= fdr
  Config$ThrFre = SeqRepetition.Frequency
  Config$beta = beta
  Config$nBins = nbins
  Config$X = Target.Displacement
  Config$DrawLogo = DrawLogo
  Config$AssignGenes = AssignGenes
  Config$CooDis = Coord.Displacement
  Config$ThrFDR = ThrFDR
  Config$MetHist = MetHist
  Config$ScnTpe = Scan.Type
  Config$DistDif = Dist.Difference
  Config$InputFormat = Input.Format                                                                
  Config$SignalFile = SignalFile
  Config$cutoff = cutoff_value  
  
  ###Check if any of the Config's elements has been set to NULL
  if (length(Config)<NumConfigElem) stop("A parameter has been incorrectly introduced")
  
  print(Config)
  ####Principal Function
  print("Let's go")
  ListFDR = DNAMethylationMotifFinding(Config)
  return(ListFDR)
}


