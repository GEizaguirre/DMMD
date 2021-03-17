
OptimalDisplacement = DiscoverOptimalMotifDisplacement<-function(Species,
                                                                 Growing.Mode="Central",
                                                                 Clust.Metric="Cosine",
                                                                 InputDataFileName="ffipsc1911",
                                                                 Num.CPU=10,
                                                                 Prone.Threshold=0.95,
                                                                 Resistant.Threshold=0.15,
                                                                 Min.MotifLength=6,
                                                                 Max.MotifLength=15,
                                                                 SeqRepetition.Frequency=5,
                                                                 Target.Displacement=0,
                                                                 Coord.Displacement=1,
                                                                 FullInputDataFile,
                                                                 PathReferenceSequence,
                                                                 PathOutput=NA,
                                                                 VecDisplacements=seq(40,100,1),
                                                                 MotifTarget="CG",
                                                                 MetHist=FALSE,
                                                                 Scan.Type='ss',
                                                                 Dist.Difference='mw',
                                                                 Input.Format="PB3Seq",
                                                                 bedMethylFile="ENCFF232JGO.bed.gz"){
  
  
  # Scan.Type: way of scanning POM.
  #   'ss': Scan POMs against the selected subsets of resistant and prone words. Appropriate when 
  #     bith subsets have similar number of sequences for each length.
  #   'sr': Scan prone POMs against the prone subset and a bunch of random sequences from the resistant
  #     total sequences (same number of sequences as the prone subset). Vice versa for resistant POMs.
  #   'tt': Scan prone and resistant POMs against the total set of prone and resistant sequences. 
  #     Great memory and computational cost.
  


Species <- "homosapiens"
Growing.Mode <- "Central"
Clust.Metric <- "Cosine"
InputDataFileName <- "ffipsc1911"
Num.CPU <- 50
Prone.Threshold <- 0.95
Resistant.Threshold <- 0.15
Min.MotifLength <- 8
Max.MotifLength <- 35
SeqRepetition.Frequency <- 5
Target.Displacement <- 0
Coord.Displacement <- 1
FullInputDataFile="/media/horus/8TB-HD/German/Dropbox/10_Methylation_Call"
# TODO: this is an auxiliar path.
PathReferenceSequence="/media/horus/8TB-HD/German/Dropbox/REFGENOME/hg19/fasta"
# TODO: Original output.
# PathOutput="/home/german/GermanEizaguirre/DMMD-ScriptsAndData/Scripts/GOUT"
# TODO: Ouput for testing.
PathOutput="/media/horus/8TB-HD/German/Dropbox/GOUTforNewScan"
# Displacement in the description file
# VecDisplacements <- seq(1,100,1)
# Original displacement
VecDisplacements <- seq(1,200,20)
MetHist=TRUE
Scan.Type='sr'
Dist.Difference='mw'

  #Parsing
  if(!is.numeric(VecDisplacements)) stop("VecDisplacements must be a numeric vector")
  if(length(VecDisplacements)<2) stop("VecDisplacements must be a numeric vector")

  # Output preparation
  PathOutput = paste(PathOutput, "/DMMD_GOUT", sep="")
  
  #Initialization
  
  ThrFDR=-Inf
  # Draw mmotifs' logos.
  DrawLogo=TRUE
  # Associate motifs with genes.
  AssignGenes=TRUE

  VecMeaFDR<-rep(-Inf,length(VecDisplacements))
  j=1
  
  for(i in VecDisplacements){
    
    #CG.Displacement=i
    ListFDR = DiscoverMotif(Species,
                            Growing.Mode=Growing.Mode,
                            Clust.Metric=Clust.Metric,
                            Num.CPU=Num.CPU,
                            Prone.Threshold=Prone.Threshold,
                            Resistant.Threshold=Resistant.Threshold,
                            Min.MotifLength=Min.MotifLength,
                            Max.MotifLength=Max.MotifLength,
                            SeqRepetition.Frequency=SeqRepetition.Frequency,
                            Target.Displacement=i,
                            Coord.Displacement=Coord.Displacement,
                            DrawLogo=DrawLogo,
                            AssignGenes=AssignGenes,
                            FullInputDataFile,
                            PathReferenceSequence,
                            PathOutput=paste(PathOutput, "Disp", i, sep = ""),
                            ThrFDR=ThrFDR,
                            MotifTarget=MotifTarget,
                            MetHist=MetHist,
                            Scan.Type=Scan.Type,
                            Dist.Difference=Dist.Difference,
                            Input.Format=Input.Format,
                            SignalFile =bedMethylFile)
    
    
    VecFDR=unlist(ListFDR)
    MeaFDR=mean(VecFDR)
    
    #Actualizar
    VecMeaFDR[j]=MeaFDR
    
    ThrFDR=max(ThrFDR,MeaFDR)
    
    j=j+1
  }
  save(VecMeaFDR,file="VecMeaFDR28-08-2019.RData")
  #Calculate optimal value
  OptimalFDR=max(VecMeaFDR)
  iMax=which(VecMeaFDR==max(VecMeaFDR))
  OptimalDisplacement=VecDisplacements[iMax]
  
  return(OptimalDisplacement)
}

