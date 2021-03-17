GenMetHst=function(Config, CooMetFor){
  
  # Generates an histogram of the methylation rate of the annotated CG sites.
  
  # Bind CG sites for all chromosomes.
  HstBnd <- do.call(rbind, CooMetFor)
  HstBndVct <- as.vector(HstBnd[,2])
  # Generate the histogram in a png image.
  HstPth <- file.path(dirname(Config$PathOutput), "CGMethylHistogram.png")
  print(paste("Generating CG methylation histogram at",HstPth))
  png(file=file.path(dirname(Config$PathOutput), "CGMethylHistogram.png"))
  
  # Adjust margins.
  par(mar = c(5, 5, 5, 2.5))
  hist(HstBndVct, breaks=seq(0,1,0.05), 
       main="Histogram of methylation rates of the CG sites.", 
       xlab="Methylation rate",
       ylab="Frequency of CG sites")
  
  not <- dev.off()
}

#CooMov=function(Config,CooMetFor,sense="for"){
CooMov=function(Config,CooMetFor){
  
  # Displaces targets' coordinates to adjust the later
  # calculations.
  
  CooMetForUpd=list()
  #print(paste("DEV ",length(CooMetFor)))
  for(i in 1:Config$NumChr){
    print(i)
    CooMetForMat=CooMetFor[[i]]
    # Sum Config$CooDis positions to each target's coordinate.
    # (Config$CooDis is 1 by default).
    CooMetForMat[c("ColCoo")] <- lapply(CooMetForMat[c("ColCoo")], function(x) x+ Config$CooDis)
    #if (sense=="for") CooMetForMat[c("ColCoo")] <- lapply(CooMetForMat[c("ColCoo")], function(x) x+ Config$CooDis)
    #else CooMetForMat[c("ColCoo")] <- lapply(CooMetForMat[c("ColCoo")], function(x) if (x-Config$CooDis > 0) {x-Config$CooDis} else {1})
    CooMetForUpd[[i]]=CooMetForMat
  }
  return(CooMetForUpd)
}

ReadFasta=function(Config){
  
  # Reads chromosomes' fasta files and saves the sequences in
  # a R structure.
  
  SeqChr=list()
  
  for(i1 in 1:Config$NumAutosomes){
   
    NamFilFas=paste("chr",i1,".fa",sep="")
    # print(NamFilFas)
    FulNamFilFas=file.path(Config$DirFas,NamFilFas)
    seqs=read.fasta(FulNamFilFas,as.string=TRUE)
    SeqChr[[i1]]=seqs[[1]]
  }
  
  for(i2 in Config$Allosomes){
    
    i1=i1+1
    NamFilFas=paste("chr",i2,".fa",sep="")
    # print(NamFilFas)
    FulNamFilFas=file.path(Config$DirFas,NamFilFas)
    seqs=read.fasta(FulNamFilFas,as.string=TRUE)
    SeqChr[[i1]]=seqs[[1]]
  }
  return(SeqChr)
}
 
 

DicWord=function(Config,CooMet,SeqChr,sense="for"){
  
  # Creates  [word, appaerance indexes] dictionaries based on
  # the annotation files and the chromosome FASTA files.
  # Indexes point to the lines of the coordinate files created
  # at the C function CooChr.

  DicWordRes=list()
  # SeqMetFreW will contain grouped sequences with methylation rates and filtered by frequency.
  DicWordRes$SeqMetFreW=list()
  # SeqMetTot will contain all sequences with their methylation rates, classified as 'prone' or
  # 'resistant' to methylation.
  DicWordRes$SeqTot=list()
  DicWordRes$SeqTot$SeqPrn=list()
  DicWordRes$SeqTot$SeqRes=list()
  
  #temporaldir=tempdir()
  temporaldir=Config$PathOutput
  CooFilDir=paste(temporaldir,"/CoordinateFiles",sep="")
  if(file.exists(CooFilDir)==FALSE) dir.create(CooFilDir)
  
  # For each possible length of word.
  for(w in Config$w_min:Config$w_max){
    
    DicWord=list()
    MetRati=list()
    CooVec=list()
    ChrVec=list()
  
    # For each chromosome.
    for(i1 in 1:Config$NumChr){
    
      # Get targets' coordinates and sequence for each chromosome.
      CooMetMat = CooMet[[i1]]
      SeqChrMat = SeqChr[[i1]]
      LenSeqChrMat = nchar(SeqChrMat)
      
      # Get total number of coordinates.
      nrowCooMetMat=nrow(CooMetMat)
      print(paste("Chr: ",i1,", Length: ",2*w+2, sep=""))
      # SeqDicMet will contain words in its first column and methylation frequency of the target related to
      # each word in the second column.
      if (sense=="for") SeqDicMet = .Call("SeqDic",CooMetMat,w,nrowCooMetMat,SeqChrMat,LenSeqChrMat,Config$GrowingMode,Config$X,1)
      else SeqDicMet = .Call("SeqDic",CooMetMat,w,nrowCooMetMat,SeqChrMat,LenSeqChrMat,Config$GrowingMode,-Config$X,0)
        
      #print(nrowCooMetMat)
      #print(length(SeqDicMet[[1]]))
      
      # Delete words that exceeded the limit.
      WrdMat <- matrix(SeqDicMet[[1]])
      WrdMatLst <- as.list(WrdMat[,1])
      IndOutBnd <- as.integer(c(which(WrdMatLst=="no"), which(WrdMatLst=="rval")))
      if (length(IndOutBnd)!=0) { WrdMatLst=WrdMatLst[-IndOutBnd] }
      WrdMat <- as.matrix(as.character(WrdMatLst))
      
      # Delete methylation rates for these words.
      
      MetMat <- matrix(SeqDicMet[[2]])
      MetMatLst <- as.list(MetMat[,1])
      if (length(IndOutBnd)!=0) { MetMatLst=MetMatLst[-IndOutBnd] }
      MetMat <- as.matrix(as.double(MetMatLst))
      
      # Obtained words.
      DicWord[[i1]] = WrdMat
      # Methylation rates of the targets associated to the words.
      MetRati[[i1]] = MetMat
      # Coordinates of the targets associated to the words.
      if (length(IndOutBnd)!=0) { CooVec[[i1]] = matrix(CooMetMat[[1]][-IndOutBnd]) }
      else { CooVec[[i1]] = matrix(CooMetMat[[1]]) }
      # Structure to create the column indicating the chromosome.
      #ChrVec[[i1]] = matrix(rep(i1,(nrowCooMetMat-length(IndOutBnd))))
      ChrVec[[i1]] = matrix(rep(i1,length(WrdMat)))
      #print(paste("dicword",nrow(WrdMat),"nrowCooMetMat",nrowCooMetMat,"ChrVec[[]]",nrow(ChrVec[[i1]])))
    }
    
    
    # Bind the content of each of the four structures
    # into a unique column.
    MatSeqColl = do.call(rbind,DicWord)
    Methyl = do.call(rbind,MetRati)
    Coo = do.call(rbind,CooVec)
    Chr = do.call(rbind,ChrVec)

    # save(MatSeqColl, file=paste("testDicWordCluster",Config$X,sep=""))     
    #print(paste("MatSeqColl:", nrow(MatSeqColl), "Methyl:", nrow(Methyl),"Coo:",nrow(Coo),"Chr:",nrow(Chr)))
    
    #print(paste("MatSeqColl:", nrow(MatSeqColl), "Methyl:", nrow(Methyl),"Coo:",nrow(Coo),"Chr:",nrow(Chr)))
    # Create a unique information table.
    MatSeqMet = data.table(MatSeqColl,Methyl,Coo,Chr)
    setnames(MatSeqMet, c("V1","V2","V3","V4"))
    # Set words as the key for the dictionary.
    setkey(MatSeqMet,"V1")
   
    # Separate again the structure into columns.
    MatSeqCollapse=MatSeqMet$V1
    MethylationOrdered=MatSeqMet$V2
    CooOrdered=MatSeqMet$V3
    ChrOrdered=MatSeqMet$V4

    #print(object.size(MatSeqCollapse))   
    #print(nrow(MatSeqCollapse))

    SeqMetCooChr = data.table(MatSeqCollapse,MethylationOrdered)
    #save(SeqMetCooChr, MatSeqCollapse, MethylationOrdered, file=paste("testDicWordCluster",Config$X,sep=""))
    # Group the word and methylation rate  pairs' table by words (get number of instances of each word).
    SeqFre = SeqMetCooChr[, .N ,by = MatSeqCollapse]
    SeqFre = as.data.frame(SeqFre)
    
    #print(paste("nrow Seqfre:",nrow(SeqFre),"length SeqFre:",length(SeqFre),"size SeqFre:",object.size(SeqFre)))
 
    #print(paste("before frequency threshold", nrow(SeqFre[,1])))
    # Select those elements from the data frame with methylation rate above the threshold (Get original indexes .
    IndFreAboThr=which(SeqFre[,2]>=Config$ThrFre)
    NumSelSeq=length(IndFreAboThr)
    #print(paste("selected sequences",length(IndFreAboThr)))
    SeqMetCooChr = as.data.frame(MatSeqMet)

    CumSum=cumsum(SeqFre[,2])
    SeqMetCooChr[,3]=as.numeric(SeqMetCooChr[,3])
    SeqMetCooChr[,4]=as.numeric(SeqMetCooChr[,4])
    # Create coordinate files for the words.
    .Call("CooChr",SeqMetCooChr,CumSum,w,CooFilDir)
    
    # Join words and methylation rates.
    SeqMetMea <- structure(list(Seq = MatSeqCollapse,
                                Methyl = MethylationOrdered),
                           .Names = c("Seq", "Methyl"), class = c("data.table","data.frame"))
    
    # Calculate mean methylation rate for each sequence.
    colstoavg <- names(SeqMetMea[,2,with=F])
    SeqMetMea.mean <- SeqMetMea[,lapply(.SD,mean,na.rm=TRUE),by=SeqMetMea$Seq,.SDcols=colstoavg]
    
    # Join grouped sequences with mean methylation rates.
    SeqMetFre=as.data.frame(cbind(SeqMetMea.mean,SeqFre[,2]))
    names(SeqMetFre)[3]="Freq"
    names(SeqMetFre)[1]="Seq"
  
    SeqMetFreFilt=SeqMetFre[which(SeqMetFre[,3] >= Config$ThrFre), ]
    SeqMetFreFilt=cbind(SeqMetFreFilt,IndFreAboThr)
    names(SeqMetFreFilt)[4]="Index"
    
    #print(paste("nrow SeqMetFreFilt:", nrow(SeqMetFreFilt)))
    #print(format(object.size(SeqMetFreFilt), units="Gb"))
    #save(SeqMetFreFilt, file=paste("testDicWordCluster",Config$X,sep=""))
    DicWordRes$SeqMetFreW[[w]]=SeqMetFreFilt
    #print(format(object.size(DicWordRes), units="Gb"))
    
    print(paste(NumSelSeq,"sequences selected for equal representation."))
    # Select methylation prone/resistant sequences for scanning.
#    IndMetPrn <- which(Methyl>=Config$MethProne)
#    IndMetRes <- which(Methyl<=Config$MethResis)
#    print(paste("Total sequences for prone:",length(IndMetPrn),"Total sequences for resistant:",length(IndMetRes)))
#    SeqPrn <- MatSeqColl[IndMetPrn]
#    SeqRes <- MatSeqColl[IndMetRes]
#    DicWordRes$SeqTot$SeqPrn[[w]] <- sample(SeqPrn, NumSelSeq)
#    DicWordRes$SeqTot$SeqRes[[w]] <- sample(SeqRes, NumSelSeq)
    
    PrnSelectedIndexes = which(SeqMetFreFilt$Methyl >= Config$MethProne )
    PrnSeqs = SeqMetFreFilt$Seq[PrnSelectedIndexes]
    PrnFreqs = SeqMetFreFilt$Freq[PrnSelectedIndexes]
    PrnSeqsMult <- foreach(i=1:length(PrnSeqs), .combine=c) %do% {
      rep(PrnSeqs[i], PrnFreqs[i])
    }
  
    DicWordRes$SeqTot$SeqPrn[[w]] <- PrnSeqsMult
    
    ResSelectedIndexes = which(SeqMetFreFilt$Methyl <= Config$MethResis )
    ResSeqs = SeqMetFreFilt$Seq[ResSelectedIndexes]
    ResFreqs = SeqMetFreFilt$Freq[ResSelectedIndexes]
    ResSeqsMult <- foreach(i=1:length(ResSeqs), .combine=c) %do% {
      rep(ResSeqs[i], ResFreqs[i])
    }
    
    DicWordRes$SeqTot$SeqRes[[w]] <- ResSeqsMult


  
    #print(object.size(SeqPrn))
    #print(object.size(SeqRes))
    #print(paste("IndMetPrn:",length(IndMetPrn),"IndMetRes:",length(IndMetRes),"SeqPrn",length(SeqPrn),"SeqRes:",length(SeqRes)))
    rm(SeqPrn)
    rm(SeqRes)
    gc()
  }
  return(DicWordRes)
}

Rev=function(Config,SeqMetFreW){
  
  for(w in Config$w_min:Config$w_max){

    Seq=SeqMetFreW[[w]]$Seq
    Methyl=SeqMetFreW[[w]]$Methyl
    Freq=SeqMetFreW[[w]]$Freq
    Index=SeqMetFreW[[w]]$Index
    
    RevSeq=.Call("Reverse",Seq,w)
    
    SeqMetFre=data.frame(RevSeq,Methyl,Freq,Index)
    SeqMetFre$RevSeq=as.character(SeqMetFre$RevSeq)
    names(SeqMetFre)[1]="Seq"
    SeqMetFreW[[w]]=SeqMetFre
  }
  return(SeqMetFreW)
}

RevTot=function(Config,Seqs){
  
  for(w in Config$w_min:Config$w_max){
    

    Seq=Seqs[[w]]
    
    RevSeq=.Call("Reverse",Seq,w)
    
    Seqs[[w]]=as.character(RevSeq)
  }
  return(Seqs)
}

DelGaps=function(Config,SeqMetFreW){
  
  # Delete words with gaps.
  # Gaps are represented with "n" characters.
  
  Gap="n"
  
  for(w in Config$w_min:Config$w_max){
    
    Seq=SeqMetFreW[[w]]$Seq
    Methyl=SeqMetFreW[[w]]$Methyl
    Freq=SeqMetFreW[[w]]$Freq
    Index=SeqMetFreW[[w]]$Index
    # Detect positions of words that contain "n".
    IndFou=grepl(Gap, Seq, ignore.case = TRUE)
    IndRm=which(IndFou==TRUE)
    
    # Delete detected positions.
    if(length(IndRm)!=0) {
      Seq=Seq[-IndRm]
      Methyl=Methyl[-IndRm]
      Freq=Freq[-IndRm]
      Index=Index[-IndRm]
      SeqMetFre=data.frame(Seq,Methyl,Freq,Index)
      SeqMetFre$Seq=as.character(SeqMetFre$Seq)
      SeqMetFreW[[w]]=SeqMetFre
    }
  }
  return(SeqMetFreW)
}

DelGapsTot=function(Config,Seqs){
  
  # Delete words with gaps.
  # Gaps are represented with "n" characters.
  
  Gap="n"
  
  for(w in Config$w_min:Config$w_max){
    
    Seq=Seqs[[w]]
    # Detect positions of words that contain "n".
    IndFou=grepl(Gap, Seq, ignore.case = TRUE)
    IndRm=which(IndFou==TRUE)
    
    # Delete detected positions.
    if(length(IndRm)!=0) {
      Seq=Seq[-IndRm]
      Seq=as.character(Seq)
      Seqs[[w]]=Seq
    }
  }
  return(Seqs)
}

ProneMet=function(Config,SeqMetFreW){
  
  # Select those words with methylation rate above
  # the desired methylation proneness threshold.
  
  SeqMetFreProne=list()
  
  for(w in Config$w_min:Config$w_max){
    
    SeqMetFreProne[[w]]=SeqMetFreW[[w]][which(SeqMetFreW[[w]][,2] >= Config$MethProne), ]
  }
  return(SeqMetFreProne)
}

ResisMet=function(Config,SeqMetFreW){
  
  # Select those words with methylation rate under
  # the desired methylation resistance threshold.
  
  SeqMetFreResis=list()
  
  for(w in Config$w_min:Config$w_max){
    
    SeqMetFreResis[[w]]=SeqMetFreW[[w]][which(SeqMetFreW[[w]][,2] <= Config$MethResis), ]
  }
  return(SeqMetFreResis)
}

FreqVec=function(Config,SeqMetFre){
  
  FreVecW=list()
  
  for(w in Config$w_min:Config$w_max){
    
    SeqMetFreMat=SeqMetFre[[w]]
    FreVec=SeqMetFreMat[,3]
    FreVec=matrix(FreVec,length(FreVec),1)
    FreWVec=matrix(0,length(FreVec),2*w+2)
    
    if (length(FreWVec)>0){FreWVec=t(sapply(1:length(FreVec),function(i) rep(FreVec[i,],2*w+2)))}
    else FreWVec = NULL;
    
    
    FreVecW[[w]]=FreWVec 
  }
  return(FreVecW)
}

######POM


max_word_limit = 46000

ReduceWords=function(Config, SeqMetFreWFreVecW){
  SeqMetFreWFreVecWNew = list()
  for (w in Config$w_min:Config$w_max){
     if (length(SeqMetFreWFreVecW[[w]]$SeqW)>max_word_limit){
        freq_list = SeqMetFreWFreVecW[[w]]$X1
        ind_out <- head(sort.list(freq_list), n=(length(SeqMetFreWFreVecW[[w]]$seqW)-max_word_limit))
        SeqMetFreWFreVecWNew[[w]] <- SeqMetFreWFreVecW[[w]][-ind_out,]
     }
    else {
      SeqMetFreWFreVecWNew[[w]] <- SeqMetFreWFreVecW[[w]]
    }
  }
  return(SeqMetFreWFreVecWNew)
}



CreatePom=function(SeqSep,FreWVec,w){
  
  PomList=list()
  
  for(i in 1:nrow(SeqSep)){
    
    POM = matrix(0, 4, 2*w+2)
    
    for(j in 1:ncol(SeqSep)){
      
      POM[SeqSep[i,j],j]=FreWVec[i,j]
    }
    
    PomList[[i]]=POM         
  }    
  
  return (PomList)
}
POM=function(Config,SeqMetFreWFreVecW){
  
  PomW=list()
  
  for(w in Config$w_min:Config$w_max){
    
    print(w)
   
    if (!is.null(SeqMetFreWFreVecW[[w]]) && (length(SeqMetFreWFreVecW[[w]]$MetW)>0)){
    SeqMetFreFreVec=SeqMetFreWFreVecW[[w]]
    FreWVec=SeqMetFreFreVec[,5:ncol(SeqMetFreFreVec)]
    FreWVec=setNames(FreWVec,rep(" ",length(FreWVec)))
    FreWVec=as.matrix(FreWVec)
    
    SeqMetFreFreVec=data.frame(lapply(SeqMetFreFreVec, as.character), stringsAsFactors=FALSE)
    
    Split=strsplit(SeqMetFreFreVec$SeqW,NULL,fixed=TRUE)
    
    SeqSep=t(sapply(1:length(Split), function(x) unlist(Split[[x]])))
    
    SeqSep[SeqSep == "a"] = 1
    SeqSep[SeqSep == "c"] = 2
    SeqSep[SeqSep == "g"] = 3
    SeqSep[SeqSep == "t"] = 4
    
    class(SeqSep)='numeric'

    
    pom=CreatePom(SeqSep,FreWVec,w)
    PomW[[w]]=pom
    }
  }
  return(PomW)
}

POMVec=function(Config,PomW){
  
  PomWVec=list()
  
  for(w in Config$w_min:Config$w_max){
    
    if (!is.null(PomW[[w]])) {
    
    Pom=PomW[[w]]
    PomVec=matrix(0,length(Pom),4*(2*w+2))
    
    PomVec=t(sapply(1:length(Pom),function(i) as.vector(Pom[[i]])))
    }
    
    PomWVec[[w]]=PomVec
  } 
  return(PomWVec)
}

PomClu=function(Config,PomWVec,SeqMetFreWFreVecW){
  
  PomClu=list()
  CutOff=seq(0.05,0.95,0.05)
  
  
  for(w in Config$w_min:Config$w_max){
    
    # #Log
    # line <- paste("PomClu for ", w)
    # write(line,file=Config$LogFile,append=TRUE)
    
    if (!is.null(SeqMetFreWFreVecW[[w]])){
      
      if (length(SeqMetFreWFreVecW[[w]]$IndW)!=1){
        SeqMetFreIndFreVec=SeqMetFreWFreVecW[[w]]$IndW
        print(paste(length(SeqMetFreIndFreVec), "indexes"))
        SilVec=rep(0,length(CutOff))
        i3=1
        
        PomMat=PomWVec[[w]]     
        pommat=lapply(seq_len(nrow(PomMat)), function(i) PomMat[i,])
        nrow=length(pommat)
      
        Len=4*(2*w+2)
        
        Distance=.Call("DissimilarityMatrix", pommat, nrow, Len, Config$MetricMode)
      
        # #Log
        # line <- paste("DissimilarityMatrix done for ", w)
        # write(line,file=Config$LogFile,append=TRUE)
        
        size=nrow
        DistMat=vec2dist(Distance, size)
        rm(Distance)
        gc()
        hcl=hclust(DistMat,method="complete")
        MaxHei=max(hcl$height)
        MinHei=min(hcl$height)
   	
      	if ( is.na(Config$cutoff) ){
      		for(c in CutOff){
        
          	cut=(MaxHei-MinHei)*c + MinHei
          	CutTre=cutree(hcl, h = cut)
          	NumUni=length(unique(CutTre))
          
          	if(NumUni >= 2 & NumUni <= nrow-1){
            		SilCoef=silhouette(CutTre, DistMat)
            		SilAvgWid=summary(SilCoef)$avg.width
          	}
          	else{
            		SilAvgWid=-1
          	}
          
          	SilVec[i3]=SilAvgWid
          	i3=i3+1
      		}
      		MaxSilIdx=which(SilVec==max(SilVec))
            		CutPoint=CutOff[MaxSilIdx]
      		CutPointDef=max(CutPoint)
      
      	}
      	else {
      		CutPointDef = Config$cutoff
      	}
        #save(SilVec, file=paste("SilhouetteVec",w,"17-07-2019-0117-07-2019-01.RData", sep = ""))
        rm(DistMat)
        gc()
  
        CutHei=(MaxHei-MinHei)*CutPointDef + MinHei
        CutTreDef=cutree(hcl, h = CutHei)
        NumClu=max(CutTreDef)
        rm(hcl)
        gc()
        Idx=rep(0,NumClu)
        POM=list()
        Ind=list()
      
        for (c in 1:NumClu){
          
          IndEleClu=which(CutTreDef==c)
          
          if(length(IndEleClu)==1) Idx[c]=c
          
          # Sum the POM of the members of each cluster
          Mat=rbind(PomMat[IndEleClu,])
          PomSums=colSums(Mat)
          PomSum=matrix(PomSums,nrow=4,ncol=2*w+2)
          POM[[c]]=PomSum
          
          IndSum=rbind(SeqMetFreIndFreVec[IndEleClu])
          Ind[[c]]=IndSum
        }
    
        IdxRm=which(Idx!=0)
        if (length(IdxRm)>0){
          POM=POM[-IdxRm]
          Ind=Ind[-IdxRm]
        }
        
        # Control for motif lengths with only 1 POM (not very probable)
      } else if (length(SeqMetFreWFreVecW[[w]]$IndW)==1) {
        POM=list()
  
        mat<-matrix(vector(), nrow=4, ncol=w*2+2)
        for (i in seq(1, (w*2+2)*4, 4)){
          mat[i%%4, i%/%4+1] = PomWVec[[w]] [1, i]
          mat[i%%4+1, i%/%4+1] = PomWVec[[w]] [1, i+1]
          mat[i%%4+2, i%/%4+1] = PomWVec[[w]] [1, i+2]
          mat[i%%4+3, i%/%4+1] = PomWVec[[w]] [1, i+3]
        }
  
        POM[[1]] = mat
        Ind=list()
        Ind[1] = SeqMetFreWFreVecW[[w]]$IndW[1] 
      }
      PomClu[[w]]=list(POM,Ind)
    }
  }
  return(PomClu)
}

GetLen=function(Config,PomW){
  
  Len=rep(0,(Config$w_max-Config$w_min))  
  
  for(w in Config$w_min:Config$w_max){
   if (!is.null(PomW[[w]])){
   if(length(PomW[[w]][[1]]) != 0 ) Len[w]=w
   }
    else {
      Len[w] = 0
    }
  }
 Len=which(Len!=0)
 return (Len)
}

PMV=function(Config,PomW,LenMotif){
  
  PMV=list()

  for(w in LenMotif){
  
    if (!is.null(PomW[[w]])){
    PmvW=list()
    
    for(i in 1:length(PomW[[w]][[1]])){
      PmvW[[i]] = apply(PomW[[w]][[1]][[i]],2,max)
    }
    PMV[[w]]=PmvW
    }
  }
  return(PMV)
}

Normalization=function(x){
  
  Tot = sum(x)
  return(x / Tot); 
}

PWM=function(Config,PomW,LenMotif){
  
  PWM=list()
  
  for(w in LenMotif){
    
    if (!is.null(PomW[[w]])){
    
    Pwm=list()
    
    for(i in 1:length(PomW[[w]][[1]])){
      Pwm[[i]] = apply(PomW[[w]][[1]][[i]], 2, Normalization) 
      rownames(Pwm[[i]])[1:4]<-c("A","C","G","T") 
    } 
    
    PWM[[w]]=Pwm
    }
  }
  return(PWM)
}

WeiPOM=function(Config,PomW,LenMotif){
  
  WeiPOM=list()
  
  for(w in LenMotif){
    
    if (!is.null(PomW[[w]])){
    Pom=PomW[[w]][[1]]
    Weipom=list()
    
    for(i in 1:length(Pom)){
      SumPOM=apply(Pom[[i]],2,sum)
      SumSumPOM=sum(SumPOM)
      Weipom[[i]]=SumPOM/SumSumPOM
    }
    
    WeiPOM[[w]]=Weipom
    }
  }
  return(WeiPOM)
}

WeiLogPMV=function(Config,WeiPOM,PMV,LenMotif){
  
  WeiLogPMV=list()
  
  for(w in LenMotif){
    
    if (!is.null(WeiPOM[[w]]) && !is.null(PMV[[w]])){
    WeiPom=WeiPOM[[w]]
    Pmv=PMV[[w]]
    WeiLogpmv=list()
    
    for(i in 1:length(Pmv)){
      WeiLogpmv[[i]]=WeiPom[[i]]*log(Pmv[[i]]+ Config$beta)  
    }
    WeiLogPMV[[w]]=WeiLogpmv
    }
  }
  return(WeiLogPMV)
}

WeiLogPOM=function(Config,WeiPOM,POM,LenMotif){
  
  WeiLogPOM=list()
  
  for(w in LenMotif){
    
    if (!is.null(WeiPOM[[w]]) && !is.null(POM[[w]])){
    WeiPom=WeiPOM[[w]]
    Pom=POM[[w]][[1]]
    WeiLogpom=list()
    
    for(i in 1:length(WeiPom)){
      WeiLogpom[[i]]=WeiPom[[i]]*log(Pom[[i]]+ Config$beta)  
    }  
    WeiLogPOM[[w]]=WeiLogpom
    }
  }
  return(WeiLogPOM)
}

###########
IndBasSeq=function(SeqMetFreW){

  #
  # IndBasSeq
  # Auxiliary function. Passes a frame of nucleotidic sequences into
  # numeric format.
  # a -> 0
  # c -> 1
  # g -> 2
  # t -> 3
  
  
  Split=strsplit(SeqMetFreW$Seq,NULL,fixed=TRUE)
  SeqSep=t(sapply(1:length(Split), function(x) unlist(Split[[x]])))
  
  SeqSep[SeqSep == "a"] = 0
  SeqSep[SeqSep == "c"] = 1
  SeqSep[SeqSep == "g"] = 2
  SeqSep[SeqSep == "t"] = 3
  
  class(SeqSep)='numeric'
  
  SeqSep=lapply(seq_len(nrow(SeqSep)), function(i) SeqSep[i,])
  
  return(SeqSep) 
}

IndBasSeqTot=function(Seqs){
  
  #
  # IndBasSeq
  # Auxiliary function. Passes a frame of nucleotidic sequences into
  # numeric format.
  # a -> 0
  # c -> 1
  # g -> 2
  # t -> 3
  
  
  Split=strsplit(Seqs,NULL,fixed=TRUE)
  SeqSep=t(sapply(1:length(Split), function(x) unlist(Split[[x]])))
  
  SeqSep[SeqSep == "a"] = 0
  SeqSep[SeqSep == "c"] = 1
  SeqSep[SeqSep == "g"] = 2
  SeqSep[SeqSep == "t"] = 3
  
  class(SeqSep)='numeric'
  
  SeqSep=lapply(seq_len(nrow(SeqSep)), function(i) SeqSep[i,])
  
  return(SeqSep) 
}

ScanSeqs=function(WeiLogPOM,WeiLogPWV,NumSeq,LenMot) {
  
  #
  # Calculates the affinity of each word from a group of them
  # to a motif of the same length according to the POM of the motif.
  # WeiLogPOM: POM (Position Ocurrence Matrix) of the motif.
  # WeiLogPMV: PMV (Position Weight Matrix)  of the motif - POM's normalization.
  # NumSeq: Bunch of sequences to be analyzed.
  # LenMot: Lenght of the sequences and the motif (must be the same).
  

  # Resulting structure with scores results for each motif.
  
  ScoRes = rep(NA, length(NumSeq))
  
  nSeq <- length(NumSeq)
  
  WeiLogPomElem <- WeiLogPOM[[1]]
  WeiLogPWVElem <- WeiLogPWV[[1]]
  
  if (nSeq>0){
  for (i1 in 1:nSeq) {
    
    NumSeqElem <- NumSeq[[i1]]
    Sco <- 0.0
    
    for (i2 in 1:LenMot) {
      
      # Get the nucleotide in the current position of the sequence.
      IndBas <- NumSeqElem[i2] + 1
      # Get the score for the nucleotide in the motif.
      Num <- WeiLogPomElem[IndBas,i2]
      # Get the normalization value for the position.
      Den <- WeiLogPWVElem[i2]
      CurSco <- Num - Den
      Sco <- Sco + CurSco
      
    }
    
    ScoRes[i1] <- Sco
    
  }
  }
  
  class(ScoRes) <- 'numeric'
  return (ScoRes)
}

ScanSeqsFreq=function(WeiLogPOM,WeiLogPWV,NumSeq,LenMot, FreqVec) {
  
  #
  # Calculates the affinity of each word from a group of them
  # to a motif of the same length according to the POM of the motif.
  # WeiLogPOM: POM (Position Ocurrence Matrix) of the motif.
  # WeiLogPMV: PMV (Position Weight Matrix)  of the motif - POM's normalization.
  # NumSeq: Bunch of sequences to be analyzed.
  # LenMot: Lenght of the sequences and the motif (must be the same).
  
  # Resulting structure with scores results for each motif.
  
  ScoRes = rep(NA, sum(FreqVec))
  
  nSeq <- length(NumSeq)
  
  WeiLogPomElem <- WeiLogPOM[[1]]
  WeiLogPWVElem <- WeiLogPWV[[1]]
  
  PosCounter <- 1
  if (nSeq>0){
  for (i1 in 1:nSeq) {
    
    NumSeqElem <- NumSeq[[i1]]
    Sco <- 0.0
    
    for (i2 in 1:LenMot) {
      
      # Get the nucleotide in the current position of the sequence.
      IndBas <- NumSeqElem[i2] + 1
      # Get the score for the nucleotide in the motif.
      Num <- WeiLogPomElem[IndBas,i2]
      # Get the normalization value for the position.
      Den <- WeiLogPWVElem[i2]
      CurSco <- Num - Den
      Sco <- Sco + CurSco
      
    }
    
    #ScoRes[i1] <- Sco
    #ScoRes <- c(ScoRes, rep(Sco, FreqVec[i1]))
    ScoRes[PosCounter:(PosCounter + FreqVec[i1] -1)] <- Sco
    PosCounter <- PosCounter + FreqVec[i1]
  }
  }
  
  class(ScoRes) <- 'numeric'
  return (ScoRes)
}

ScanSeqsPar=function(WeiLogPOM,WeiLogPWV,NumSeq,LenMot) {
  
  #
  # Calculates the affinity of each word from a group of them
  # to a motif of the same length according to the POM of the motif.
  # WeiLogPOM: POM (Position Ocurrence Matrix) of the motif.
  # WeiLogPMV: PMV (Position Weight Matrix)  of the motif - POM's normalization.
  # NumSeq: Bunch of sequences to be analyzed.
  # LenMot: Lenght of the sequences and the motif (must be the same).
  
  # Resulting structure with scores results for each motif.
  
  nSeq <- length(NumSeq)
  
  WeiLogPomElem <- WeiLogPOM[[1]]
  WeiLogPWVElem <- WeiLogPWV[[1]]
  registerDoMC(Config$nCPU)
  
  ScoRes = foreach (i1=1:nSeq, .combine=c) %dopar%{
    
    NumSeqElem <- NumSeq[[i1]]
    Sco <- 0.0
    
    for (i2 in 1:LenMot) {
      
      # Get the nucleotide in the current position of the sequence.
      IndBas <- NumSeqElem[i2] + 1
      # Get the score for the nucleotide in the motif.
      Num <- WeiLogPomElem[IndBas,i2]
      # Get the normalization value for the position.
      Den <- WeiLogPWVElem[i2]
      CurSco <- Num - Den
      Sco <- Sco + CurSco
      
    }
    
    #ScoRes[i1] <- Sco
    Sco
    
  }
  
  #class(ScoRes) <- 'numeric'
  return (ScoRes)
}

intervals=function(Pts, TtlLng){
  
  # Generates a list of intervals from a total value.
  # Pts: Number of intervals.
  # TtlLng: Total value.
  
  RngSiz=round(TtlLng/Pts,0)
  Res=list()
  PrvEnd<-0
  for(i in 1:Pts){
    Str<-round((i-1)*RngSiz+1,0)
    if (i==Pts) End<-TtlLng
    else End<-round((i)*RngSiz,0)
    if (Str== PrvEnd) Str=Str+1
    
    # Each interval's start will be in the first column
    # and the end in the second column.
    Res<-rbind(Res, c(Str,End))
    PrvEnd<-End
  }
  return(Res)
}

SelSamples=function(PopSeqs, RefSeqs, Lengths){
  
  # Generates random samples of sequences from the population structure,
  # The size of the sample is the same as the number of sequences for each length 
  # in the reference sequences' structure.
  
  ResSmp <- list()
  for (i in Lengths){
    ResSmp[[i]] <- sample(PopSeqs[[i]], min(length(PopSeqs[[i]]), sum(RefSeqs[[i]][[3]])))
  }
  return(ResSmp)
  
}

SelSamples2=function(RefSeqs_match, RefSeqs_unmatch, Lengths){
  
  # Generates random samples of sequences from the population structure,
  # The size of the sample is the same as the number of sequences for each length 
  # in the freference sequences' structure.
  
  ResSmp <- list()
  for (i in Lengths){
    
	if ( length( RefSeqs_match[[i]] ) <= length( RefSeqs_unmatch[[i]] ) ){
		ResSmp[[i]] <- RefSeqs_match[[i]]
	}
	else {
		ResSmp[[i]] <- sample( RefSeqs_match[[i]], length( RefSeqs_unmatch[[i]] ) )
	}
    
  }
  return(ResSmp)
  
}

#############

Scan=function(Config,SeqMetFreW,WeiLogPMV,WeiLogPOM,LenMotif){
  
  # Result.
  ScoPerLen=list()
  
  for(w in LenMotif){
    
    if (!is.null(SeqMetFreW[[w]]) && !is.null(WeiLogPOM[[w]]) && !is.null(WeiLogPMV[[w]])){
    #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    
    #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    
    # All sequences for the current length.
    SeqMetFre=SeqMetFreW[[w]]
    # Sequences are passed to numeric format.
    NumSeq=IndBasSeq(SeqMetFre)
    LenMot=2*w+2
    ScoDisPerMot=list()
    
    for(i in 1:length(WeiLogPOM[[w]])){
      
      ScoDis=ScanSeqs(WeiLogPOM[[w]][i],WeiLogPMV[[w]][i],NumSeq,LenMot)
      
      # Scan sequences against POM in cpp
      ScoDis=scan_seqs_c(length(NumSeq), LenMot, NumSeq, WeiLogPOM[[w]][i], WeiLogPMV[[w]][i])
      
      class(ScoDis) <- 'numeric'
      
      bincounts=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$counts
      binbreaks=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$breaks
      binbreaks=binbreaks[1:(length(binbreaks)-1)]
      ScoDatFram=data.frame(bincounts,binbreaks)
      ScoDisPerMot[[i]]=ScoDatFram
      gc()
    }
    ScoPerLen[[w]]=ScoDisPerMot
    }
  }
  return (ScoPerLen)
  
}

ScanPar=function(Config,SeqMetFreW,WeiLogPMV,WeiLogPOM,LenMotif){
  
  # 
  # Parallel Scanning using only R code.
  # 
  
  ScoPerLen=list()
  registerDoMC(Config$nCPU)
  
  # For each length.
  for (w in LenMotif) {
    
    if (!is.null(SeqMetFreW[[w]]) && !is.null(WeiLogPOM[[w]]) && !is.null(WeiLogPMV[[w]])){
    # #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    # 
    # #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    
    # Every sequence for this length.
    SeqMetFre=SeqMetFreW[[w]]
    # Sequences passed to numerical format.
    NumSeq=IndBasSeq(SeqMetFre)
    
    LenMot=2*w+2
    
    ScoDisPerMot=foreach(i=1:length(WeiLogPOM[[w]]), .combine=c) %dopar% {
      
      # #Log
      # line <- paste("Calling Scan for i=", i, " NumSeq=", NumSeq, " LenMot=", LenMot, sep="")
      # write(line,file=Config$LogFile,append=TRUE)
      
      # #Log
      # line <- "Call Scan"
      # write(line,file=Config$LogFile,append=TRUE)
    
      ScoDis=ScanSeqs(WeiLogPOM[[w]][i],WeiLogPMV[[w]][i],NumSeq,LenMot)
      
      # Extract binding counts.
      bincounts=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$counts
      # Extract binding breaks.
      binbreaks=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$breaks
      binbreaks=binbreaks[1:(length(binbreaks)-1)]
      ScoDatFram=data.frame(bincounts,binbreaks)
      ScoDatFram

    }
    ScoPerLen[[w]]=ScoDisPerMot
    }
  }
  
  # Structure adptation.
  Res <- list()
  
  for ( i in 1:length(ScoPerLen)){
    
    ParRes <- list()
    AuxLst <- ScoPerLen[[i]]
    
    j <- 1
    z <- 1
    
    while ( j <= length(AuxLst)){
      
      bincounts <- AuxLst[j]$bincounts
      binbreaks <- AuxLst[j+1]$binbreaks
      ParRes[[z]] <- data.frame(bincounts, binbreaks)
      
      j = j + 2
      z = z + 1
    }
    
    if ( length(ParRes) == 0 ) Res[[i]] <- NULL
    else Res[[i]] <- ParRes
    
  }
  
  return (Res)
}

ScanParFrag=function(Config,SeqMetFreW,WeiLogPMV,WeiLogPOM,LenMotif){
  
  # 
  # Parallel Scanning using only R code.
  # For reducing the size of the objects returned by fork childs,
  # parallelization is done over seuqential chunks of the total
  # number of POMs per motif length.
  # 
  
  ScoPerLen=list()
  registerDoMC(Config$nCPU)
  
  # For each length.
  for (w in LenMotif) {
    

    ScoDisPerMot=list()
    if (!is.null(SeqMetFreW[[w]]) && !is.null(WeiLogPOM[[w]]) && !is.null(WeiLogPMV[[w]])){
      print(paste("Scanning length", w, "number of POMs", length(WeiLogPOM[[w]])))
      # #Log
      # line <- paste("Scan for LenMotif ", w, sep="")
      # write(line,file=Config$LogFile,append=TRUE)
      # 
      # #Log
      # line <- paste("Scan for LenMotif ", w, sep="")
      # write(line,file=Config$LogFile,append=TRUE)
      
      # Every sequence for this length.
      SeqMetFre=SeqMetFreW[[w]]
      # Sequences passed to numerical format.
      NumSeq=IndBasSeq(SeqMetFre)
      # Get frequencies.
      FreqVec = SeqMetFre[[3]]
      
      LenMot=2*w+2
      
      # The number of chunks will be the number of digits
      # of the number of POMs.
      CnkNum <- nchar(trunc(length(WeiLogPOM[[w]])))
      
      CnkLst<-intervals(CnkNum, length(WeiLogPOM[[w]]))

      
      for (j in 1:CnkNum){
        
        Str<-unlist(CnkLst[j,1])
        End<-unlist(CnkLst[j,2])
        
        print(paste("Fragment",Str,"-",End))
        
        ScoDisPerMotPar=foreach(i=Str:End, .combine=c) %dopar% {
          
          # #Log
          # line <- paste("Calling Scan for i=", i, " NumSeq=", NumSeq, " LenMot=", LenMot, sep="")
          # write(line,file=Config$LogFile,append=TRUE)
          
          # #Log
          # line <- "Call Scan"
          # write(line,file=Config$LogFile,append=TRUE)
          
          ScoDis=ScanSeqsFreq(WeiLogPOM[[w]][i],WeiLogPMV[[w]][i],NumSeq,LenMot, FreqVec)
          #ScoDis=ScanSeqs(WeiLogPOM[[w]][i],WeiLogPMV[[w]][i],NumSeq,LenMot)
          
          # Extract binding counts.
          bincounts=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$counts
          # Extract binding breaks.
          binbreaks=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$breaks
          binbreaks=binbreaks[1:(length(binbreaks)-1)]
          ScoDatFram=data.frame(bincounts,binbreaks)
          ScoDatFram
          
        }
        ScoDisPerMot=c(ScoDisPerMot,ScoDisPerMotPar)
        
      }
      ScoPerLen[[w]]=ScoDisPerMot
    }
  }
  
  # Structure adptation.
  Res <- list()
  
  for ( i in 1:length(ScoPerLen)){
    ParRes <- list()
    if (!is.null(ScoPerLen[[i]]) && length(ScoPerLen[[i]])>0){
      AuxLst <- ScoPerLen[[i]]
      
      j <- 1
      z <- 1
      
      while ( j <= length(AuxLst)){
        
        bincounts <- AuxLst[j]$bincounts
        binbreaks <- AuxLst[j+1]$binbreaks
        ParRes[[z]] <- data.frame(bincounts, binbreaks)
        
        j = j + 2
        z = z + 1
      }
    }
    
    if ( length(ParRes) == 0 ) Res[[i]] <- NULL
    else Res[[i]] <- ParRes
    
  }
  
  return (Res)
}

POMsToVector = function(Config, WeiLogPOM, LenMotif, width) {
  POMMat = matrix (nrow = LenMotif, ncol = width * 4)
  for (i in 1:LenMotif){
    POMMat[i,] = as.vector(WeiLogPOM[i][[1]])
  }
  return (as.vector(t(POMMat)))
}

PMVsToVector = function(Config, WeiLogPMV, LenMotif, width) {
  PMVMat = matrix (nrow = LenMotif, ncol = width)
  for (i in 1:LenMotif){
    PMVMat[i,] = as.vector(WeiLogPMV[i][[1]])
  }
  return(as.vector(t(PMVMat)))
}

IndBasSeqVector=function(seqs){
  
  #
  # IndBasSeq
  # Auxiliary function. Passes a frame of nucleotidic sequences into
  # numeric format.
  # a -> 0
  # c -> 1
  # g -> 2
  # t -> 3
  
  
  Split=strsplit(seqs,NULL,fixed=TRUE)
  SeqSep=t(sapply(1:length(Split), function(x) unlist(Split[[x]])))
  
  SeqSep[SeqSep == "a"] = 0
  SeqSep[SeqSep == "c"] = 1
  SeqSep[SeqSep == "g"] = 2
  SeqSep[SeqSep == "t"] = 3
  
  class(SeqSep)='numeric'
  
  #SeqSep=lapply(seq_len(nrow(SeqSep)), function(i) SeqSep[i,])
  
  return(as.vector(t(SeqSep))) 
}

ScanFast = function(Config,Seqs,WeiLogPMV,WeiLogPOM,LenMotif) {
  ScoPerLen=list()
  
  for (w in LenMotif) {
    
    if (!is.null(Seqs[[w]]) && !is.null(WeiLogPOM[[w]]) && !is.null(WeiLogPMV[[w]])){
      
      print(paste("Scanning length", w, "number of POMs", length(WeiLogPOM[[w]])))
      
      # Every sequence for this length.
      Seq = Seqs[[w]]
      # Sequences passed to numerical format.
      Seq = IndBasSeqVector(Seq)
      nPOMs = length(WeiLogPOM[[w]])
      nSeqs = length(Seqs[[w]])
      
      LenMot = 2 * w + 2
      
      POMvec <- POMsToVector(Config, WeiLogPOM[[w]], nPOMs, 2 * w + 2)
      PMVvec <- PMVsToVector(Config, WeiLogPMV[[w]], nPOMs, 2 * w + 2)
      
      res = .Call( "scanPOMs", POMvec, PMVvec, nPOMs, Seq, nSeqs, 2 * w + 2, Config$nBins, 1 )
      
      ScoDisPerMot = list()
      ScoDisPerMot <- lapply(seq(0,(nPOMs-1)), 
                             function(i) data.frame(res[((i*Config$nBins)+1):((i+1)*Config$nBins)], 
                                                    res[((nPOMs*Config$nBins)+(i*Config$nBins+1)):((nPOMs*Config$nBins)+((i+1)*Config$nBins))]))

      for ( i in 1:nPOMs )
        colnames(ScoDisPerMot[[i]]) <- c("bincounts","binbreaks")
      
      ScoPerLen[[w]] = ScoDisPerMot
    }
  }
  return (ScoPerLen)
}


ScanParFragTot=function(Config,Seqs,WeiLogPMV,WeiLogPOM,LenMotif){
  
  # 
  # Parallel Scanning using only R code.
  # For reducing the size of the objects returned by fork childs,
  # parallelization is done over seuqential chunks of the total
  # number of POMs per motif length.
  # 
  
  ScoPerLen=list()
  registerDoMC(Config$nCPU)
  
  # For each length.
  for (w in LenMotif) {
    
    if (!is.null(Seqs[[w]]) && !is.null(WeiLogPOM[[w]]) && !is.null(WeiLogPMV[[w]])){
    print(paste("Scanning length", w, "number of POMs", length(WeiLogPOM[[w]])))
    # #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    # 
    # #Log
    # line <- paste("Scan for LenMotif ", w, sep="")
    # write(line,file=Config$LogFile,append=TRUE)
    
    # Every sequence for this length.
    Seq=Seqs[[w]]
    # Sequences passed to numerical format.
    NumSeq=IndBasSeqTot(Seq)
    
    LenMot=2*w+2
    
    # The number of chunks will be the number of digits
    # of the number of POMs.
    CnkNum <- nchar(trunc(length(WeiLogPOM[[w]])))
    
    CnkLst<-intervals(CnkNum, length(WeiLogPOM[[w]]))
    
    ScoDisPerMot=list()
    
    for (i in 1:CnkNum){
      
      Str<-unlist(CnkLst[i,1])
      End<-unlist(CnkLst[i,2])
      
      print(paste("Fragment",Str,"-",End))
      
      ScoDisPerMotPar=foreach(i=Str:End, .combine=c) %dopar% {
        
        # #Log
        # line <- paste("Calling Scan for i=", i, " NumSeq=", NumSeq, " LenMot=", LenMot, sep="")
        # write(line,file=Config$LogFile,append=TRUE)
        
        # #Log
        # line <- "Call Scan"
        # write(line,file=Config$LogFile,append=TRUE)
        
        ScoDis=ScanSeqs(WeiLogPOM[[w]][i],WeiLogPMV[[w]][i],NumSeq,LenMot)
        
        # Extract binding counts.
        bincounts=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$counts
        # Extract binding breaks.
        binbreaks=hist(ScoDis, breaks =Config$nBins, plot = FALSE)$breaks
        binbreaks=binbreaks[1:(length(binbreaks)-1)]
        ScoDatFram=data.frame(bincounts,binbreaks)
        ScoDatFram
        
      }
      ScoDisPerMot=c(ScoDisPerMot,ScoDisPerMotPar)
      
    }
    ScoPerLen[[w]]=ScoDisPerMot
    }
  }
  
  # Structure adptation.
  Res <- list()
  
  for ( i in 1:length(ScoPerLen)){
    
    ParRes <- list()
    if (!is.null(ScoPerLen[[i]]) && length(ScoPerLen[[i]])>0){
      AuxLst <- ScoPerLen[[i]]
      
      j <- 1
      z <- 1
      
      while ( j <= length(AuxLst)){
        
        bincounts <- AuxLst[j]$bincounts
        binbreaks <- AuxLst[j+1]$binbreaks
        ParRes[[z]] <- data.frame(bincounts, binbreaks)
        
        j = j + 2
        z = z + 1
      }
    }
    
    if ( length(ParRes) == 0 ) Res[[i]] <- NULL
    else Res[[i]] <- ParRes
    
  }
  
  return (Res)
}

GetLen2=function(Config,MotfAftTest,LenMotif){
  
  Len=rep(0,length(LenMotif))  
  
  for(w in LenMotif){
    
    if (!is.null(MotfAftTest[[w]])){
    if(length(MotfAftTest[[w]][[1]]) != 0 ) Len[w]=w
    }
    else {Len[w]=0}
  }
  Len=which(Len!=0)
  return (Len)
}

ExtractFDR=function(Config,MotAftTes,LenMotif){
  
  # Extract FDR column from the structure list.
  
  #Log
  #line <- paste("Starting ExtractFDR ", MotAftTes, " ", LenMotif)
  #write(line,file=Config$LogFile,append=TRUE)
  
  FDR=list()
  i=0
  for(w in LenMotif){
    if (!is.null(MotAftTes[[w]]) && (length(MotAftTes[[w]])>=5)){
    i=i+1
    FDR[[i]]=MotAftTes[[w]][[5]]
    }
  }
  
  #Log
  line <- "Ending ExtractFDR "
  write(line,file=Config$LogFile,append=TRUE)
  
  return(FDR)
}

GetLen3=function(Config,MotifSelec,LenMotif){
  
  Len=rep(0,length(LenMotif))  
  
  for(w in LenMotif){
    if(!is.null(MotifSelec[[w]])){
    if(length(MotifSelec[[w]][[1]]) != 0 ) Len[w]=w
    }
    else { Len[w] = w }
  }
  Len=which(Len!=0)
  return (Len)
}

GetCooChr=function(Config,FinMot,LenMotif){
  
  CooChr=list()
  #temporaldir=tempdir()
  temporaldir=Config$PathOutput
  CooFilDir=paste(temporaldir,"/CoordinateFiles",sep="")
  
  for(w in LenMotif){
    
    if (!is.null(FinMot[[w]]) && (length(FinMot[[w]])>=2)){
    print(paste("length",w))
    #Log
    # line <- paste("Starting GetCooChr LenMotif ", w)
    # write(line,file=Config$LogFile,append=TRUE)
    
    # All indexes from each selected motif.
    Indlist=FinMot[[w]][[2]]
    # Number of indexes for each selected motif.
    NumIndPerSeq=sapply(1:length(Indlist),function(i) length(Indlist[[i]]))
    CumSumIndPerSeq=cumsum(NumIndPerSeq)
    VecRow=unlist(Indlist)
    IndOrder=order(VecRow)
    SortVecRow=sort(VecRow)
    LenVecRow=length(VecRow)
    
    #Log
    # line <- paste("Calling readCooChrFile LenMotif", w)
    # write(line,file=Config$LogFile,append=TRUE)
    
    #Log
    # line <- paste("Calling readCooChrFile LenMotif ", w, sep="")
    # write(line,file="readCooChrFile.logs",append=TRUE)
    

    CooChrLis=.Call("readCooChrFile",SortVecRow,LenVecRow,w,CooFilDir)
    CooChrLis=CooChrLis[IndOrder]
    
    bindlist=function(i){
      
      if(i==1) cbind(CooChrLis[1:CumSumIndPerSeq[i]])
      else cbind(CooChrLis[(CumSumIndPerSeq[i-1]+1):CumSumIndPerSeq[i]])   
    }
    CooChrListPerSeq=lapply(1:length(CumSumIndPerSeq),bindlist)
    CooChr[[w]]=CooChrListPerSeq
    }
  }
  return (CooChr)
}

MatchMotifToGene=function(Config,seq_gene,CooChr,LenMotif){
  
  GenePerChr=list()
  
  for(i in 1:Config$NumChr){
    
    ind <- which(seq_gene$Chr == i)
    SeqChrFirVec <- seq_gene[ind,]
    SeqChrSecVec <- SeqChrFirVec
    SeqChrFirVec$Chr <- NULL
    SeqChrFirVec$PosEnd <- NULL
    SeqChrSecVec$Chr <- NULL
    SeqChrSecVec$PosIni <-NULL
    SeqChrFirVec <- setNames(SeqChrFirVec,rep(" ", 3))
    SeqChrSecVec <- setNames(SeqChrSecVec,rep(" ", 3))
    SeqChr <- rbind(SeqChrFirVec,SeqChrSecVec)
    SeqChr[,2] <- as.character(SeqChr[,2])
    GenePerChr[[i]] <- SeqChr
  } 
  
  MotPerLenGeneID=list()
  
  for(w in LenMotif){
    GeneID=list()
    MotGeneID=list()
    if (!is.null(CooChr[[w]]) && (length(CooChr[[w]])>0)){
      for(i in 1:length(CooChr[[w]])){
            VecChrCoo=unlist(CooChr[[w]][[i]][,1])
            MatChrCoo=t(matrix(VecChrCoo,nrow=2,ncol=length(VecChrCoo)/2))
            GeneVec=list()
            for(i2 in 1:nrow(MatChrCoo)){
                In=MatChrCoo[i2,1]
                if (!is.na(In) & In > 0){
                  Distance=abs(GenePerChr[[In]][,1]-MatChrCoo[i2,2])
                  IndMin=which(Distance==min(Distance))
                  GeneVec[[i2]]=c(GenePerChr[[In]][IndMin,2],GenePerChr[[In]][IndMin,3])
                }
                else{
                  #Log
                  # line <- paste(" In gene matching, length", w, "chr", i, "row", i2, "failed") 
                  # write(line,file=Config$LogFile,append=TRUE)
                }
            }
            MotGeneID[[i]]=GeneVec
      }
    }
    MotPerLenGeneID[[w]]=MotGeneID
  }
  return (MotPerLenGeneID)
}
 
LogoScanPlot=function(Config,MotifSelec,LenMotif,TypeMotif){

  Mode=switch(TypeMotif, 
           ForProne = "_p_PM_", 
           RevProne = "_n_PM_",
           ForResis = "_p_RM_",
           RevResis = "_n_RM_")

  DifResult=list()
  	
  for(w in LenMotif){
    
    if (!is.null(MotifSelec[[w]]) && (length(MotifSelec[[w]])>=4)) {
      DifResult[[w]]=list()
      j = 1
      PathLogosDir=Config$Logos
      PathScanDir=Config$Scan
    
      DistProne=MotifSelec[[w]][[3]]
      DistResis=MotifSelec[[w]][[4]]
    
      Pwm=MotifSelec[[w]][[1]]
      names = c(1:length(DistProne))
    
      for(i in 1:length(DistProne)){
    
        ###Scan Plots
        
        ScoDatFraPro=DistProne[[i]]
        ScoDatFraRes=DistResis[[i]]
        p=unlist(sapply(1:length(ScoDatFraPro[,1]), function(i) rep(ScoDatFraPro[,2][i],ScoDatFraPro[,1][i])))
        r=unlist(sapply(1:length(ScoDatFraRes[,1]), function(i) rep(ScoDatFraRes[,2][i],ScoDatFraRes[,1][i])))
      
        dat1 <- data.frame(Set = factor(rep(c("Prone"), each=length(p))), Score = p)
        
        dat2 <- data.frame(Set = factor(rep(c("Resistant"), each=length(r))), Score = r)
      
        dat=rbind(dat1,dat2)
        cdat <- ddply(dat, "Set", summarise, rating.mean=mean(Score))
        DifResult[[w]][[j]] = cdat
        j = j + 1 

	if(Config$DrawLogo==TRUE){
     

        mypathScan <- file.path(PathScanDir, paste(Config$MotifTarget,"_", Config$DatasetName , Mode, Config$X,"_", 2*w+2, "_", names[i],".png",sep=""))

        # TODO: ajustar.
        # Fondo transparente.
        # Leyenda red/green por prone/resistant
        # FDR de resistant a otro color
        # Colores más vivos.
        # Info extra de secuencias en el html (nº secuencias por logo)
        ggplot(dat, aes(x=Score, fill=Set,)) +
          geom_histogram(aes(y=..count..),binwidth=.5, alpha=.3, position="identity") +
          xlab("Matching Score") +
          ylab("Frequency") +
          theme(axis.title.x = element_text(colour = "gray"),
                axis.title.y = element_text(colour = "gray")) +
          theme(panel.background = element_rect(fill='white', colour='gray')) +
          scale_fill_manual("",values = c('green', 'red')) +
          geom_vline(data=cdat[1,] ,aes(xintercept=rating.mean, colour="red"),
                     linetype="longdash", size=1.4) +
          geom_vline(data=cdat[2,] ,aes(xintercept=rating.mean, colour="green"),
                     linetype="longdash", size=1.4)
        
        ggsave(filename = mypathScan,width=5,height=5)
        
        ###Logo Plots
        mypathLogo <- file.path(PathLogosDir, paste(Config$MotifTarget,"_", Config$DatasetName, Mode, Config$X,"_", 2*w+2,"_", names[i],".png",sep=""))
        
        png(filename = mypathLogo)
        pwm=makePWM(Pwm[[i]])
        seqLogo(pwm, ic.scale = FALSE)
        dev.off()
	}
      }
  }
  }
  return(DifResult)
}

ReportHtml=function(Config,MotifSelec,GeneID,LenMotif,TypeMotif){
  
  PathHtml=switch(TypeMotif, 
              ForProne = Config$ForProneHtml, 
              RevProne = Config$RevProneHtml,
              ForResis = Config$ForResisHtml,
              RevResis = Config$RevResisHtml)
  
  ModePlot=switch(TypeMotif, 
              ForProne = "_p_PM_", 
              RevProne = "_n_PM_",
              ForResis = "_p_RM_",
              RevResis = "_n_RM_")
  
  InitHtml=c("<TABLE BORDER=7 CELLPADDING=10>"
             ,"<TR>"
             ,"<TH>Motif ID</TH>"
             ,"<TH>Motif Logo</TH>"
             ,"<TH>Scanning of Binding Energy Plot</TH>"
             ,"<TH>Gene Targets</TH>"
             ,"</TR>")
  
  
  FinHtml=c("</TABLE>"
            ,"</body>"
            ,"</html>")
  
  pF=sink(PathHtml)
  cat(InitHtml,pF,sep = "\n")
  
  for(w in LenMotif){
  
    if(!is.null(MotifSelec[[w]]) && (length(MotifSelec[[w]])>0)){
      Len=length(MotifSelec[[w]][[1]])
      names=c(1:Len)
    
      for(i in 1:Len){

        cat("<TR>\n")
        cat("<TH><font size=1>")
        # Name of the motif.
        cat(paste(Config$MotifTarget,"_",Config$DatasetName,ModePlot,Config$X,"_",2*w+2,"_",names[i],sep=""))
        cat("</TH>\n")
        cat("<TH><IMG SRC=.././Figures/Logos/")
        cat(paste(Config$MotifTarget,"_",Config$DatasetName,ModePlot,Config$X,"_",2*w+2,"_",names[i],".png",sep=""))
        cat(" width=400 ALT=logo></TH>\n")
        cat("<TH><IMG SRC=.././Figures/ScanPlot/")
        cat(paste(Config$MotifTarget,"_",Config$DatasetName,ModePlot,Config$X,"_",2*w+2,"_",names[i],".png",sep=""))
        cat(" width=400 ALT=scan></TH>\n")
        cat("<TH><font size=1>")
        
          
        for(i2 in 1:length(GeneID[[w]][[i]])){
            
            GeID=unlist(GeneID[[w]][[i]][i2])
            if (length(GeID)!=0){
              NumAs <- length(GeID)/2
              for (j in 1:NumAs){
                PosNam <- j
                PosCod <- j + NumAs
                if (!is.na(GeID[PosNam]) && !is.na(GeID[PosCod])){
                  cat(paste("<a href=http://www.ncbi.nlm.nih.gov/gene/",as.numeric(GeID[PosCod]),"><i>",GeID[PosNam],sep=""))
                  cat("</i></a>\n") 
                }
                else{
                  #Log
                  # line <- paste(" Motif", paste(Config$MotifTarget,"_",Config$InputDataFileName,ModePlot,Config$X,"_",2*w+2,"_",names[i],sep=""), " no genes.") 
                  # write(line,file=Config$LogFile,append=TRUE)
                }
              }
              
            }
        }
        cat("</TH>\n")
        cat("</TR>\n")
      }
    }
  }
  cat(FinHtml,pF,sep="\n")
  sink()
}

IndexHtmlReport=function(Config){

  pF=sink(Config$IndexHtml)
  
  IniIndexHtml=c("<html><head><title>Index DNA Methylation motifs</title>"
                     ,"<!-- DOCNAME: DNA methylation motifs -->"
                     ,"<!-- CHUNKNAME: Index -->"
                     ,"<!-- CHAPNAME: Index -->"
                     ,"<!-- HEADSTUFF -->"
                     ,"</head>"
                     ,"<body bgcolor=#ffffff>"
                     ,"<a name=998339>"
                     ,"<!-- NAVBARTOP -->"
                     ,"<table border=0 width=100% cellpadding=0 cellspacing=0><tr>"
                     ,"<td valign=baseline bgcolor=ffe4b0><b>DNA methylation motifs</b></td>"
                     ,"<td valign=baseline bgcolor=ffe4b0 align=right>"
                     ,"<a href=ReaFil.htm><img src=b_next.gif border=0></a></td></tr></table>"
                     ,"<br>"
                     ,"<p>Gene expression regulation is gated by the promoter methylation states 
enabling or disabling transcription factors to exert their fine tuning regulation. 
The known DNA methylation/unmethylation mechanisms are sequence unspecific but different cell types with same genome have different methylomes. 
Thus, additional cell-type specific processes bringing specificity to the methylation/unmethylation mechanisms are needed. 
Searching for such processes, 
we investigate whether a relation between CpGs methylation and CpGs sequence context exists. 
We focus on the cellular reprogramming problem because its mechanism is poorly understood, 
but it is hypothesized to be based on a crosstalk between genetic and epigenomic networks.
We demonstrated that the CpGs methylation states are influenced by the sequence context surrounding the CpGs. 
We used such property to develop a CpG methylation motif discovery algorithm. 
The motifs found discriminate between hypomethylated and hypermethylated regions.  
The algorithm can be used in any NGS based methylomics study."
                     ,"</p>"
                     ,"<p>We term the DNA methylation motifs centered in each CpG as CpGMM.
This page helps to navigate through the differerent DNA methylation 
motifs found by our algorithm in the context of human pluripotent cells.</p>"
                     ,"<!-- H2 --><font size=+2 color=990000><b>All found CpGMMs</b></font><br></p>"
                     ,"<p>"
                     ,"Here we show the logos all the found motifs.
We show the distribution of their binding affinity to resistant- and prone-methytlation regions."
                     ,"</p>")
  
  FinIndexHtml=c("<br>"
            ,"<!-- Copyright 2012 Marcos J. Arauzo-Bravo. -->"
            ,"<!-- Last updated: Fri Oct 5 22:05:19 2012 -->"
            ,"</body></html>")
  
  cat(IniIndexHtml,pF,sep = "\n")
  
  cat("<p><font size=+1 color=990000><b>CpGMMs</b></font><br>\n")
  cat("<ul>\n")
  cat("<li><a href=./DocHtml/")
  cat(paste(Config$MotifTarget,"_",Config$DatasetName,"_disp",Config$X,"_p_prone.html",sep=""))
  cat(">Forward prone-methylation CpGMMs</a>.\n")
  cat("<li><a href=./DocHtml/")
  cat(paste(Config$MotifTarget,"_",Config$DatasetName,"_disp",Config$X,"_p_resistant.html",sep=""))
  cat(">Forward resistant-methylation CpGMMs</a>.\n")
  cat("<li><a href=./DocHtml/")
  cat(paste(Config$MotifTarget,"_",Config$DatasetName,"_disp",Config$X,"_n_prone.html",sep=""))
  cat(">Reverse prone-methylation CpGMMs</a>.\n")
  cat("<li><a href=./DocHtml/")
  cat(paste(Config$MotifTarget,"_",Config$DatasetName,"_disp",Config$X,"_n_resistant.html",sep=""))
  cat(">Reverse resistant-methylation CpGMMs</a>.\n")
  cat("<ul>\n")
  
  cat(FinIndexHtml,pF,sep="\n")
    
  sink()  
}


## Main function
ListFDR =DNAMethylationMotifFinding=function(Config){
  
  print(paste("lambda",Config$lambda))
  #Filter .md file
  
  t1 <- Sys.time()
  print("Filtering and reading .md file")
  Seq_Gene = filtreadmdfile(Config)
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "Filtering done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #Read targets' Coordinates
  print(paste("Reading ",Config$MotifTarget," coordinates", sep=""))
  t1 <- Sys.time()
  CooMet = ReadCooMet(Config)
  t2 <- Sys.time()
  print(t2-t1)
  
  CooMetFor = CooMet[[1]]
  CooMetRev = CooMet[[2]]
  
  #Log
  # line <- paste(Config$MotifTarget," Coord read done",sep="")
  # write(line,file=Config$LogFile,append=TRUE)
  
  # Optionally generate methylation histogram.
  #TODO: Methylation histogram from CooMetFor and CooMetRev
  if (Config$MetHist == TRUE) {
    GenMetHst(Config, CooMetFor)
  }

  
  #####Displace one unit the coordinate from forward strand
  # CooMetForUpd = CooMov(Config,CooMetFor, "for")
  CooMetForUpd = CooMov(Config,CooMetFor)
  # CooMetRevUpd = CooMov(Config,CooMetRev, "rev")
  rm(CooMetFor, CooMet)
  gc()
  #Log
  # line <- "Displace done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####Read Fasta 
  t1 <- Sys.time()
  print("Reading fasta")
  SeqChr = ReadFasta(Config)
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "Fasta read done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #save(Seq_Gene,
  #     CooMetForUpd,
  #     SeqChr,
  #     Config,
  #     file=paste("beforeDicWord28-08-2019Disp",Config$X,".RData",sep=""))
  
  #####Compilation of CpG Word Dictionaries
  t1 <- Sys.time()
  print("Generating word dictionary - forward")
  SeqMetFreFor = DicWord(Config,CooMetForUpd,SeqChr, "for")
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "Forward compilation done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  t1 <- Sys.time()
  print("Generating word dictionary - reverse")
  SeqMetFreRev = DicWord(Config,CooMetRev,SeqChr, "rev")
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "Reverse compilation done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #print(paste("SeqMetFreFor size:",format(object.size(SeqMetFreFor), units="Gb")))
  #save(SeqMetFreFor,
  #     file=paste("SeqMetFreFor28-08-2019Disp",Config$X,".RData",sep=""))
  
  # Get all sequences for scanning.
  SeqTotFor = SeqMetFreFor$SeqTot
  SeqTotRev = SeqMetFreRev$SeqTot
  # Get sequences for POM generation.
  SeqMetFreFor = SeqMetFreFor$SeqMetFreW
  SeqMetFreRev = SeqMetFreRev$SeqMetFreW
  
  # Reverse of total sequences.
  SeqTot = list()
  SeqTot$SeqPrnFor = DelGapsTot(Config, SeqTotFor$SeqPrn)
  SeqTot$SeqResFor = DelGapsTot(Config, SeqTotFor$SeqRes)
  SeqTot$SeqPrnRev = DelGapsTot(Config, SeqTotRev$SeqPrn)
  SeqTot$SeqResRev = DelGapsTot(Config, SeqTotRev$SeqRes)
  
  save(SeqTot, file=paste("SeqTotSaving",Config$X,".RData",sep=""))
  
  # save(SeqMetFreFor,
  #      SeqMetFreRev,
  # Config,
  # file=paste("afterDicWord28-08-2019Disp",Config$X,".RData",sep=""))
  
  rm(SeqTot)
  gc()

  #####Extract DicWord,SeqFre,calculate nLines 
  #print("Reverse")
  #SeqMetFreRev = Rev(Config,SeqMetFreFor)
  rm(SeqChr)
  gc()
  #Log
  # line <- "Extract, calculate done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####Delete Gaps from CpG Word Dictionaries
  print("Deleting Gaps")
  t1 <- Sys.time()
  SeqMetFreFor=DelGaps(Config,SeqMetFreFor)
  SeqMetFreRev=DelGaps(Config,SeqMetFreRev)
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "Gap delection done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####Classification into two subsets based on their methylation ratio
  print("Classification")
  t1 <- Sys.time()
  SeqMetFreForProne = ProneMet(Config,SeqMetFreFor)
  SeqMetFreForResis = ResisMet(Config,SeqMetFreFor)
  SeqMetFreRevProne = ProneMet(Config,SeqMetFreRev)
  SeqMetFreRevResis = ResisMet(Config,SeqMetFreRev)
  t2 <- Sys.time()
  print(t2-t1)  
  #Log
  # line <- "Classification done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  ##Fusion of CpG word dictionaries
  print("Fusion")
  t1 <- Sys.time()
  print("Forward prone")
  FreWVecForProne = FreqVec(Config,SeqMetFreForProne)
  SeqMetFreWFreVecWForProne = Fusion(Config,SeqMetFreForProne,FreWVecForProne)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Forward resistant")
  FreWVecForResis = FreqVec(Config,SeqMetFreForResis)
  SeqMetFreWFreVecWForResis = Fusion(Config,SeqMetFreForResis,FreWVecForResis)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse prone")
  FreWVecRevProne = FreqVec(Config,SeqMetFreRevProne)
  SeqMetFreWFreVecWRevProne = Fusion(Config,SeqMetFreRevProne,FreWVecRevProne)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse resistant")
  FreWVecRevResis = FreqVec(Config,SeqMetFreRevResis)
  SeqMetFreWFreVecWRevResis = Fusion(Config,SeqMetFreRevResis,FreWVecRevResis)
  t2 <- Sys.time()
  print(t2-t1)
  

#   save (SeqMetFreForProne,
# 	SeqMetFreForResis,
# 	SeqMetFreRevProne,
# 	SeqMetFreRevResis,
# 	file = paste("beforeFusion08-05-2020Disp",Config$X,".RData",sep="") )	
  
  # save (SeqMetFreWFreVecWForProne,
  #       SeqMetFreWFreVecWForResis,
  #       SeqMetFreWFreVecWRevProne,
  #       SeqMetFreWFreVecWRevResis,
  #       file = paste("afterFusion08-05-2020Disp",Config$X,".RData",sep="") )

  #Log
  # line <- "Fusion done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####POM, POMVec
  print("POM")
  
  # Position Ocurrence matrix.
  
  # save(Config,
  #      SeqMetFreWFreVecWForProne,
  #      SeqMetFreWFreVecWForResis,
  #      SeqMetFreWFreVecWRevProne,
  #      SeqMetFreWFreVecWRevResis,
  #      file=paste("beforePOM28-08-2019Disp",Config$X,".RData",sep=""))
  
  print("Generating POMs")
  t1 <- Sys.time()
  SeqMetFreWFreVecWForProne = ReduceWords(Config, SeqMetFreWFreVecWForProne)
  SeqMetFreWFreVecWForResis = ReduceWords(Config, SeqMetFreWFreVecWForResis)
  SeqMetFreWFreVecWRevProne = ReduceWords(Config, SeqMetFreWFreVecWRevProne)
  SeqMetFreWFreVecWRevResis = ReduceWords(Config, SeqMetFreWFreVecWRevResis)

  POMForProne = POM(Config,SeqMetFreWFreVecWForProne)
  POMForResis = POM(Config,SeqMetFreWFreVecWForResis)
  POMRevProne = POM(Config,SeqMetFreWFreVecWRevProne)
  POMRevResis = POM(Config,SeqMetFreWFreVecWRevResis)
  
  POMVecForProne = POMVec(Config,POMForProne)
  POMVecForResis = POMVec(Config,POMForResis)
  POMVecRevProne = POMVec(Config,POMRevProne)
  POMVecRevResis = POMVec(Config,POMRevResis)
  t2 <- Sys.time()
  print(t2-t1)
    
  # #Log
  # line <- "POM done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####Hierarchical Clustering
  print("Clustering")
  t1 <- Sys.time()
  print("Forward prone")
  POMCluForProne = PomClu(Config,POMVecForProne,SeqMetFreWFreVecWForProne)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Forward resistant")
  POMCluForResis = PomClu(Config,POMVecForResis,SeqMetFreWFreVecWForResis)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse prone")
  POMCluRevProne = PomClu(Config,POMVecRevProne,SeqMetFreWFreVecWRevProne)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse resistant")
  POMCluRevResis = PomClu(Config,POMVecRevResis,SeqMetFreWFreVecWRevResis)
  t2 <- Sys.time()
  print(t2-t1)
  
  # save(Config,
  #      POMCluForProne,
  #      POMCluForResis,
  #      POMCluRevProne,
  #      POMCluRevResis,
  #      file=paste("afterClust28-08-2019Disp",Config$X,".RData",sep=""))  

  #Log
  # line <- "HClustering done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  ####Get lengths of the existing motifs
    
  LenMotifForProne = GetLen(Config,POMCluForProne)
  LenMotifForResis = GetLen(Config,POMCluForResis)
  LenMotifRevProne = GetLen(Config,POMCluRevProne)
  LenMotifRevResis = GetLen(Config,POMCluRevResis)
  #Log
  # line <- "GetLengths done"
  # write(line,file=Config$LogFile,append=TRUE)
  
    
  #####PMV,PWM,WeiPOM,WeiLogPMV,WeiLogPOM
  print("Weighted POMs")
  t1 <- Sys.time()
  PMVForProne = PMV(Config,POMCluForProne,LenMotifForProne)
  PMVForResis = PMV(Config,POMCluForResis,LenMotifForResis)
  PMVRevProne = PMV(Config,POMCluRevProne,LenMotifRevProne)
  PMVRevResis = PMV(Config,POMCluRevResis,LenMotifRevResis)
  
  PWMForProne = PWM(Config,POMCluForProne,LenMotifForProne)
  PWMForResis = PWM(Config,POMCluForResis,LenMotifForResis)
  PWMRevProne = PWM(Config,POMCluRevProne,LenMotifRevProne)
  PWMRevResis = PWM(Config,POMCluRevResis,LenMotifRevResis)
  
  WeiPOMForProne = WeiPOM(Config,POMCluForProne,LenMotifForProne)
  WeiPOMForResis = WeiPOM(Config,POMCluForResis,LenMotifForResis)
  WeiPOMRevProne = WeiPOM(Config,POMCluRevProne,LenMotifRevProne)
  WeiPOMRevResis = WeiPOM(Config,POMCluRevResis,LenMotifRevResis)
  
  WeiLogPMVForProne = WeiLogPMV(Config,WeiPOMForProne,PMVForProne,LenMotifForProne)
  WeiLogPMVForResis = WeiLogPMV(Config,WeiPOMForResis,PMVForResis,LenMotifForResis)
  WeiLogPMVRevProne = WeiLogPMV(Config,WeiPOMRevProne,PMVRevProne,LenMotifRevProne)
  WeiLogPMVRevResis = WeiLogPMV(Config,WeiPOMRevResis,PMVRevResis,LenMotifRevResis)
  
  WeiLogPOMForProne = WeiLogPOM(Config,WeiPOMForProne,POMCluForProne,LenMotifForProne)
  WeiLogPOMForResis = WeiLogPOM(Config,WeiPOMForResis,POMCluForResis,LenMotifForResis)
  WeiLogPOMRevProne = WeiLogPOM(Config,WeiPOMRevProne,POMCluRevProne,LenMotifRevProne)
  WeiLogPOMRevResis = WeiLogPOM(Config,WeiPOMRevResis,POMCluRevResis,LenMotifRevResis)
  t2 <- Sys.time()
  print(t2-t1)
    
  #Log
  # line <- "Weighted POMs done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #####Scan
  print("Scanning")

  
  #   'ss': Scan POMs against the selected subsets of resistant and prone words. Appropriate when 
  #     bith subsets have similar number of sequences for each length.
  if (Config$ScnTpe=='ss') {
    print("ss Scanning.")
    ScoDisHigMetForProne = ScanParFrag(Config,SeqMetFreForProne,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #Log
    # line <- "First scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetForResis = ScanParFrag(Config,SeqMetFreForResis,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #Log
    # line <- "Second scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetRevProne = ScanParFrag(Config,SeqMetFreRevProne,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    #Log
    # line <- "Third scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetRevResis = ScanParFrag(Config,SeqMetFreRevResis,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    #Log
    # line <- "Fourth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetForProne = ScanParFrag(Config,SeqMetFreForProne,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #Log
    # line <- "Fifth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetForResis = ScanParFrag(Config,SeqMetFreForResis,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #Log
    # line <- "Sixth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetRevProne = ScanParFrag(Config,SeqMetFreRevProne,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #Log
    # line <- "Seventh scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetRevResis = ScanParFrag(Config,SeqMetFreRevResis,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #Log
    # line <- "Eighth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
  }
  
  #   'sr': Scan prone POMs against the prone subset and a bunch of random sequences from the resistant
  #     total sequences (same number of sequences as the prone subset). Vice versa for resistant POMs.
  if(Config$ScnTpe=='sr') {
    print("sr scanning.")
    
    load(file=paste("SeqTotSaving",Config$X,".RData",sep=""))
    
    #Select random sampling of the total bunch of sequences.
    #RanSeqForRes <- SelSamples(SeqTot$SeqResFor, SeqMetFreForProne, LenMotifForProne)
    #RanSeqForPrn <- SelSamples(SeqTot$SeqPrnFor, SeqMetFreForResis, LenMotifForResis)
    #RanSeqRevRes <- SelSamples(SeqTot$SeqResRev, SeqMetFreRevProne, LenMotifRevProne)
    #RanSeqRevPrn <- SelSamples(SeqTot$SeqPrnRev, SeqMetFreRevResis, LenMotifRevResis)
    
    print("Selecting samples")
    t1 <- Sys.time()
    RanSeqForRes <- SelSamples2(SeqTot$SeqResFor, SeqTot$SeqPrnFor, LenMotifForResis)
    RanSeqForPrn <- SelSamples2(SeqTot$SeqPrnFor, SeqTot$SeqResFor, LenMotifForProne)
    RanSeqRevRes <- SelSamples2(SeqTot$SeqResRev, SeqTot$SeqPrnRev, LenMotifRevResis)
    RanSeqRevPrn <- SelSamples2(SeqTot$SeqPrnRev, SeqTot$SeqResRev, LenMotifRevProne)
    t2 <- Sys.time()
    print(t2-t1)

    rm(SeqTot)
    gc()
    
    # ScoDisHigMetForProne = ScanParFrag(Config,SeqMetFreForProne,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #ScoDisHigMetForProne = ScanParFragTot(Config,RanSeqForPrn,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    t1 <- Sys.time()
    print("POM: Forward Prone - Seqs: Forward Prone")
    ScoDisHigMetForProne = ScanFast(Config,RanSeqForPrn,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Forward Resistant - Seqs: Forward Prone")
    #Log
    # line <- "First scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisHigMetForResis = ScanParFragTot(Config,RanSeqForRes,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #ScoDisHigMetForResis = ScanParFragTot(Config,RanSeqForPrn,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    ScoDisHigMetForResis = ScanFast(Config,RanSeqForPrn,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Reverse Prone - Seqs: Reverse Prone")
    #Log
    # line <- "Second scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisHigMetRevProne = ScanParFrag(Config,SeqMetFreRevProne,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    # ScoDisHigMetRevProne = ScanParFragTot(Config,RanSeqRevPrn,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    ScoDisHigMetRevProne = ScanFast(Config,RanSeqRevPrn,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Reverse Resistant - Seqs: Reverse Prone")
    #Log
    # line <- "Third scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisHigMetRevResis = ScanParFragTot(Config,RanSeqRevRes,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    #ScoDisHigMetRevResis = ScanParFragTot(Config,RanSeqRevPrn,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    ScoDisHigMetRevResis = ScanFast(Config,RanSeqRevPrn,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Forward Prone - Seqs: Forward Resistant")
    #Log
    # line <- "Fourth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisLowMetForProne = ScanParFragTot(Config,RanSeqForPrn,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #ScoDisLowMetForProne = ScanParFragTot(Config,RanSeqForRes,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    ScoDisLowMetForProne = ScanFast(Config,RanSeqForRes,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Forward Resistant - Seqs: Forward Resistant")
    #Log
    # line <- "Fifth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisLowMetForResis = ScanParFrag(Config,SeqMetFreForResis,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #ScoDisLowMetForResis = ScanParFragTot(Config,RanSeqForRes,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    ScoDisLowMetForResis = ScanFast(Config,RanSeqForRes,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Reverse Prone - Seqs: Reverse Resistant")
    #Log
    # line <- "Sixth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisLowMetRevProne = ScanParFragTot(Config,RanSeqRevPrn,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #ScoDisLowMetRevProne = ScanParFragTot(Config,RanSeqRevRes,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    ScoDisLowMetRevProne = ScanFast(Config,RanSeqRevRes,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
    print("POM: Reverse Resistant - Seqs: Reverse Resistant")
    #Log
    # line <- "Seventh scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    # ScoDisLowMetRevResis = ScanParFrag(Config,SeqMetFreRevResis,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    # ScoDisLowMetRevResis = ScanParFragTot(Config,RanSeqRevRes,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    ScoDisLowMetRevResis = ScanFast(Config,RanSeqRevRes,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #Log
    # line <- "Eighth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)

  }
  
  #   'tt': Scan prone and resistant POMs against the total set of prone and resistant sequences. 
  #     Great memory and computational cost.
  if(Config$ScnTpe=='tt') {
    load(file=paste("SeqTotSaving",Config$X,".RData",sep=""))
    ScoDisHigMetForProne = ScanParFragTot(Config,SeqTot$SeqPrnFor,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #Log
    # line <- "First scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetForResis = ScanParFragTot(Config,SeqTot$SeqResFor,WeiLogPMVForProne,WeiLogPOMForProne,LenMotifForProne)
    #Log
    # line <- "Second scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetRevProne = ScanParFragTot(Config,SeqTot$SeqPrnRev,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    #Log
    # line <- "Third scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisHigMetRevResis = ScanParFragTot(Config,SeqTot$SeqResRev,WeiLogPMVRevProne,WeiLogPOMRevProne,LenMotifRevProne)
    #Log
    # line <- "Fourth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetForProne = ScanParFragTot(Config,SeqTot$SeqPrnFor,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #Log
    # line <- "Fifth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetForResis = ScanParFragTot(Config,SeqTot$SeqResFor,WeiLogPMVForResis,WeiLogPOMForResis,LenMotifForResis)
    #Log
    # line <- "Sixth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetRevProne = ScanParFragTot(Config,SeqTot$SeqPrnRev,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #Log
    # line <- "Seventh scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    ScoDisLowMetRevResis = ScanParFragTot(Config,SeqTot$SeqResRev,WeiLogPMVRevResis,WeiLogPOMRevResis,LenMotifRevResis)
    #Log
    # line <- "Eighth scanning done"
    # write(line,file=Config$LogFile,append=TRUE)
    rm(SeqTot)
    gc()
  }
  
  file.remove(paste("SeqTotSaving",Config$X,".RData",sep=""))
  
  
  #Log
  # line <- "All scanning done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  # #Log
  # line <- "Scanning done"
  # write(line,file=Config$LogFile,append=TRUE)
  # 
  ####First Filter:FDR
  print("FDR")
  
  FP="ForProne"
  RP="RevProne"
  FR="ForResis"
  RR="RevResis"
  
  #context saving
  # save(Config, 
  #       ScoDisHigMetForProne,
  #       ScoDisHigMetRevProne,
  #       ScoDisHigMetForResis,
  #       ScoDisHigMetRevResis,
  #       ScoDisLowMetForProne,
  #       ScoDisLowMetForResis,
  #       ScoDisLowMetRevProne,
  #       ScoDisLowMetRevResis,
  #       PWMForProne,
  #       PWMRevProne,
  #       PWMForResis,
  #       PWMRevResis,
  #       POMCluForProne,
  #       POMCluRevProne,
  #       POMCluForResis,
  #       POMCluRevResis,
  #       LenMotifForProne,
  #       LenMotifRevProne,
  #       LenMotifForResis,
  #       LenMotifRevResis,
  #       FP,
  #       RP,
  #       FR,
  #       RR,
  #       Seq_Gene,
  #       file=paste("beforeFDR28-08-2019Disp",Config$X,".RData",sep=""))
  
  # save.image(file = "DMMD2019.beforeFDR17-07-2019-01.RData")
  #load(""DMMD2019.beforeFDR17-07-2019-01.RData")
  
  t1 <- Sys.time()
  print("Forward Prone POMs")
  MotAftTesForProne=FDR(Config,ScoDisHigMetForProne,ScoDisLowMetForProne,PWMForProne,POMCluForProne,LenMotifForProne,FP)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse Prone POMs")
  #Log
  # line <- "fdr 1 done "
  # write(line, file=Config$LogFile, append=TRUE)
  MotAftTesRevProne=FDR(Config,ScoDisHigMetRevProne,ScoDisLowMetRevProne,PWMRevProne,POMCluRevProne,LenMotifRevProne,RP)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Forward Resistant POMs")
  #Log
  # line <- "fdr 2 done "
  # write(line, file=Config$LogFile, append=TRUE)
  MotAftTesForResis=FDR(Config,ScoDisHigMetForResis,ScoDisLowMetForResis,PWMForResis,POMCluForResis,LenMotifForResis,FR)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse Resistant POMs")
  #Log
  # line <- "fdr 3 done "
  # write(line, file=Config$LogFile, append=TRUE)
  MotAftTesRevResis=FDR(Config,ScoDisHigMetRevResis,ScoDisLowMetRevResis,PWMRevResis,POMCluRevResis,LenMotifRevResis,RR)
  t2 <- Sys.time()
  print(t2-t1)
  #Log
  # line <- "fdr 4 done "
  # write(line, file=Config$LogFile, append=TRUE)
  
  # FDRForProne = ExtractFDR(Config,MotAftTesForProne,LenMotifForProne)
  # FDRRevProne = ExtractFDR(Config,MotAftTesRevProne,LenMotifRevProne)
  # FDRForResis = ExtractFDR(Config,MotAftTesForResis,LenMotifForResis)
  # FDRRevResis = ExtractFDR(Config,MotAftTesRevResis,LenMotifRevResis)
  # 
  # #Log
  # line <- "ExtractFDR done "
  # write(line, file=Config$LogFile, append=TRUE)
  # 
  # ListFDR = list(FDRForProne,FDRRevProne,FDRForResis,FDRRevResis)
  # 
  # VecFDR=unlist(ListFDR)
  # MeaFDR=mean(VecFDR)
  
  # TODO: retree was here.
  
  #Log
  # line <- "FDR done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  ####Get lengths of the existing motifs
    
  LenMotif2ForProne = GetLen2(Config,MotAftTesForProne,LenMotifForProne)
  LenMotif2ForResis = GetLen2(Config,MotAftTesForResis,LenMotifForResis)
  LenMotif2RevProne = GetLen2(Config,MotAftTesRevProne,LenMotifRevProne)
  LenMotif2RevResis = GetLen2(Config,MotAftTesRevResis,LenMotifRevResis)
  
  #Log
  # line <- "Lengths2 done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  #context saving
  # save(Config,
  #      LenMotif2ForProne,
  #      LenMotif2RevProne,
  #      LenMotif2ForResis,
  #      LenMotif2RevResis,
  #      FP,
  #      RP,
  #      FR,
  #      RR,
  #      Seq_Gene,
  #      MotAftTesForProne,
  #      MotAftTesRevProne,
  #      MotAftTesForResis,
  #      MotAftTesRevResis,
  #      file=paste("beforeKolmogorov28-08-2019Disp",Config$X,".RData",sep=""))
  #save.image(file = "DMMD2019.beforeKolmogorov17-07-2019-01.RData")
  #load("DMMD2019.beforeKolmogorov17-07-2019-01.RData)
  
  #####Second Filter: Kolmogorov-Test
  print("Kolmogorov")
  t1 <- Sys.time()
  print("Forward Prone POMs")
  MotifSelecForProne=CheckDifferencesInDistributions(Config,MotAftTesForProne,LenMotif2ForProne,FP)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse Prone POMs")
  MotifSelecRevProne=CheckDifferencesInDistributions(Config,MotAftTesRevProne,LenMotif2RevProne,RP)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Forward Resistant POMs")
  MotifSelecForResis=CheckDifferencesInDistributions(Config,MotAftTesForResis,LenMotif2ForResis,FR)
  t2 <- Sys.time()
  print(t2-t1)
  t1 <- Sys.time()
  print("Reverse Resistant POMs")
  MotifSelecRevResis=CheckDifferencesInDistributions(Config,MotAftTesRevResis,LenMotif2RevResis,RR)
  
  #Log
  # line <- "Kolmogorov done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  # Added:
  FDRForProne = ExtractFDR(Config,MotifSelecForProne,LenMotif2ForProne)
  FDRRevProne = ExtractFDR(Config,MotifSelecRevProne,LenMotif2RevProne)
  FDRForResis = ExtractFDR(Config,MotifSelecForResis,LenMotif2ForResis)
  FDRRevResis = ExtractFDR(Config,MotifSelecRevResis,LenMotif2RevResis)
  
  #Log
  # line <- "ExtractFDR done "
  # write(line, file=Config$LogFile, append=TRUE)
  
  ListFDR = list(FDRForProne,FDRRevProne,FDRForResis,FDRRevResis)
  
  VecFDR=unlist(ListFDR)
  MeaFDR=mean(VecFDR)
  # Added: end
  
  ###Get lengths of existing motifs
  
  LenMotif3ForProne = GetLen3(Config,MotifSelecForProne,LenMotif2ForProne)
  LenMotif3RevProne = GetLen3(Config,MotifSelecRevProne,LenMotif2RevProne)
  LenMotif3ForResis = GetLen3(Config,MotifSelecForResis,LenMotif2ForResis)
  LenMotif3RevResis = GetLen3(Config,MotifSelecRevResis,LenMotif2RevResis)
  
  #Log
  # line <- "Lengths3 done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  # save(Config,
  # Seq_Gene,
  # MeaFDR,
  # FP,
  # RP,
  # FR,
  # RR,
  # ListFDR,
  # MotifSelecForProne,
  # LenMotif3ForProne,
  # MotifSelecRevProne,
  # LenMotif3RevProne,
  # MotifSelecForResis,
  # LenMotif3ForResis,
  # MotifSelecRevResis,
  # LenMotif3RevResis,
  # file=paste("beforeAssign28-08-2019Disp",Config$X,".RData",sep=""))
  
  if(Config$AssignGenes==TRUE){

    ##Get CooChr
    print("Getting motif coordinates")
    t1 <- Sys.time()
    CooChrForProne=GetCooChr(Config,MotifSelecForProne,LenMotif3ForProne)
    #Log
    # line <- "First getCooChr done"
    # write(line,file=Config$LogFile,append=TRUE)
    CooChrRevProne=GetCooChr(Config,MotifSelecRevProne,LenMotif3RevProne)
    #Log
    # line <- "Second getCooChr done"
    # write(line,file=Config$LogFile,append=TRUE)
    CooChrForResis=GetCooChr(Config,MotifSelecForResis,LenMotif3ForResis)
    #Log
    # line <- "Third getCooChr done"
    # write(line,file=Config$LogFile,append=TRUE)
    CooChrRevResis=GetCooChr(Config,MotifSelecRevResis,LenMotif3RevResis)
    t2 <- Sys.time()
    print(t2-t1)

    #Log
    # line <- "Fourth getCooChr done"
    # write(line,file=Config$LogFile,append=TRUE)
    # 
    # #Log
    # line <- "GetCooChr done"
    # write(line,file=Config$LogFile,append=TRUE)
    
    ###Match Motifs to Genes
    print("Assigning genes to motifs")
    t1 <- Sys.time()
    GeneIDForProne=MatchMotifToGene(Config,Seq_Gene,CooChrForProne,LenMotif3ForProne)
    GeneIDRevProne=MatchMotifToGene(Config,Seq_Gene,CooChrRevProne,LenMotif3RevProne)
    GeneIDForResis=MatchMotifToGene(Config,Seq_Gene,CooChrForResis,LenMotif3ForResis)
    GeneIDRevResis=MatchMotifToGene(Config,Seq_Gene,CooChrRevResis,LenMotif3RevResis)
    t2 <- Sys.time()
    print(t2-t1)
    #Log
    # line <- "MatchMotifs done"
    # write(line,file=Config$LogFile,append=TRUE)
    
    print("Gene assignment ended")
    
  }
  
  #if(!is.na(Config$ThrFDR)){
  #  
  #  if(MeaFDR>Config$ThrFDR){
  #    
  #    Config$DrawLogo=TRUE
  #    Config$AssignGenes=TRUE
  #    
  #    if(file.exists(Config$PathOutput)) {
  #      
  #      #delete old structure of directories
  #      unlink(Config$PathOutput,recursive=TRUE)
  #      
  #      #create new structure of directories
  #      dir.create(Config$PathOutput)
  #      dir.create(Config$PathInputDataFileDir)
  #      dir.create(Config$PathFiguresDir)
  #      dir.create(Config$Logos)
  #      dir.create(Config$Scan)
  #      dir.create(Config$HtmlDir)
  #      file.create(Config$IndexHtml)
  #      file.create(Config$ForProneHtml)
  #      file.create(Config$RevProneHtml)
  #      file.create(Config$ForResisHtml)
  #      file.create(Config$RevResisHtml)
  #    }
  #  }
  
  
  #####Logos and Histograms
  
    print("Drawing logos")

    DifResultForProne = LogoScanPlot(Config,MotifSelecForProne,LenMotif3ForProne,FP)
    DifResultRevProne = LogoScanPlot(Config,MotifSelecRevProne,LenMotif3RevProne,RP)
    DifResultForResis = LogoScanPlot(Config,MotifSelecForResis,LenMotif3ForResis,FR)
    DifResultRevResis = LogoScanPlot(Config,MotifSelecRevResis,LenMotif3RevResis,RR)

  if(Config$DrawLogo==TRUE){
    
    ###Html Reports
  
    ReportHtml(Config,MotifSelecForProne,GeneIDForProne,LenMotif3ForProne,FP)
    ReportHtml(Config,MotifSelecRevProne,GeneIDRevProne,LenMotif3RevProne,RP)
    ReportHtml(Config,MotifSelecForResis,GeneIDForResis,LenMotif3ForResis,FR)
    ReportHtml(Config,MotifSelecRevResis,GeneIDRevResis,LenMotif3RevResis,RR)
  
    ###Index Html Report
  
    IndexHtmlReport(Config)
  }
  
  #Log
  # line <- "LogosandHistograms done"
  # write(line,file=Config$LogFile,append=TRUE)
  
  
  # Added:
  # Save motifs' scores, POMs and FDR values.
  save(Config,
       MotifSelecForProne,
       MotifSelecRevProne,
       MotifSelecForResis,
       MotifSelecRevResis,
       GeneIDForProne,
       GeneIDRevProne,
       GeneIDForResis,
       GeneIDRevResis,
       SeqMetFreForProne,
       SeqMetFreForResis,
       SeqMetFreRevProne,
       SeqMetFreRevResis,
       CooChrForProne,
       CooChrRevProne,
       CooChrForResis,
       CooChrRevResis,
       DifResultForProne,
       DifResultRevProne,
       DifResultForResis,
       DifResultRevResis,
       file=paste(Config$PathOutput,"/resultsDisp",Config$X,".RData",sep=""))
       
  print("Finished and data saved.")
  
  return(ListFDR)
}

