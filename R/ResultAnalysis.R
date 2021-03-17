ExtractFDRMerge=function(MotAftTes){
  
  # Extract FDR column from the structure list.
  
  FDR=list()
  i=0
  if(length(MotAftTes)!=0){
    #for(w in 1:min(length(MotAftTes),30)){
    for(w in 1:length(MotAftTes)){
      i=i+1
      if(!is.null(MotAftTes[[w]])){
        FDR[[i]]=MotAftTes[[w]][[5]]
      }
    }
  }
  
  return(FDR)
}

ExtractBndMerge=function(MotAftTes){
  
  # Extract bind scores' columns from the structure list.
  
  Bnd=list()
  Bnd$Prn = list()
  Bnd$Res = list()
  i=1
  if(length(MotAftTes)!=0){
    #for(w in 1:min(length(MotAftTes),30)){
    for(w in 1:length(MotAftTes)){
      if(!is.null(MotAftTes[[w]])){
        Bnd$Prn[[i]]=MotAftTes[[w]][[3]]
        Bnd$Res[[i]]=MotAftTes[[w]][[4]]
      }
      i = i + 1
    }
  }
  
  return(Bnd)
}

dataIsVoidFdr=function(FdrStr){
  res = TRUE
  #print(paste("length ", length(FdrStr$Data)))
  for(i in seq(1, length(FdrStr$Data), 1)){
    #print(length(FdrStr$Data[[i]]))
    if (length(FdrStr$Data[[i]])>1) {
      found = FALSE
      for (val in FdrStr$Data[[i]]){
        if ((!is.null(val)) && (length(val)>0)){
          found = TRUE 
          break
        }
      }
      if (found){
        res = FALSE
        break
      }
    }
  }
  return(res)
}

HMDispDFDR=function(FdrStr, Type, cell_name){
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  #FdrHstMat <- matrix(nrow=length(seq(0,0.05,0.005))-1)
  FdrHstMat <- matrix(nrow=length(seq(0, 0.05, 0.0025))-1)
  #FdrHstMat <- matrix(nrow=length(seq(0,0.02,0.001))-1)
  #FdrHstMat <- matrix(nrow=length(seq(0,0.00175,0.0025))-1)
  
  DspNms <- list()
  for(w in 1:length(FdrStr)) {
    
    #print(paste("Length ", w))
    if (!dataIsVoidFdr(FdrStr[[w]])){
      
      #print(FdrStr[[w]]$Data)
      FdrValVec <- unlist(FdrStr[[w]]$Data)
      # FdrValVec = FdrValVec[-which(FdrValVec>0.02)]
      TotCnt <- length(FdrValVec)
    
      #print(paste("TotCnt ", TotCnt))
      # Classify fdr values in breaks using a histogram.
      #if (TotCnt>0) {
        #FdrHst <- hist(FdrValVec, breaks=seq(0,0.05,0.005), plot = FALSE)
        FdrHst <- hist(FdrValVec, breaks=seq(0,0.05,0.0025), plot = FALSE)
        
        # Get counts for each FDR value range.
        FdrApp <- FdrHst$counts
      
      #}
      
      
      
      # Calculate the frequency
      FdrApp <- unlist(lapply(FdrApp, function(x) x/TotCnt ))
      
      #dim(FdrApp) <- c(length(seq(0,0.05,0.005))-1, 1)
      dim(FdrApp) <- c(length(seq(0,0.05,0.0025))-1, 1)
      # dim(FdrApp) <- c(length(seq(0,0.02,0.001))-1, 1)
      # Append values to the structure.
      DspNms <- c(DspNms, FdrStr[[w]]$Dsp)
      FdrHstMat <- cbind(FdrHstMat, FdrApp)
    }
  } 
  
  FdrHstMat <- FdrHstMat[,-1]
  RngNms <- "[0, 0.0025]"
  #for (i in seq(0.0025,0.0475,0.0025)){
  for (i in seq(0.0025,0.0475,0.0025)) {
    Nam <- paste("(",i,", ",(i+0.0025),"]",sep="")
    Nam <- paste("(",i,", ",(i+0.0025),"]",sep="")
    RngNms <- c(RngNms, Nam)
  }
  # RngNms <- "[0, 0.005]"
  # for (i in seq(0.005,0.045,0.005)){
  #   Nam <- paste("(",i,", ",(i+0.005),"]",sep="")
  #   Nam <- paste("(",i,", ",(i+0.005),"]",sep="")
  #   RngNms <- c(RngNms, Nam)
  # }
  
  
  
  rownames(FdrHstMat) <- RngNms
  colnames(FdrHstMat) <- DspNms
  
  FdrHstMat.m = melt(FdrHstMat)
  
  p <- (ggplot(FdrHstMat.m, aes(Var2, Var1)) + geom_tile(aes(fill = value),
                                                         colour = "white") + scale_fill_gradient(low="lightyellow1", high="red3")) + 
    #theme_void() + # Empty theme without axis lines and texts
    theme(
      # panel.background = element_rect(fill = "transparent", colour = NA),
      # plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_x_continuous(position = "top", expand = c(0, 0)) +
    xlab("Number of bases of displacement from target") +
    ylab("Ranges of FDR values") +
    # ggtitle(paste("Mean FDR values' frequency for each displacement (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    ggtitle(paste(cell_name, " (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 10), size=12),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 30, l = 10), size=12),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 20, l = 0), size=12),
          legend.title = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0), size=12)) +
    labs(fill = "Frequency") +
    theme(text=element_text(family="Arial", face="bold", size=12))
  
  #print(p)
  ggsave(paste(cell_name, "_FDR_", TypNme, ".png"), plot = last_plot(), width=7, height=5, dpi=300)
  
  return(TRUE)
  
}

HMFDRDispLog=function(FdrStr, option, Type){
  
  # Options:
  # H: expand high FDR values.
  # L: Expand low FDR values.
  # HL: Expan high and low values.
  
  option=tolower(option)
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  
  # Transform data to logarithm.
  LogFdrLst <- list()
  DspNms <- list()
  
  maxFDR <- -Inf
  
  for(i in 1:length(FdrStr)) {
    
        
    if (!dataIsVoidFdr(FdrStr[[i]])){
      FdrValVec <- unlist(FdrStr[[i]]$Data)
      FdrValVec <- unlist(lapply(FdrValVec, function(x) (log(x)*-1)))
      FdrValVec <- lapply(FdrValVec, function(x) if (x=="Inf") 1500 else x)
      LogFdrLst[[i]] <- FdrValVec
      DspNms <- c(DspNms, FdrStr[[i]]$Dsp)
      maxFDR <- max(maxFDR, max(unlist(FdrValVec)))
    }
  }
  
  if (option=="h"){
    # Generate breaks.
    if (maxFDR < 10){
      BrkLst <- c(0,seq(3.25,ceiling(maxFDR)+1,0.25))
    } 
    else{
      BrkLst <- c(0, seq(3.25,10,0.25), ceiling(maxFDR)+1)
    }
  }
  else{
    if (option=="l"){
      BrkLst <- c(0, 10, seq(50,ceiling(maxFDR)+1,50))
    }
    else{
      BrkLst <- c(0, seq(3.25, 5, 0.25), 10, seq(50,ceiling(maxFDR)+1,50))
    }
  }
  
  
  # Generate the matrix for the heatmap.
  FdrHstMat <- matrix(nrow=length(BrkLst)-1)
  #print(BrkLst)
  
  for(w in 1:length(LogFdrLst)) {
    
    #print(paste("Length ", w))
    if ((!is.null(LogFdrLst[[w]])) & (length(LogFdrLst[[w]])>0)){
      
      #print(length(LogFdrLst[[w]]))
      FdrValVec <- unlist(LogFdrLst[[w]])
      TotCnt <- length(FdrValVec)
      
      # Classify fdr values in breaks using a histogram.
      #save(FdrValVec, file="FdrValVec.RData")
      FdrHst <- hist(FdrValVec, breaks=BrkLst, plot = FALSE)
      
      # Get counts for each FDR value range.
      FdrApp <- FdrHst$counts
      
      # Calculate the frequency
      FdrApp <- unlist(lapply(FdrApp, function(x) x/TotCnt ))
      
      dim(FdrApp) <- c(length(BrkLst)-1, 1)
      # Append values to the structure.
      FdrHstMat <- cbind(FdrHstMat, FdrApp)
    }
    # else {
    #   empty_list <- vector(mode = "list", length = length(BrkLst)-1)
    #   dim(empty_list) <- c(length(BrkLst)-1, 1)
    #   FdrHstMat <- cbind(FdrHstMat, empty_list)
    # }
  } 
  
  BrkLst[length(BrkLst)]=BrkLst[length(BrkLst)]-1
  FdrHstMat <- FdrHstMat[,-1]
  
  if (option=="h"){
    RngNms <- "(3, 3.25]"
    for (i in 2:(length(BrkLst)-1)){
      
      if (i==(length(BrkLst)-1))  { Nam <- paste("(",BrkLst[i],", ",ceiling(maxFDR),"]",sep="") }
      else { Nam <- paste("(",BrkLst[i],", ",(BrkLst[i]+0.25),"]",sep="") }
      RngNms <- c(RngNms, Nam)
    }
  }
  else{
    if (option=="l"){
      RngNms <- "(3, 10]"
      for (i in 2:(length(BrkLst)-1)){
        Nam <- paste("(",BrkLst[i],", ",(BrkLst[i]+0.25),"]",sep="")
        RngNms <- c(RngNms, Nam)
      }
    }
    else{
      RngNms <- "(3, 3.25]"
      for (i in 2:(length(BrkLst)-1)){
        Nam <- paste("(",BrkLst[i],", ",(BrkLst[i+1]),"]",sep="")
        RngNms <- c(RngNms, Nam)
      }
    }
  }
  
  rownames(FdrHstMat) <- RngNms
  
  print(paste("Displacements' length", length(DspNms)))
  print(paste("Matrix dimension ", dim(FdrHstMat)))
  colnames(FdrHstMat) <- DspNms
  
  FdrHstMat.m = melt(FdrHstMat)
  p <- (ggplot(FdrHstMat.m, aes(Var2, Var1)) + geom_tile(aes( fill = value),
                                                         colour = "white") + scale_fill_gradient(low = "lightyellow",
                                                                                                 high = "red3")) +
    #theme_void() + # Empty theme without axis lines and texts
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    xlab("Number of bases of displacement from target") +
    ylab("-log(μFDR)") +
    ggtitle(paste("Logarithmic FDR values' frequency for each displacement (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
          legend.title = element_text(margin = margin(t = 0, r = 0, b = 10, l = 0))) +
    labs(fill = "Frequency") +
    theme(text=element_text(family="Times New Roman", face="bold", size=12))
  
  print(p)
  
  return(TRUE)
  
}

HSDispFdr <- function(FdrStr, CntStr, Type, nPomLimUp=4000){
  
  # Type:
  # "p": prone motifs.
  # "r": resistant motifs.
  # "pr": prone and resistant motifs.
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  
  MeaFdr <- list()
  DspNms <- list()
  FdrVar <- list()
  FdrCntLst <- list()
  
  j = 1 
  for(w in 1:length(FdrStr)) {
    
    #print(paste("length", Idx))
    if (!dataIsVoidFdr(FdrStr[[w]])){
      
      #print(FdrStr[[w]]$Data)
      FdrValVec <- unlist(FdrStr[[w]]$Data)
      FdrCntLst[j] <- switch(tolower(Type),
                             "pr" = sum(CntStr[[w]]$NumForPrn, CntStr[[w]]$NumForRes, CntStr[[w]]$NumRevPrn, CntStr[[w]]$NumRevRes),
                             "p" = sum(CntStr[[w]]$NumForPrn, CntStr[[w]]$NumRevPrn),
                             "r" = sum(CntStr[[w]]$NumForRes, CntStr[[w]]$NumRevRes))
      
      
      MeaFdr[j] <- mean(FdrValVec)
      FdrVar[j] <- var(FdrValVec)
      #print(paste("Idx", Idx))
      #print(BndValVec)
      # print(paste("Idx", Idx, "mean", mean(BndValVec), "var", var(unlist(BndValVec))))
      DspNms[j] <- FdrStr[[w]]$Dsp
      j = j + 1
      
    }
  }
  
  fdr.mt <- cbind(unlist(DspNms), unlist(MeaFdr), unlist(FdrVar), unlist(FdrCntLst))
  #bnd.df <- as.data.frame(bnd.mt)
  fdr.df <- as.data.frame(fdr.mt)
  #print(bnd.df)
  # 
  # p <- (ggplot(data=bnd.df) +
  #    geom_bar(stat="identity", aes(x=bnd.df$displacement, y=bnd.df$scores)))
  # p<-ggplot(data=bnd.df, aes(x=V1, y=V2)) +
  #   geom_bar(stat="identity", fill="indianred4")+
  #   geom_errorbar(aes(ymin=V2-V3, ymax=V2+V3), width=.8,
  #                 position=position_dodge(.9)) +
  #   theme_minimal() +
  #   theme(axis.title.y = element_text(margin = margin(t = 5, r = 20, b = 0, l = 10)),
  #         axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 10, l = 0)),
  #         plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 25, l = 0))) +
  #   xlab("Number of bases of displacement from target") +
  #   ylab("Mean of the difference of binding scores") +
  #   ggtitle("Difference between resistant and prone binding scores against motif displacement", subtitle = waiver())
  
  POMUpLim <- nPomLimUp
  #POMUpLim <- 500
  
  # POMBoLim <- 0
  POMBoLim <- 0
  
  FdrUpLim <- 0.0045
  #BndUpLim <- 1.5
  FdrLoLim <- 0
  
  p<-ggplot(fdr.df) +
    geom_bar(mapping=aes(x = V1, y = V2*POMUpLim/FdrUpLim), stat="identity", fill="mediumpurple", color="mediumpurple4") +
    geom_errorbar(aes(x=V1, ymin=(V2-V3)*(POMUpLim/FdrUpLim), ymax=(V2+V3)*(POMUpLim/FdrUpLim)), width=1.6,
                  position=position_dodge(.9)) +
    geom_line(aes(x=V1, y=V4), size=1, color="orange1") +
    scale_y_continuous(name = "Number of filtered POMs", limits = c(0,POMUpLim),
                       sec.axis = sec_axis(~ . * FdrUpLim / POMUpLim , name = "Mean FDR value")) +
    theme_minimal() +
    xlab("Number of bases of displacement from target") +
    theme(axis.title.y.left = element_text(margin = margin(t = 5, r = 10, b = 0, l = 10), colour = "darkorange4"),
          axis.line.y.left = element_line(color="darkorange4"),
          axis.text.y.left = element_text(color="darkorange4"),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
          axis.title.y.right = element_text(margin = margin(t = 5, r = 10, b = 0, l = 10), colour = "purple4"),
          axis.line.y.right = element_line(color="purple4"),
          axis.text.y.right = element_text(color="purple4"),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 20, l = 0))) +
    # ggtitle(paste("Difference between resistant and prone\nbinding scores against motif displacement \n(", TypNme, " motifs, " , senseName, ")", sep=""), subtitle = waiver())+
    ggtitle(paste("IPSC (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    theme(text=element_text(family="Arial", face="bold", size=12))
  
  
  print(p)
  
  return(TRUE)
}

HSDispBnd <- function(BndStr, CntStr, Type, cell_name, nPomLimUp=6000, direction = 0, max_dif = 0.6){
  
  # Type:
  # "p": prone motifs.
  # "r": resistant motifs.
  # "pr": prone and resistant motifs.
  
  if (direction==0){
    senseName <- "positive and negative strands"
  }
  else {
    if (direction<0){
      senseName <- "negative strand"
    }
    else {
      senseName <- "positive strand"
    }
    
  }
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  
  MeaBnd <- list()
  DspNms <- list()
  BndVar <- list()
  BndCntLst <- list()
  
  j = 1 
  for(Idx in 1:length(BndStr)) {
    
    #print(paste("length", Idx))
    if (!is.null(BndStr[[Idx]]) && (length(BndStr[[Idx]]$Data)>0)){
      
      if (direction==0){
        BndValVec <- unlist(BndStr[[Idx]]$Data)
        BndCntLst[j] <- switch(tolower(Type),
                               "pr" = sum(CntStr[[Idx]]$NumForPrn, CntStr[[Idx]]$NumForRes, CntStr[[Idx]]$NumRevPrn, CntStr[[Idx]]$NumRevRes),
                               "p" = sum(CntStr[[Idx]]$NumForPrn, CntStr[[Idx]]$NumRevPrn),
                               "r" = sum(CntStr[[Idx]]$NumForRes, CntStr[[Idx]]$NumRevRes))
      }
      if (direction>0){
        BndValVec <- switch(tolower(Type),
                               "pr" = unlist(c(BndStr[[Idx]]$Data[[1]], BndStr[[Idx]]$Data[[3]])),
                               "p" = unlist(c(BndStr[[Idx]]$Data[[1]])),
                               "r" = unlist(c(BndStr[[Idx]]$Data[[1]])))
        BndCntLst[j] <- switch(tolower(Type),
                               "pr" = sum(CntStr[[Idx]]$NumForPrn, CntStr[[Idx]]$NumForRes),
                               "p" = sum(CntStr[[Idx]]$NumForPrn),
                               "r" = sum(CntStr[[Idx]]$NumForRes))
      }
      if (direction<0){
        BndValVec <- switch(tolower(Type),
                            "pr" = unlist(c(BndStr[[Idx]]$Data[[2]], BndStr[[Idx]]$Data[[4]])),
                            "p" = unlist(c(BndStr[[Idx]]$Data[[2]])),
                            "r" = unlist(c(BndStr[[Idx]]$Data[[2]])))
        BndCntLst[j] <- switch(tolower(Type),
                               "pr" = sum(CntStr[[Idx]]$NumRevPrn, CntStr[[Idx]]$NumRevRes),
                               "p" = sum(CntStr[[Idx]]$NumRevPrn),
                               "r" = sum(CntStr[[Idx]]$NumRevRes))
      }
      MeaBnd[j] <- mean(BndValVec)
      BndVar[j] <- var(BndValVec)
      #print(paste("Idx", Idx))
      #print(BndValVec)
      # print(paste("Idx", Idx, "mean", mean(BndValVec), "var", var(unlist(BndValVec))))
      DspNms[j] <- BndStr[[Idx]]$Dsp
      j = j + 1
      
    }
  }
  
  bnd.mt <- cbind(unlist(DspNms), unlist(MeaBnd), unlist(BndVar), unlist(BndCntLst))
  #bnd.df <- as.data.frame(bnd.mt)
  bnd.df <- as.data.frame(bnd.mt)
  #print(bnd.df)
  # 
  # p <- (ggplot(data=bnd.df) +
  #    geom_bar(stat="identity", aes(x=bnd.df$displacement, y=bnd.df$scores)))
  # p<-ggplot(data=bnd.df, aes(x=V1, y=V2)) +
  #   geom_bar(stat="identity", fill="indianred4")+
  #   geom_errorbar(aes(ymin=V2-V3, ymax=V2+V3), width=.8,
  #                 position=position_dodge(.9)) +
  #   theme_minimal() +
  #   theme(axis.title.y = element_text(margin = margin(t = 5, r = 20, b = 0, l = 10)),
  #         axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 10, l = 0)),
  #         plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 25, l = 0))) +
  #   xlab("Number of bases of displacement from target") +
  #   ylab("Mean of the difference of binding scores") +
  #   ggtitle("Difference between resistant and prone binding scores against motif displacement", subtitle = waiver())
  
  POMUpLim <- nPomLimUp
  #POMUpLim <- 500
  
  # POMBoLim <- 0
  POMBoLim <- 0
  
  BndUpLim <- max_dif
  #BndUpLim <- 1.5
  BndLoLim <- 0
  
  p<-ggplot(bnd.df) +
    geom_bar(mapping=aes(x = V1, y = V2*POMUpLim/BndUpLim), stat="identity", fill="mediumpurple", color="mediumpurple4") +
    geom_errorbar(aes(x=V1, ymin=(V2-V3)*(POMUpLim/BndUpLim), ymax=(V2+V3)*(POMUpLim/BndUpLim)), width=1.6,
                  position=position_dodge(.9)) +
    geom_line(aes(x=V1, y=V4), size=1, color="orange1") +
    scale_y_continuous(name = "Number of filtered POMs", limits = c(0,POMUpLim),
                       sec.axis = sec_axis(~ . * BndUpLim / POMUpLim , name = "Mean of the difference of matching scores")) +
    theme_minimal() +
    xlab("Number of bases of displacement from target") +
    theme(axis.title.y.left = element_text(margin = margin(t = 5, r = 10, b = 0, l = 10), colour = "darkorange4"),
          axis.line.y.left = element_line(color="darkorange4"),
          axis.text.y.left = element_text(color="darkorange4"),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
          axis.title.y.right = element_text(margin = margin(t = 5, r = 10, b = 0, l = 10), colour = "purple4"),
          axis.line.y.right = element_line(color="purple4"),
          axis.text.y.right = element_text(color="purple4"),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 20, l = 0))) +
    # ggtitle(paste("Difference between resistant and prone\nbinding scores against motif displacement \n(", TypNme, " motifs, " , senseName, ")", sep=""), subtitle = waiver())+
    ggtitle(paste(cell_name, " (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
              theme(text=element_text(family="Arial", face="bold", size=12))
  
  ggsave(paste(cell_name, "_CNTBND_", TypNme, ".png"), plot = last_plot())
  #print(p)
  
  return(TRUE)
}


HSDispBndNotCount <- function(BndStr, Type, cell_name, direction = 0, max_dif = 0.6){
  
  # Type:
  # "p": prone motifs.
  # "r": resistant motifs.
  # "pr": prone and resistant motifs.
  
  if (direction==0){
    senseName <- "positive and negative strands"
  }
  else {
    if (direction<0){
      senseName <- "negative strand"
    }
    else {
      senseName <- "positive strand"
    }
    
  }
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  
  MeaBnd <- list()
  DspNms <- list()
  BndVar <- list()
  BndCntLst <- list()
  
  j = 1 
  for(Idx in 1:length(BndStr)) {
    
    #print(paste("length", Idx))
    if (!is.null(BndStr[[Idx]]) && (length(BndStr[[Idx]]$Data)>0)){
      
      if (direction==0){
        BndValVec <- unlist(BndStr[[Idx]]$Data)
      }
      if (direction>0){
        BndValVec <- switch(tolower(Type),
                            "pr" = unlist(c(BndStr[[Idx]]$Data[[1]], BndStr[[Idx]]$Data[[3]])),
                            "p" = unlist(c(BndStr[[Idx]]$Data[[1]])),
                            "r" = unlist(c(BndStr[[Idx]]$Data[[1]])))
      }
      if (direction<0){
        BndValVec <- switch(tolower(Type),
                            "pr" = unlist(c(BndStr[[Idx]]$Data[[2]], BndStr[[Idx]]$Data[[4]])),
                            "p" = unlist(c(BndStr[[Idx]]$Data[[2]])),
                            "r" = unlist(c(BndStr[[Idx]]$Data[[2]])))
      }
      MeaBnd[j] <- mean(BndValVec)
      BndVar[j] <- var(BndValVec)
      #print(paste("Idx", Idx))
      #print(BndValVec)
      # print(paste("Idx", Idx, "mean", mean(BndValVec), "var", var(unlist(BndValVec))))
      DspNms[j] <- BndStr[[Idx]]$Dsp
      j = j + 1
      
    }
  }
  
  bnd.mt <- cbind(unlist(DspNms), unlist(MeaBnd), unlist(BndVar))
  #bnd.df <- as.data.frame(bnd.mt)
  bnd.df <- as.data.frame(bnd.mt)
  #print(bnd.df)
  # 
  # p <- (ggplot(data=bnd.df) +
  #    geom_bar(stat="identity", aes(x=bnd.df$displacement, y=bnd.df$scores)))
  # p<-ggplot(data=bnd.df, aes(x=V1, y=V2)) +
  #   geom_bar(stat="identity", fill="indianred4")+
  #   geom_errorbar(aes(ymin=V2-V3, ymax=V2+V3), width=.8,
  #                 position=position_dodge(.9)) +
  #   theme_minimal() +
  #   theme(axis.title.y = element_text(margin = margin(t = 5, r = 20, b = 0, l = 10)),
  #         axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 10, l = 0)),
  #         plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 25, l = 0))) +
  #   xlab("Number of bases of displacement from target") +
  #   ylab("Mean of the difference of binding scores") +
  #   ggtitle("Difference between resistant and prone binding scores against motif displacement", subtitle = waiver())
  
  #POMUpLim <- 500
  
  # POMBoLim <- 0
  
  BndUpLim <- max_dif
  #BndUpLim <- 1.5
  BndLoLim <- 0
  
  p<-ggplot(bnd.df) +
    geom_bar(mapping=aes(x = V1, y = V2), stat="identity", fill="mediumpurple", color="mediumpurple4") +
    geom_errorbar(aes(x=V1, ymin=(V2-V3), ymax=(V2+V3)), width=1.6,
                  position=position_dodge(.9)) +
    scale_y_continuous( name = "Mean of the difference of matching scores", limits = c(BndLoLim,BndUpLim)) +
    theme_minimal() +
    xlab("Number of bases of displacement from target") +
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0)),
          axis.title.y.left = element_text(margin = margin(t = 5, r = 10, b = 0, l = 10), colour = "purple4"),
          axis.line.y.left = element_line(color="purple4"),
          axis.text.y.left = element_text(color="purple4"),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 15, r = 0, b = 20, l = 0))) +
    # ggtitle(paste("Difference between resistant and prone\nbinding scores against motif displacement \n(", TypNme, " motifs, " , senseName, ")", sep=""), subtitle = waiver())+
    ggtitle(paste(cell_name, " (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    theme(text=element_text(family="Arial", face="bold", size=12))
  
  ggsave(paste(cell_name, "_CNTBND_", TypNme, ".png"), plot = last_plot())
  #print(p)
  
  return(TRUE)
}


HMDispBnd <- function(BndStr, Type){
  
  TypNme <- switch(tolower(Type),
                   "p"="Prone",
                   "r"="Resistant",
                   "pr"="Prone and resistant")
  
  MaxBnd <- 0
  MinBnd <- Inf
  
  MeaBnd <- list()
  DspNms <- list()
  BndVar <- list()
  
  # Find maximum and minum matching score difference values.
  for(Idx in 1:length(BndStr)) {
    
    if (!is.null(BndStr[[Idx]])  && (length(BndStr[[Idx]]$Data)>0)){
      
      MaxBnd <- max(max(unlist(BndStr[[Idx]]$Data)),MaxBnd)
      MinBnd <- min(min(unlist(BndStr[[Idx]]$Data)),MinBnd)
      
    }
  }

  MaxBnd <- ceiling_dec(MaxBnd, 1)

  MinBndAux <- round(MinBnd, 1)
  while (MinBndAux>MinBnd) MinBndAux = MinBndAux - 0.05
  MinBnd <- MinBndAux
  BrkLst <- seq(MinBnd, MaxBnd, 0.05)
  
  
  RngNms <- paste("[",MinBnd,", ",BrkLst[2],"]",sep="")
  for (i in 2:(length(BrkLst)-1)){
    Nam <- paste("(",BrkLst[i],", ",BrkLst[i+1],"]",sep="")
    RngNms <- c(RngNms, Nam)
  }
  
  BndHstMat <- matrix(nrow=length(BrkLst)-1)
  
  # Genertae frequencies for each displacement.
  for(Idx in 1:length(BndStr)) {
    
    #print(Idx)
    if (!is.null(BndStr[[Idx]]) && (length(BndStr[[Idx]]$Data)>0)){
      
      BndValVec <- unlist(BndStr[[Idx]]$Data)
      TotCnt <- length(BndValVec)
      
      # Classify fdr values in breaks using a histogram.
      BndHst <- hist(BndValVec, breaks=BrkLst, plot = FALSE)
      
      # Get counts for each FDR value range.
      BndApp <- BndHst$counts
      
      # Calculate the frequency
      BndApp <- unlist(lapply(BndApp, function(x) x/TotCnt ))
      
      dim(BndApp) <- c(length(BrkLst)-1, 1)
      # Append values to the structure.
      DspNms <- c(DspNms, BndStr[[Idx]]$Dsp)
      BndHstMat <- cbind(BndHstMat, BndApp)
    }
  }
  BndHstMat <- BndHstMat[,-1]
  rownames(BndHstMat) <- RngNms
  colnames(BndHstMat) <- DspNms
  
  BndHstMat.m = melt(BndHstMat)
  
  
  p <- (ggplot(BndHstMat.m, aes(Var2, Var1)) + geom_tile(aes(fill = value),
                                                         colour = "white") + scale_fill_gradient(low = "lightyellow",
                                                                                                 high = "red3")) +
    
    #theme_void() + # Empty theme without axis lines and texts
    theme(
      # panel.background = element_rect(fill = "transparent", colour = NA),
      # plot.background = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA)
    ) +
    scale_x_continuous(position = "top", expand = c(0, 0)) +
    xlab("Number of bases of displacement from target") +
    ylab("Difference between prone and resistant matching scores") +
    ggtitle(paste("IPSC (", TypNme, " motifs)", sep=""), subtitle = waiver()) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 10), size=12),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 20, l = 10), size=12),
          plot.title = element_text(hjust = 0.5, margin = margin(t = 0, r = 0, b = 20, l = 0), size=12),
          legend.title = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0), , size=12)) +
          labs(fill = "Frequency") +
          # theme(text=element_text(family="Times New Roman", face="bold", size=12))
          theme(text=element_text(face="bold", size=12))
  
  print(p)
  
  return(TRUE)
}

# HMBnd <- function(BndStr) {
#   
#   
#   
#   
#   
#   
#   
# }

CountPOMs <- function(GlbEst) {
  
  CntStr=list()
  
  for (Idx in 1:length(GlbEst)){
    
    
    CntStr[[Idx]] <- list()
    DspData <-  GlbEst[[Idx]]$Data
    CntStr[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
    
    CntStr[[Idx]]$NumForPrn <- 0
    CntStr[[Idx]]$NumForRes <- 0
    CntStr[[Idx]]$NumRevPrn <- 0
    CntStr[[Idx]]$NumRevRes <- 0
    
    
    if(length(DspData[[1]])!=0){
      for (w in 1:length(DspData[[1]])){
        CntStr[[Idx]]$NumForPrn =  CntStr[[Idx]]$NumForPrn + length(DspData[[1]][[w]][[1]])
      }
    }
    
    if(length(DspData[[3]])!=0){
      for (w in 1:length(DspData[[3]])){
        CntStr[[Idx]]$NumForRes =  CntStr[[Idx]]$NumForRes + length(DspData[[3]][[w]][[1]])
      }
    }
    
    if(length(DspData[[2]])!=0){
      for (w in 1:length(DspData[[2]])){
        CntStr[[Idx]]$NumRevPrn =  CntStr[[Idx]]$NumRevPrn + length(DspData[[2]][[w]][[1]])
      }
    }
    
    if(length(DspData[[4]])!=0){
      for (w in 1:length(DspData[[4]])){
        CntStr[[Idx]]$NumRevRes =  CntStr[[Idx]]$NumRevRes + length(DspData[[4]][[w]][[1]])
      }
    }
    
  }
  
  return(CntStr)
  
}

CountPOMsPrn <- function(GlbEst) {
  
  CntStr=list()
  
  for (Idx in 1:length(GlbEst)){
    
    CntStr[[Idx]] <- list()
    DspData <-  GlbEst[[Idx]]$Data
    CntStr[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
    
    CntStr[[Idx]]$NumForPrn <- 0
    CntStr[[Idx]]$NumRevPrn <- 0
    
    
    if(length(DspData[[1]])!=0){
      for (w in 1:length(DspData[[1]])){
        CntStr[[Idx]]$NumForPrn =  CntStr[[Idx]]$NumForPrn + length(DspData[[1]][[w]][[1]])
      }
    }
    
    if(length(DspData[[2]])!=0){
      for (w in 1:length(DspData[[2]])){
        CntStr[[Idx]]$NumRevPrn =  CntStr[[Idx]]$NumRevPrn + length(DspData[[2]][[w]][[1]])
      }
    }
  }
  
  return(CntStr)
  
}

CountPOMsRes <- function(GlbEst) {
  
  CntStr=list()
  
  for (Idx in 1:length(GlbEst)){
    
    CntStr[[Idx]] <- list()
    DspData <-  GlbEst[[Idx]]$Data
    CntStr[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
    
    CntStr[[Idx]]$NumForRes <- 0
    CntStr[[Idx]]$NumRevRes <- 0
    
    
    if(length(DspData[[3]])!=0){
      for (w in 1:length(DspData[[3]])){
        CntStr[[Idx]]$NumForRes =  CntStr[[Idx]]$NumForRes + length(DspData[[3]][[w]][[1]])
      }
    }
    
    if(length(DspData[[4]])!=0){
      for (w in 1:length(DspData[[4]])){
        CntStr[[Idx]]$NumRevRes =  CntStr[[Idx]]$NumRevRes + length(DspData[[4]][[w]][[1]])
      }
    }
    
  }
  
  return(CntStr)
  
}

CountPOMsIdx <- function(MotifInfo, Dsp) {
  
  CntStr=list()
    
    
    CntStr$Dsp <- Dsp
    
    CntStr$NumForPrn <- 0
    CntStr$NumForRes <- 0
    CntStr$NumRevPrn <- 0
    CntStr$NumRevRes <- 0
    
    
    if(length(MotifInfo[[1]])!=0){
      for (w in 1:length(MotifInfo[[1]])){
        CntStr$NumForPrn =  CntStr$NumForPrn + length(MotifInfo[[1]][[w]][[1]])
      }
    }
    
    if(length(MotifInfo[[3]])!=0){
      for (w in 1:length(MotifInfo[[3]])){
        CntStr$NumForRes =  CntStr$NumForRes + length(MotifInfo[[3]][[w]][[1]])
      }
    }
    
    if(length(MotifInfo[[2]])!=0){
      for (w in 1:length(MotifInfo[[2]])){
        CntStr$NumRevPrn =  CntStr$NumRevPrn + length(MotifInfo[[2]][[w]][[1]])
      }
    }
    
    if(length(MotifInfo[[4]])!=0){
      for (w in 1:length(MotifInfo[[4]])){
        CntStr$NumRevRes =  CntStr$NumRevRes + length(MotifInfo[[4]][[w]][[1]])
      }
    }
    
  
  return(CntStr)
  
}

CountPOMsPrnIdx <- function(MotifInfo, Dsp) {
  
  CntStr=list()
  
    DspData <- MotifInfo
    CntStr$Dsp <-Dsp
    
    CntStr$NumForPrn <- 0
    CntStr$NumRevPrn <- 0
    
    
    if(length(DspData[[1]])!=0){
      for (w in 1:length(DspData[[1]])){
        CntStr$NumForPrn =  CntStr$NumForPrn + length(DspData[[1]][[w]][[1]])
      }
    }
    
    if(length(DspData[[2]])!=0){
      for (w in 1:length(DspData[[2]])){
        CntStr$NumRevPrn =  CntStr$NumRevPrn + length(DspData[[2]][[w]][[1]])
      }
    }
  
  return(CntStr)
  
}

CountPOMsResIdx <- function(MotifInfo, Dsp) {
  
  CntStr=list()
  
    DspData <- MotifInfo
    CntStr <- list()
    CntStr$Dsp <- Dsp
    
    CntStr$NumForRes <- 0
    CntStr$NumRevRes <- 0
    
    
    if(length(DspData[[1]])!=0){
      for (w in 1:length(DspData[[1]])){
        CntStr$NumForRes =  CntStr$NumForRes + length(DspData[[1]][[w]][[1]])
      }
    }
    
    if(length(DspData[[2]])!=0){
      for (w in 1:length(DspData[[2]])){
        CntStr$NumRevRes =  CntStr$NumRevRes + length(DspData[[2]][[w]][[1]])
      }
    }
    
  
  return(CntStr)
  
}

# TODO:
# change 30's limit in merge and extracts.


mergeDMMDResults=function(DspLst, InpPath, GenPrefRdt, subDir=TRUE, info_to_save=NULL){
  
  # mergeDMMDResults.
  # 
  # The resulting structure, contains:
  # - The displacement value with the best average FDR.
  # - FDR values, POMs and PWMs of the selected motifs.
  # 
  # Also it creates an html file with the following data:
  # - Bar plot of mean FDR values against displacements.
  # - Bar plot of selected motifs against displacements.
  # - Histogram of motif frequency for FDR values' ranges.
  # - Scatter plot of FDR values against motif lengths.
  # - Bar chart of motif frecuency against motif lengths.
  # - Graphical representation of the selected motifs'
  #   logos', scan plots and their related genes.
  
  if ( (!is.null(info_to_save)) && !is.vector(info_to_save) ){
    info_to_save = c(info_to_save)
  } 
  # info_to_save = c("fdr","cnt","bnd")
  
  MrgRes=list()
  
  # Global structure.
  GlbEst=list()
  
  # Displacements with loaded data.
  LodDsp=list()
  
  if(!dir.exists(InpPath)) stop(paste("The input directory ", InpPath, " does not exist.", sep=""))
  
  
  
  # start structure
  maxDsp <- max(DspLst)
  minDsp <- min(DspLst)
  progression_value <- ( if (length(DspLst)==0){
    0
  }
  else { 
    if (DspLst[1]<DspLst[2]) {
    DspLst[2] - DspLst[1]
    }
    else {
    DspLst[1] - DspLst[2]
    }
  } )
  DspLst <- seq(minDsp, maxDsp, progression_value)
  
  Idx <- 1
  for(Dsp in DspLst){
   GlbEst[[Idx]] = list()
   GlbEst[[Idx]]$Dsp <- Dsp
   GlbEst[[Idx]]$Data <- list()
   Idx <- Idx + 1
  }
  
  
  # Load displacements' RData files.
  print("Loading data for each displacement.")
  Idx <- 1
  
  if ( ( is.null(info_to_save ) ) || ( "cnt" %in% info_to_save ) ) {
    print("Counting total POMs")
    CntStr <- list()
    CntStrPrn <- list()
    CntStrRes <- list()
  }
  
  for(Dsp in DspLst){
    
    if (subDir) {
      NxtFle<-file.path(InpPath, paste("GOUTDisp",Dsp,sep=""), paste(GenPrefRdt, Dsp, ".RData", sep=""))      
    }
    else {
      NxtFle<-file.path(InpPath, paste(GenPrefRdt, Dsp, ".RData", sep=""))
    } 

    if (file.exists(NxtFle)){
      print(paste("Loading displacement ", Dsp, ".", sep=""))
      load(NxtFle)
      
      if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
        FDRForProne <- ExtractFDRMerge(MotifSelecForProne)
        FDRRevProne <- ExtractFDRMerge(MotifSelecRevProne)
        FDRForResis <- ExtractFDRMerge(MotifSelecForResis)
        FDRRevResis <- ExtractFDRMerge(MotifSelecRevResis)
      }
      
      
      if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
        GlbEst[[Idx]]$Data[[5]] <- FDRForProne
        GlbEst[[Idx]]$Data[[6]] <- FDRRevProne
        GlbEst[[Idx]]$Data[[7]] <- FDRForResis
        GlbEst[[Idx]]$Data[[8]] <- FDRRevResis
      }
      if ( ( is.null(info_to_save ) ) || ( "bnd" %in% info_to_save ) ) {
        GlbEst[[Idx]]$Data[[9]] <- DifResultForProne
        GlbEst[[Idx]]$Data[[10]] <- DifResultRevProne
        GlbEst[[Idx]]$Data[[11]] <- DifResultForResis
        GlbEst[[Idx]]$Data[[12]] <- DifResultRevResis
      }
       
      # inverse_displacement_indx = Idx
      # if (Dsp != 0) {
      #   for (j in 1:length(GlbEst)){
      #     if (GlbEst[[j]]$Dsp == - Dsp) {
      #       inverse_displacement_indx = j
      #       rm(MotifSelecRevProne)
      #       rm(MotifSelecRevResis)
      #       break
      #     }
      #   }
      # }
      # print(paste("inverse displacement", inverse_displacement_indx))
    
      
      # if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
      #   GlbEst[[inverse_displacement_indx]]$Data[[6]] <- FDRRevProne
      #   GlbEst[[inverse_displacement_indx]]$Data[[8]] <- FDRRevResis
      # }
      # if ( ( is.null(info_to_save ) ) || ( "bnd" %in% info_to_save ) ) {
      #   GlbEst[[inverse_displacement_indx]]$Data[[10]] <- DifResultRevProne
      #   GlbEst[[inverse_displacement_indx]]$Data[[12]] <- DifResultRevResis
      # }
      
      if ( ( is.null(info_to_save ) ) || ( "cnt" %in% info_to_save ) ) {
        
        
        
        if (subDir) NxtFleAux<-file.path(InpPath, paste("GOUTDisp",-Dsp,sep=""), paste(GenPrefRdt, -Dsp, ".RData", sep=""))
        else NxtFleAux<-file.path(InpPath, paste(GenPrefRdt, -Dsp, ".RData", sep=""))
     
        
        CntStr[[Idx]] <- CountPOMsIdx(list(MotifSelecForProne, MotifSelecRevProne, MotifSelecForResis, MotifSelecRevResis), Dsp)
        CntStrPrn[[Idx]] <- CountPOMsPrnIdx(list(MotifSelecForProne, MotifSelecRevProne), Dsp)
        CntStrRes[[Idx]] <- CountPOMsResIdx(list(MotifSelecForResis, MotifSelecRevResis), Dsp)
      }
      
      rm(MotifSelecForProne)
      rm(MotifSelecRevProne)
      rm(MotifSelecForResis)
      rm(MotifSelecRevResis)
      if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
        rm(FDRForProne)
        rm(FDRRevProne)
        rm(FDRForResis)
        rm(FDRRevResis)
      }
      gc()
      
      print(paste(as.numeric(object.size(GlbEst))/(1024**2), " MB"))
      Idx <- Idx + 1
      
    }
    else{
      print(paste("The input file ", NxtFle, " could not be found.", sep=""))
      print("Its results will not be included.")
    }
  }
  
  if (length(GlbEst)!=0){
    
    ThrFDR = Inf
    MinDsp = GlbEst[[1]]$Idx
    VecMeaFDR = list()
    # Structure with FDR values for each displacement for the following heatmap.
    FdrStr=list()
    FdrStrPrn=list()
    FdrStrRes=list()
    BndStr=list()
    BndStrPrn=list()
    BndStrRes=list()
    
    # Calculate each displacement's mean FDR and get the displacement value with the best FDR.
    
    for (Idx in 1:length(GlbEst)){
      
      DspData <-  GlbEst[[Idx]]$Data
      
      if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
        FdrStr[[Idx]] <- list()
        FdrStrPrn[[Idx]] <- list()
        FdrStrRes[[Idx]] <- list()

        
        
        FdrStr[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
        FdrStrPrn[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
        FdrStrRes[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp

        
        # FDR values' stadistics.
        print(paste("Computing FDR values for displacement",(GlbEst[[Idx]]$Dsp)))
        ListFDR = list(DspData[[5]],DspData[[6]],DspData[[7]],DspData[[8]])
        FdrStr[[Idx]]$Data = ListFDR
        FdrStrPrn[[Idx]]$Data = list(DspData[[5]],DspData[[6]])
        FdrStrRes[[Idx]]$Data = list(DspData[[7]],DspData[[8]])
        VecFDR = unlist(ListFDR)
        #print(VecFDR)
        MeaFDR = mean(VecFDR)
        VecMeaFDR[Idx] = MeaFDR
        if (MeaFDR<ThrFDR) MinDsp = GlbEst[[Idx]]$Dsp
        ThrFDR = min(ThrFDR,MeaFDR)
      }
      
      # Binding values' stadistics.
      print(paste("Computing binding scores for displacement",(GlbEst[[Idx]]$Dsp)))
      
      if ( ( is.null(info_to_save ) ) || ( "bnd" %in% info_to_save ) ) {
        
        BndStr[[Idx]] <- list()
        BndStrPrn[[Idx]] <- list()
        BndStrRes[[Idx]] <- list()
        
        BndStr[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
        BndStrPrn[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
        BndStrRes[[Idx]]$Dsp <- GlbEst[[Idx]]$Dsp
        
        BndForProne = DspData[[9]]
        BndRevProne = DspData[[10]]
        BndForResis = DspData[[11]]
        BndRevResis = DspData[[12]]
        
        registerDoMC(Config$nCPU)
        
        LstBndForPrn=list()
        if (length(BndForProne)!=0){
          
          LstBndForPrn = foreach (w=1:length(BndForProne), .combine=c) %dopar% {
            BndForProneVec = BndForProne[[w]]
            
            PList=list()
            j = 1
            if (length(BndForProneVec)!=0){
              for (i in 1:length(BndForProneVec)){
                PList[j] = abs(BndForProneVec[[j]]$rating.mean[1] - BndForProneVec[[j]]$rating.mean[2])
                j = j + 1
              }
            }
            PList
          }
        }
        
        LstBndRevPrn=list()
        if (length(BndRevProne)!=0){
          
          LstBndRevPrn = foreach (w=1:length(BndRevProne), .combine=c) %dopar% {
            BndRevProneVec = BndRevProne[[w]]
            
            PList=list()
            j = 1
            if (length(BndRevProneVec)!=0){
              for (i in 1:length(BndRevProneVec)){
                PList[j] = abs(BndRevProneVec[[j]]$rating.mean[1] - BndRevProneVec[[j]]$rating.mean[2])
                j = j + 1
              }
            }
            PList
          }
        }
        
        LstBndForRes=list()
        if (length(BndForResis)!=0){
          
          LstBndForRes = foreach (w=1:length(BndForResis), .combine=c) %dopar% {
            BndForResisVec = BndForResis[[w]]
            
            PList=list()
            j = 1
            if (length(BndForResisVec)!=0){
              for (i in 1:length(BndForResisVec)){
                PList[j] = abs(BndForResisVec[[j]]$rating.mean[1] - BndForResisVec[[j]]$rating.mean[2])
                j = j + 1
              }
            }
            PList
          }
        }
        
        LstBndRevRes=list()
        if (length(BndRevResis)!=0){
          
          LstBndRevRes = foreach (w=1:length(BndRevResis), .combine=c) %dopar% {
            BndRevResisVec = BndRevResis[[w]]
            
            PList=list()
            j = 1
            if (length(BndRevResisVec)!=0){
              for (i in 1:length(BndRevResisVec)){
                PList[j] = abs(BndRevResisVec[[j]]$rating.mean[1] - BndRevResisVec[[j]]$rating.mean[2])
                j = j + 1
              }
            }
            PList
          }
        }
        
        # HERE
        BndStr[[Idx]]$Data = list()
        BndStr[[Idx]]$Data[[1]] = c(LstBndForPrn, LstBndForRes)
        BndStr[[Idx]]$Data[[2]] = c(LstBndRevPrn, LstBndRevRes)
        
        BndStrPrn[[Idx]]$Data = list()
        BndStrPrn[[Idx]]$Data[[1]] = LstBndForPrn
        BndStrPrn[[Idx]]$Data[[2]] = LstBndRevPrn
        
        BndStrRes[[Idx]]$Data = list()
        BndStrRes[[Idx]]$Data[[1]] = LstBndForRes
        BndStrRes[[Idx]]$Data[[2]] = LstBndRevRes
      }
      
      
    }
    
    # Count number of POMs for each displacement.
    
    if ( ( is.null(info_to_save ) ) || ( "fdr" %in% info_to_save ) ) {
      save(FdrStr, FdrStrPrn, FdrStrRes, file=paste("MergeStructuresFDR",Sys.Date() ,".RData", sep=""))
      print(paste("Minimum mean FDR:", ThrFDR, "for displacement", MinDsp))
    }
    if ( ( is.null(info_to_save ) ) || ( "cnt" %in% info_to_save ) ) {
      save(CntStr, CntStrPrn, CntStrRes, file=paste("MergeStructuresCNT",Sys.Date() ,".RData", sep=""))
    }
    if ( ( is.null(info_to_save ) ) || ( "bnd" %in% info_to_save ) ) {
      save(BndStr, BndStrPrn, BndStrRes, file=paste("MergeStructuresBND",Sys.Date() ,".RData", sep=""))
    }
    
    # DR(FdrStr)
    # HSDispBnd(BndStr)
    
    
    #return(VecMeaFDR)
    #return(BndStr)
    return(TRUE)
  }
  else{
    print("Could not find data for any displacement.")
    return(NULL)
  }
  
}

allTests=function(BndStr, BndStrPrn, BndStrRes,
                  CntStr, CntStrPrn, CntStrRes,
                  FdrStr, FdrStrPrn, FdrStrRes){
  # Do this at first execution?
  font_import()
  print("Calculating fdr graphs")
  HMDispDFDR(FdrStr, "pr")
  print("Calculating prone fdr graphs")
  HMDispDFDR(FdrStrPrn, "p")
  print("Calculating resistant fdr graphs")
  HMDispDFDR(FdrStrRes, "r")
  print("Calculating fdr graphs (log)")
  HMFDRDispLog(FdrStr, "h", "pr")
  print("Calculating prone fdr graphs (log)")
  HMFDRDispLog(FdrStrPrn, "h", "p")
  print("Calculating resistant fdr graphs (log)")
  HMFDRDispLog(FdrStrRes, "h", "r")
  HMFDRDispLog(FdrStr, "l", "pr")
  HMFDRDispLog(FdrStrPrn, "l", "p")
  HMFDRDispLog(FdrStrRes, "l", "r")
  HMFDRDispLog(FdrStr, "hl", "pr")
  HMFDRDispLog(FdrStrPrn, "hl", "p")
  HMFDRDispLog(FdrStrRes, "hl", "r")
  print("Calculating binding score graphs")
  HSDispBnd(BndStr, CntStr, "pr", nPomLimUp = 5500)
  print("Calculating binding score prone graphs")
  HSDispBnd(BndStrPrn, CntStrPrn, "p")
  print("Calculating binding score resistant graphs")
  HSDispBnd(BndStrRes, CntStrRes, "r")
  HMDispBnd(BndStr, "pr")
  HMDispBnd(BndStrPrn, "p")
  HMDispBnd(BndStrRes, "r")
  return(TRUE)
}

allTests=function(BndStr, BndStrPrn, BndStrRes,
                  CntStr, CntStrPrn, CntStrRes,
                  FdrStr, FdrStrPrn, FdrStrRes){
  # Do this at first execution?
  #font_import()
  print("Calculating fdr graphs")
  HMDispDFDR(FdrStr, "pr")
  print("Calculating prone fdr graphs")
  HMDispDFDR(FdrStrPrn, "p")
  print("Calculating resistant fdr graphs")
  HMDispDFDR(FdrStrRes, "r")
  print("Calculating fdr graphs (log)")
  HMFDRDispLog(FdrStr, "h", "pr")
  print("Calculating prone fdr graphs (log)")
  HMFDRDispLog(FdrStrPrn, "h", "p")
  print("Calculating resistant fdr graphs (log)")
  HMFDRDispLog(FdrStrRes, "h", "r")
  HMFDRDispLog(FdrStr, "l", "pr")
  HMFDRDispLog(FdrStrPrn, "l", "p")
  HMFDRDispLog(FdrStrRes, "l", "r")
  HMFDRDispLog(FdrStr, "hl", "pr")
  HMFDRDispLog(FdrStrPrn, "hl", "p")
  HMFDRDispLog(FdrStrRes, "hl", "r")
  print("Calculating binding score graphs")
  HSDispBnd(BndStr, CntStr, "pr", nPomLimUp = 5500)
  print("Calculating binding score prone graphs")
  HSDispBnd(BndStrPrn, CntStrPrn, "p")
  print("Calculating binding score resistant graphs")
  HSDispBnd(BndStrRes, CntStrRes, "r")
  HMDispBnd(BndStr, "pr")
  HMDispBnd(BndStrPrn, "p")
  HMDispBnd(BndStrRes, "r")
  return(TRUE)
}

testsTFG=function(FdrStrLoadPath, BndStrLoadPath, CntStrLoadPath, cell_name, pomlimuProne, pomlimupRes, maxDifProne, maxDifRes){
  
  load(FdrStrLoadPath)
  load(BndStrLoadPath)
  load(CntStrLoadPath)
  
  HMDispDFDR(FdrStrPrn, "p", cell_name)
  HMDispDFDR(FdrStrRes, "r", cell_name)
  HSDispBnd(BndStrPrn, CntStrPrn, "p", cell_name, nPomLimUp = pomlimuProne, max_dif = maxDifProne)
  HSDispBnd(BndStrRes, CntStrRes, "r", cell_name, nPomLimUp = pomlimupRes ,max_dif=maxDifRes)
}
