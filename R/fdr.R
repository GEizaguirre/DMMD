###1st Filter: FDR

FDR <- function(Config,ScoMotPro,ScoMotRes,PWMFor,POMClu,LenMotif,TypeMotif){
  
  #
  # Config: configuration parameters' structure.
  # ScoMotPro: scores for methylation prone sequences.
  # Structure: [[w]] lengths
  #                     [[Motif1]] -> [bincounts, binbreaks]
  #                     [[Motif2]] -> [bincounts, binbreaks]
  #                     ...
  #             ...
  # ScoMotRes: scores for methylation resistant sequences.
  # Structure: [[w]] lengths
  #                     [[Motif1]] -> [bincounts, binbreaks]
  #                     [[Motif2]] -> [bincounts, binbreaks]
  #                     ...
  #             ...
  # PWMFor: normalization Position Weight Matrices for motifs.
  # POMClu: Position Occurence Matrices for motifs.
  # LenMotif: lengths with annotated motifs.
  # TypeMotif: Forward|Reverse and Prone|Resistant.
  
  i=0
  # Resulting structure.
  SelMoti=list()
  FdrPerLen=list()
  
  
  for(w in LenMotif){
    
    if (!is.null(ScoMotPro[[w]]) && !is.null(ScoMotRes[[w]])){
      
      # Scores of methylation prone sequences for this length.
      ScanPro <- ScoMotPro[[w]] 
      # Scores of methylation resistant sequences for this length.
      ScanRes <- ScoMotRes[[w]]
      type_motif_numeric <- switch(TypeMotif,
                                   ForProne=0,
                                   RevProne=1,
                                   ForResis=2,
                                   RevResis=3)
      lambda <- Config$lambda
    
      # Parallelization set up.
      registerDoMC(Config$nCPU)
      
      # For each motif.
      VecFdr <- foreach(i=1:(length(ScanRes)),.combine=c) %dopar% {
        
        fdr <-tryCatch (
        {
            # Prone sequences' binding counts and breaks for this motif.
            ScoDatFraPro <- ScanPro[[i]]
            # Resistant sequences' binding counts and breaks for this motif.
            ScoDatFraRes <- ScanRes[[i]]
            
            pro_counts <- ScoDatFraPro$bincounts
            mode(pro_counts) <- "integer"
            pro_breaks <- ScoDatFraPro$binbreaks
            mode(pro_breaks) <- "numeric"
            res_counts <- ScoDatFraRes$bincounts
            mode(res_counts) <- "integer"
            res_breaks <- ScoDatFraRes$binbreaks
            mode(res_breaks) <- "numeric"
            
            if ( length(pro_counts) == 0 || length(pro_breaks) == 0 )
              fdr_ <- 1.0
            else
              # Get fdr value in cpp
              fdr_ <- fdr_c(pro_counts, pro_breaks,
                           res_counts, res_breaks,
                           lambda, type_motif_numeric)
            fdr_
          },
          error = function(e){
            1.0
          }
        )
        
        fdr
      }
      
      # We have the fdr for each motif.
      # Select those motifs with FDR value above the threshold.
      IdxRm=which(VecFdr>Config$fdr)
      
      # Remove non-desirable motifs for this length.
      if (! is.null(IdxRm)){
        PWMForUpd=PWMFor[[w]][-IdxRm]
        POMCluUpd=POMClu[[w]][[2]][-IdxRm]
        ScoMotProUpd=ScoMotPro[[w]][-IdxRm]
        ScoMotResUpd=ScoMotRes[[w]][-IdxRm]
        # Added:
        VecFDRUpd=VecFdr[-IdxRm]
      }
      #SelMoti[[w]]=list(PWMForUpd,POMCluUpd,ScoMotProUpd,ScoMotResUpd,VecFdr)
      SelMoti[[w]]=list(PWMForUpd,POMCluUpd,ScoMotProUpd,ScoMotResUpd,VecFDRUpd)
    }
  }
  
  #Log
  # line <- paste("fdr completed for ", LenMotif)
  # write(line, file=Config$LogFile, append=TRUE)
  
  # Return filtered motifs.
  return(SelMoti)  
}