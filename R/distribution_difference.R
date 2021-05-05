CheckDifferencesInDistributions=function(Config,SelMoti,LenMotif,TypeMotif){
  
  i=0
  # Specify alternative type for kolmogorov test.
  AltHyp=switch(TypeMotif, 
                ForProne = "less", 
                RevProne = "less",
                ForResis = "greater",
                RevResis = "greater")
  
  N = 10000
  # Returning structure.
  FinMot=list()
  
  gc()
  # Parallelization set up.
  registerDoMC(Config$nCPU)
  
  
  for(w in LenMotif){
    
    #Log
    if (!is.null(SelMoti[[w]]) && (length(SelMoti[[w]])>=5)) {
      
      # Extract pwm, POMs, prone sequences' scores and resistant sequences' scores.
      Pwm = SelMoti[[w]][[1]]
      Ind = SelMoti[[w]][[2]]
      ScoPro = SelMoti[[w]][[3]]
      ScoRes = SelMoti[[w]][[4]]
      # Added:
      FDRs = SelMoti[[w]][[5]]
      
      # Select motifs to be removed.
      IdxRm=foreach(i=1:length(ScoPro),.combine=c) %dopar% {
        
        ScoDatFraPro=ScoPro[[i]]
        ScoDatFraRes=ScoRes[[i]]
        ScoProVec=unlist(sapply(1:length(ScoDatFraPro[,1]), function(i) rep(ScoDatFraPro[,2][i],ScoDatFraPro[,1][i])))
        ScoResVec=unlist(sapply(1:length(ScoDatFraRes[,1]), function(i) rep(ScoDatFraRes[,2][i],ScoDatFraRes[,1][i])))
        
        # score volume reduction for faster execution and supression of ties.
        sample1 <- sample(ScoProVec, N, replace = FALSE)
        sample2 <- sample(ScoResVec, N, replace = FALSE)
        
        # Kolmogorov test or Mann-Whitney test for binding scores.
        if (Config$DistDif=="ks"){
          # out=ks.test(ScoProVec,ScoResVec,alternative=AltHyp)
          out = ks.test(sample1, sample2)
        }else{
          out = wilcox.test(sample1, sample2,alternative="two.sided",paired=FALSE)
        }
        
        if(out$p.value > Config$significance_level) i 
      }
      if (! is.null(IdxRm)){
        Pwm=Pwm[-IdxRm]
        Ind=Ind[-IdxRm]
        ScoPro=ScoPro[-IdxRm]
        ScoRes=ScoRes[-IdxRm]
        # Added:
        FDRs=FDRs[-IdxRm]
      }
      #FinMot[[w]]=list(Pwm,Ind,ScoPro,ScoRes)
      # Added:
      FinMot[[w]]=list(Pwm,Ind,ScoPro,ScoRes,FDRs)
    }
  }
  return(FinMot)
}