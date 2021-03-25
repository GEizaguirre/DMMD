#######Fusion

FuseSeq = function(Config, SeqW,MetW,FreW,IndW,FreWVec,SeqNextW,MetNextW,FreNextW,IndNextW,FreNextWVec,IndFou,IndFoUpd){
  
  for(i1 in IndFoUpd){
    
    frew=FreW[i1]
    metw=MetW[i1]
    
    Len=length(IndFou[[i1]])
    
    for(k in 1:Len){
      
      freNextw=FreNextW[IndFou[[i1]][k]]
      metNextw=MetNextW[IndFou[[i1]][k]]
      
      MetUpd=(frew*metw+freNextw*metNextw)/(frew+freNextw)
      
      MetNextW[IndFou[[i1]][k]]=MetUpd
      
      FreWVecNew=FreWVec[i1,]
      
      if (Config$GrowingMode=="C")  FreWVecNew=c(0,FreWVecNew,0)
      if (Config$GrowingMode=="R")  FreWVecNew=c(FreWVecNew,0,0)
      if (Config$GrowingMode=="L")  FreWVecNew=c(0,0,FreWVecNew)
      
      FreNextWVec[IndFou[[i1]][k],]=FreNextWVec[IndFou[[i1]][k],]+FreWVecNew
      
      FreNextW[IndFou[[i1]][k]]=FreNextW[IndFou[[i1]][k]]+frew
    }
  }
  
  Fus=list(SeqW,MetW,FreW,IndW,FreWVec,SeqNextW,MetNextW,FreNextW,IndNextW,FreNextWVec)
  return(Fus)
}

Fusion=function(Config,SeqMetFreW,FreVecW){
  
  SeqMetFreWFreVecW=list()
  
  # For elements not in the vector.
  '%ni%' = Negate('%in%') 
  
  for(w in Config$w_min:(Config$w_max-1)) {
    
    if ( !is.null(SeqMetFreW[[w]]) && nrow(SeqMetFreW[[w]]) != 0 ){
      
      # Get sequences' column.
      SeqW=SeqMetFreW[[w]][,1]
      # Get 
      MetW=SeqMetFreW[[w]][,2]
      FreW=SeqMetFreW[[w]][,3]
      IndW=SeqMetFreW[[w]][,4]
      FreWVec=FreVecW[[w]]
      
      SeqNextW=SeqMetFreW[[w+1]][,1]
      MetNextW=SeqMetFreW[[w+1]][,2]
      FreNextW=SeqMetFreW[[w+1]][,3]
      IndNextW=SeqMetFreW[[w+1]][,4]
      FreNextWVec=FreVecW[[w+1]]
      
      
      if (Config$GrowingMode=="C") {
        cut_point <- 2*w+1
        SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 2, cut_point)))
      }
      if (Config$GrowingMode=="R") SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 1, 2*(w+1))))
      if (Config$GrowingMode=="L") SeqNextWCut <- unlist(lapply(SeqNextW, function(x) substring(x, 3, 2*(w+1)+2)))
      
      mode(FreW) <- "integer"
      mode(MetW) <- "numeric"
      mode(FreWVec) <- "integer"
      mode(FreNextW) <- "integer"
      mode(MetNextW) <- "numeric"
      mode(FreNextWVec) <- "integer"
      
      # Get equal sequences in Cpp
      # IndFou <- cpp_str_sort(SeqW, SeqNextWCut)
      ind_data <- find_strings_par(SeqW, SeqNextWCut, Config$nCPU + 1)
      
      # Get index of elements that have to be fusioned.
      # IndFoUpd=which(IndFou %ni% 0)
      
      # OutFus=FuseSeq(Config, SeqW,MetW,FreW,IndW,FreWVec,SeqNextW,MetNextW,FreNextW,IndNextW,FreNextWVec,IndFou,IndFoUpd)
      grow_mode_numeric <- switch(Config$GrowingMode,
                                   C=0,
                                   L=1,
                                   R=2)
      
      # ret <- fuse_seqs_c(2*w, grow_mode_numeric, IndFou, IndFoUpd, 
      #                    FreW, MetW, FreWVec, 
      #                    FreNextW, MetNextW, FreNextWVec)
      IndFu <- ind_data$str_indexes
      IndFu_sub <- ind_data$subindexes
      IndFu_not_void <- ind_data$ind_not_void
      seqs_to_rm <- ind_data$seqs_to_rm
      
      fuse_seqs_openmp(2*w, grow_mode,IndFu, IndFu_sub, IndFu_not_void,
                       FreW, MetW, FreWVec,
                       FreNextW, MetNextW, FreNextWVec,
                       Config$nCPU)
      
      OutFus=list(SeqW,MetW,FreW,IndW,FreWVec,SeqNextW,MetNextW,FreNextW,IndNextW,FreNextWVec)
      
      SeqW=OutFus[[1]]
      MetW=OutFus[[2]]
      FreW=OutFus[[3]]
      IndW=OutFus[[4]]
      FreWVec=OutFus[[5]]

      SeqNextW=OutFus[[6]]
      MetNextW=OutFus[[7]]
      FreNextW=OutFus[[8]]
      IndNextW=OutFus[[9]]
      FreNextWVec=OutFus[[10]]
      
      if(length(seqs_to_rm!=0)){
        SeqW=SeqW[-seqs_to_rm]
        MetW=MetW[-seqs_to_rm]
        FreW=FreW[-seqs_to_rm]
        IndW=IndW[-seqs_to_rm]
        FreWVec=FreWVec[-seqs_to_rm,]
      }
      
      # Delete posible void words.
      IndVoiOut = which(SeqW=="")
      if(length(IndVoiOut!=0)){
        SeqW=SeqW[-IndVoiOut]
        MetW=MetW[-IndVoiOut]
        FreW=FreW[-IndVoiOut]
        IndW=IndW[-IndVoiOut]
        FreWVec=FreWVec[-IndVoiOut,]
      }
      
      SeqMetFreW[[w]]=data.frame(SeqW,MetW,FreW,IndW)
      FreVecW[[w]]=FreWVec
      
      SeqMetFreW[[w+1]]=data.frame(SeqNextW,MetNextW,FreNextW,IndNextW)
      FreVecW[[w+1]]=FreNextWVec
      
      #####Keep the data
      
      if (length(SeqW)==1){
        mat <- matrix(vector(), nrow = 1, ncol = 2*w+2)
        for (i in 1:(2*w+2))
          mat[1,i]=as.numeric(FreWVec[i])
        SeqMetFreWFreVecW[[w]]=data.frame(SeqW,MetW,FreW,IndW,mat)
      }
      else
      {
        SeqMetFreWFreVecW[[w]]=data.frame(SeqW,MetW,FreW,IndW,FreWVec)
      }
      
      if(w==(Config$w_max-1)){
        
        SeqMetFreWFreVecW[[Config$w_max]]=data.frame(SeqNextW,MetNextW,FreNextW,IndNextW,FreNextWVec)
        colnames(SeqMetFreWFreVecW[[Config$w_max]])[1]="SeqW"
        colnames(SeqMetFreWFreVecW[[Config$w_max]])[2]="MetW"
        colnames(SeqMetFreWFreVecW[[Config$w_max]])[3]="FreW"
        colnames(SeqMetFreWFreVecW[[Config$w_max]])[4]="IndW"
        
      }
    }
    
  }
  
  return(SeqMetFreWFreVecW)
}