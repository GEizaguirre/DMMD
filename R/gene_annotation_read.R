filtreadmdfile=function(Config){
  
  data <- NULL
  # Sometimes the scan fails
  attempt <- 1
  while( is.null(data) && attempt <= 4 ) {
    attempt <- attempt + 1
    try(
      {
        # Gets the desired annotations (not of type "Un")
        # from the gene annotations' file.
        dirmd <- paste(Config$DirFas, "/seq_gene_filt.txt", sep="")
        # Create an R-manageable file with the desired annotations.
        val <- .Call("filtmdfile", Config$NumAutosomes, Config$DirFas)
        # Create an R structure with the annotations.
        data <- read.table(dirmd)
        colnames(data) <- c("Chr","PosIni","PosEnd","GenNam","GenID") 
      }
    )
  }
  
  if (is.null(data)){
    stop(paste("filtreadmdfile: Could not read gene annotation file", paste0(Config$DirFas, "/seq_gene.md")))
  }
  
  return(data)
}
