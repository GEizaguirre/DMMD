nelem1 <- 10
nelem2 <- 10
l <- 2
Config <- list()
Config$lambda <- l
type_motif <- "RevResis"
num_motif <- 20
type_motif_numeric <- switch(type_motif,
                             ForProne=0,
                             RevProne=1,
                             ForResis=2,
                             RevResis=3)

ScanPro <- list()
ScanRes <- list()
for (i in 1:num_motif) {
  pro_reps <- sample.int(20, nelem1, replace = TRUE)
  pro_vals <- runif(nelem1)
  res_reps <- sample.int(20, nelem2, replace = TRUE)
  res_vals <- runif(nelem2)
  ScanPro[[i]] <- cbind(pro_reps, pro_vals)
  ScanRes[[i]] <- cbind(res_reps, res_vals)
}

ScanPro[[num_motif+1]] <- list()
ScanRes[[num_motif+1]] <- list()


# Parallelization set up.
registerDoMC(Config$nCPU)
lambda <- Config$lambda


# For each motif.
VecFdr <- foreach(i=1:(length(ScanRes)), .combine=c) %dopar% {
  
  # Prone sequences' binding counts and breaks for this motif.
  ScoDatFraPro=ScanPro[[i]]
  # Resistant sequences' binding counts and breaks for this motif.
  ScoDatFraRes=ScanRes[[i]]
  
  if ( length(ScoDatFraPro) == 0 || length(ScoDatFraRes) == 0 )
    fdr <- 1.0
  else
    # Get fdr value in cpp
    fdr <- fdr_c(ScoDatFraPro, ScoDatFraRes, lambda, type_motif_numeric)
  fdr
}
