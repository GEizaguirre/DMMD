library(DMMD)

print("loading data")
list.files("..")
load("../fusion_data_0.Rdata")

w = 6

Config <- list()
Config$GrowingMode <- "C"
Config$nCPU <- 3

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

# SeqW <- SeqW[1:100000]
# SeqNextWCut <- SeqNextWCut[1:100000]

# Get equal sequences in Cpp
print("finding strings")
start_time <- Sys.time()
ind_data <- find_strings_par(SeqW, SeqNextWCut, Config$nCPU + 3)
end_time <- Sys.time()

print(end_time-start_time)

print("saving image")
save.image(file="fusion_omp_test.RData")

# grow_mode_numeric <- switch(Config$GrowingMode,
#                             C=0,
#                             L=1,
#                             R=2)
# 
# IndFu <- ind_data$str_indexes
# IndFu_sub <- ind_data$subindexes
# IndFu_not_void <- ind_data$ind_not_void
# seqs_to_rm <- ind_data$seqs_to_rm

# print("fusing seqs")
# fuse_seqs_openmp(2*w, grow_mode_numeric, IndFu, IndFu_sub, IndFu_not_void,
#                  FreW, MetW, FreWVec,
#                  FreNextW, MetNextW, FreNextWVec,
#                  Config$nCPU)
# 
# print("Finished succesfully")
