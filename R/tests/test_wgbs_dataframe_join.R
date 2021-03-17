library(dplyr)
library(data.table)
library(comprehenr)
source("R/wgbs_data_read_functions.R")

config <- list()

config$NumAutosomes = 22
config$NumChromosomes = 25
config$Allosomes = c("X", "Y", "M")
config$DirDat = "/mnt/1814c93f-609b-4281-819a-1f3cfab40622/Biodonostia/Methylation_datasets/ENCODE_WGBS"


files <- list.files(config$DirDat)

if (length(files) == 0) stop("WGBS file not in directory")

# Get methylation rates in the directory
if (length(files) > 1){
  
  CooMet <- list()
  
  # As there are more than one file, we perform an outer join between each of them and calculate the average
  # methylation rate for each common methylation coordinate.
  CooMetFor <- list()
  CooMetRev <- list()
  CooMets <- to_list(for (i in length(files)) get_methylation_coordinates_bedMethyl(config, config$DirDat, files[i]))
  # Calculate methylation average for each coordinate.
  for (ci in config$NumChromosomes){
    
    dfs = to_list(for (cm in CooMets) cm[[1]][[ci]])
    CooMetFor[[ci]] <- rbindlist(dfs)[,.(ColMet=mean(ColMet)),ColCoo]
    
    dfs = to_list(for (cm in CooMets) cm[[2]][[ci]])
    CooMetRev[[ci]] <- rbindlist(dfs)[,.(ColMet=mean(ColMet)),ColCoo]
    
  }
  CooMet[[1]] <- CooMetFor
  CooMet[[2]] <- CooMetRev
  
} else {
  CooMet <- get_methylation_coordinates_bedMethyl(config, config$DirDat, files[1])
}