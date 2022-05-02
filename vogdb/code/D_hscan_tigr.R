library(here)
library(tidyverse)
library(foreach)


# this script relies on the windows subsytem for Linux (wsl)
# on which hmmer has been installed
setwd(here())

# Bacterial sigma HMM data base
hmm <- "TIGR/data/tigr_sigma_hmm/sigma.fam"

faa <- dir("vogdb/data/vog_sigma_clean", full.names = T, pattern = ".faa")

# create output directory
if (! dir.exists(here("vogdb/data/hscan_vogXtigr/"))){
  dir.create("vogdb/data/hscan_vogXtigr/")
}


foreach (f = faa) %do%{
    
    faa.name <- 
      str_remove(f,".faa.gz")%>%
      str_remove(".faa")%>%
      str_remove(".*/")
    
    
    filename <- paste0("vogdb/data/hscan_vogXtigr/",faa.name,".txt")
    
    #no threshold
    wsl <- paste("wsl hmmscan --noali --tblout", filename, hmm, f)
    shell(wsl)
  
}

# for (f in faa){
#   
#   faa.name <- 
#     str_remove(f,".faa.gz")%>%
#     str_remove(".*/")
#   
#   
#   filename <- paste0("vogdb/data/hscan_vogXtigr/",faa.name,".txt")
#   
#   #no threshold
#   wsl <- paste("wsl hmmscan --noali --tblout", filename, hmm, f)
#   shell(wsl)
#   
# }
# 
