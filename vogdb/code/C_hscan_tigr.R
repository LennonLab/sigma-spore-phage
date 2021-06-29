library(here)
library(tidyverse)


# this script relies on the windows subsytem for Linux (wsl)
# on which hmmer has been installed
setwd(here())

hmm <- "TIGR/tigr_sigma_hmm/sigma.fam"

faa <- dir("vogdb/data/vog_sigma_clean", full.names = T, pattern = ".faa.gz")


for (h in hmm){
  for (f in faa){
    
    faa.name <- 
      str_remove(f,".faa.gz")%>%
      str_remove(".*/")
    
    # hmm.name <- 
    #   str_remove(h,".HMM")%>%
    #   str_remove(".*/")
    
    filename <- paste0("vogdb/data/hscan_vogXtigr/",faa.name,".txt")
    
    # wsl <- paste("wsl hmmsearch --noali -T 20 --tblout", filename, h, f)
    # wsl <- paste("wsl hmmscan --noali -T 20 --tblout", filename, h, f)
    #no threshold
    wsl <- paste("wsl hmmscan --noali --tblout", filename, h, f)
    shell(wsl)
  }
}


