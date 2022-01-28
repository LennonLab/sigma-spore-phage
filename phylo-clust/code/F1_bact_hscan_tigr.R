library(here)
library(tidyverse)
library(seqinr)
library(foreach)

# classify bacterial sigmas in same way I did for phage sigmas

# this script relies on the windows subsytem for Linux (wsl)
# on which hmmer has been installed

setwd(here())

hmm <- "TIGR/data/tigr_sigma_hmm/sigma.fam"

faa <- here("phylo-clust/data/bacterial_sigma.faa")
# make indivdual files for each bacterial sigma in tmp folder
if (! dir.exists(here("phylo-clust/data","tmp_sigma"))){
  dir.create(here("phylo-clust/data","tmp_sigma"))
}
l.faa <- read.fasta(faa)
for (i in 1:length(l.faa)){
  file.out <- getAnnot(l.faa[[i]]) %>% 
    str_remove(">") %>% str_remove(".. .*")%>% 
    paste0(.,".faa")
  
  write.fasta(sequences = getSequence(l.faa[[i]]),
              names = getAnnot(l.faa[[i]]) %>% str_remove(">"),
              file.out = here("phylo-clust/data","tmp_sigma",file.out))
}

faa <- list.files(here("phylo-clust/data","tmp_sigma"), full.names = F) %>% 
  paste0("phylo-clust/data/tmp_sigma/",.)

#make directory for results
if (! dir.exists(here("phylo-clust/data","hscan_bactXtigr"))){
  dir.create(here("phylo-clust/data","hscan_bactXtigr"))
}

#  scan
  foreach(f = faa, .packages = "stringr") %do%
{
    
    faa.name <- 
      str_remove(f,".faa")%>%
      str_remove(".*/")
    
    # hmm.name <- 
    #   str_remove(h,".HMM")%>%
    #   str_remove(".*/")
    
    filename <- paste0("phylo-clust/data/hscan_bactXtigr/",faa.name,".txt")
    
    # wsl <- paste("wsl hmmsearch --noali -T 20 --tblout", filename, h, f)
    # wsl <- paste("wsl hmmscan --noali -T 20 --tblout", filename, h, f)
    #no threshold
    wsl <- paste("wsl hmmscan --noali --tblout", filename, hmm, f)
    shell(wsl)
}

#clean up
unlink(here("phylo-clust/data","tmp_sigma"), recursive = T)
