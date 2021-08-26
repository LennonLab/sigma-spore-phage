#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(seqinr)
library(treeio)



# Get duplicate list -------------------------------------------------------

log <- readLines(here("phylo/data/reduced_set_to_align/check_msa/check-msa.raxml.log"))
log <- log[str_detect(log, "identical")]
dups <- str_extract_all(log,"(YP_[0-9]*..-(phage|bacteria))", simplify = T) %>% as_tibble()


# Get kept sequences ------------------------------------------------------
kept <- read.phylip.seq(here("phylo/data/reduced_set_to_align/check_msa/check-msa.raxml.reduced.phy"))
kept <- names(kept)

dups$V1 %in% kept #all TRUE
dups$V2 %in% kept # all FALSE


# identify removed duplicates ---------------------------------------------


# load viral sigmas data from vog HMM analysis
load(here("vogdb/data/vog_sigma_clean_Whost.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))

d.removed <- 
  d.faa %>% 
  mutate(id = paste0(protein,"-phage")) %>% 
  filter(id %in% dups$V2) %>% 
  left_join(., dups, by = c("id" = "V2")) %>% 
  rename(dup.of = V1)
  

# Save results ------------------------------------------------------------


# select(d.faa, sp2 = sp, description2 = description, protein) %>% 
#   left_join(d.phage, . ,by = c("V2" = "protein"))
