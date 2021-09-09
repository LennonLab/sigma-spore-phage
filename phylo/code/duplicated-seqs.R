#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(seqinr)
library(treeio)



# Get duplicate list -------------------------------------------------------

log <- readLines(here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.log"))
# log <- readLines(here("phylo/data/align-trim-tree/model_test/model_test.log"))
log <- log[str_detect(log, "identical")]
dups <- str_extract_all(log,"(.P_[0-9]*..-(phage|bacteria))", simplify = T) %>% as_tibble()


# Get kept sequences ------------------------------------------------------
kept <- read.phylip.seq(here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"))
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
  
# swapping to have Bcp1 protein in MSA
#Bcp1 trimmed protein is identical to a protein from Bacillus virus BM15
pid.out <- d.removed %>% filter(str_detect(sp, "Bcp1")) %>% pull(dup.of)
pid.in <-  d.removed %>% filter(str_detect(sp, "Bcp1")) %>% pull(id)


# Save results ------------------------------------------------------------
# changing the MSA directly
phyl.in <- readLines(here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"))
phyl.out <- gsub( pid.out, pid.in, phyl.in )
cat(phyl.out,
    file=here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"),
    sep="\n")

