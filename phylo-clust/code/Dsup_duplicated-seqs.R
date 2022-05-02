#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(seqinr)
library(treeio)



# Get duplicate list -------------------------------------------------------

log <- readLines(here("phylo-clust/data/align-trim-tree/check_msa/check-msa.raxml.log"))
# log <- readLines(here("phylo/data/align-trim-tree/model_test/model_test.log"))
log <- log[str_detect(log, "identical")]
dups <- str_extract_all(log,"(.P_[0-9]*..-(phage|bacteria))", simplify = T) %>% 
  as_tibble(column_name = c("V1", "V2"))


# Get kept sequences ------------------------------------------------------
kept <- read.phylip.seq(here("phylo-clust/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"))
kept <- names(kept)

dups$V1 %in% kept #all TRUE
dups$V2 %in% kept # all FALSE


# identify removed duplicates ---------------------------------------------


# load viral sigmas data from vog HMM analysis
load(here("phylo-clust/data/faa_table_clustered.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))

d.removed <- 
  d.new %>% 
  mutate(id = paste0(protein,"-phage")) %>% 
  filter(id %in% dups$V2) %>% 
  left_join(., dups, by = c("id" = "V2")) %>% 
  rename(dup.of = V1)

# d.removed$sp
# "Pectinobacterium phage vB_PcaM_CBB" "Bacillus phage Megatron" 

# no need to do anything


# old code ------------
  
# # swapping to have Bcp1 protein in MSA
# #two of the Bcp1 trimmed proteins are identical to proteins from Bacillus virus BM15
# pid.out <- d.removed %>% filter(str_detect(sp, "Bcp1")) %>% pull(dup.of)
# pid.in <-  d.removed %>% filter(str_detect(sp, "Bcp1")) %>% pull(id)
# 
# 
# # Save results ------------------------------------------------------------
# # changing the MSA directly
# phyl.in <- readLines(here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"))
# phyl.out <- gsub( pid.out[1], pid.in[1], phyl.in )
# phyl.out <- gsub( pid.out[2], pid.in[2], phyl.in )
# cat(phyl.out,
#     file=here("phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy"),
#     sep="\n")

