#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)
library(ggtreeExtra) #https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
library(ggnewscale)


# load data ---------------------------------------------------------------

load(here("phylo-clust/data/plotted_tree.Rdata"))
load(here("phylo-clust/data/faa_table_clustered.RData"))

# load viral sigmas data from vog HMM analysis
tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")
load(here("vogdb/data/vog_sigma_clean_tigr.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
d.phage.tigr <- d.faa %>%
  filter(protein %in% d.new$protein) %>% 
  select(protein, tigr.hit) %>% 
  mutate(sig.class = case_when(
    tigr.hit %in% tigr.spore ~ "sporulation",
    tigr.hit == "no_hit" ~ "no_hit",
    TRUE ~ "other"))
rm(d.faa)

# add tigr calssification to main dataframe
d.new <- 
  left_join(d.new, d.phage.tigr)


# identify MRCAs -----------------------------------------------------------
mrca.sigSPORE <- d.rax %>% 
  filter(bs.label %in% c("sigF", "sigG", "sigE")) %>%
  pull(node) %>% 
  MRCA(d.rax, .) %>% pull(node)

mrca.sigB <- d.rax %>% 
  filter(tip.label %in% c("Bacillus_subtilis_sigB", 
                          "Mycobacterium_tuberculosis_sigF",
                          "Clostridium phage phi3626")) %>%
  pull(node) %>% 
  MRCA(d.rax, .) %>% pull(node)


# Mark sporulation sigmas from phylogeny -------------------------------------------
all.spore <- d.rax %>% 
   offspring(.node = mrca.sigSPORE, tiponly = T) %>% 
   pull(label)
    
sig.b <- d.rax %>% 
  offspring(.node = mrca.sigB, tiponly = T) %>% 
  pull(label)

sig.fgek <- all.spore[!all.spore %in% sig.b]
sig.fgek_phage <- 
  sig.fgek[str_detect(sig.fgek,"phage")] %>% 
  str_remove("-phage")

all.phage <- 
  d.rax %>% 
  filter(group == "phage") %>% 
  pull(label) %>% 
  str_remove("-phage")




# Get alignment duplicate list -------------------------------------------------------
# these were removed automatically when constructing phylogeny
log <- readLines(here("phylo-clust/data/align-trim-tree/check_msa/check-msa.raxml.log"))
log <- log[str_detect(log, "identical")]
dups <- str_extract_all(log,"(.P_[0-9]*..-(phage|bacteria))", simplify = T) %>% 
  as_tibble() %>% 
  rename(kept = V1, removed = V2)
# mark sporulation classisfication
dups <- dups %>% mutate(phylo_spor = kept %in% sig.fgek_phage)


# Mark protein cluster classifications ------------------------------------
d.phage_clstrs <- d.new %>% 
  filter(protClstr_rep) %>% 
  select(protein, protClstr)

d.phage_clstrs <- 
  d.phage_clstrs %>% 
  mutate(phylo_spor = case_when(
    protein %in% sig.fgek_phage ~ TRUE,
    protein %in% all.phage ~ FALSE,
    TRUE ~ NA))

# mark alignment duplicates
for (i in 1:nrow(dups)){
  prot <- dups$removed[i] %>% 
    str_remove("-phage")
  
  d.phage_clstrs <- d.phage_clstrs %>% 
    mutate(phylo_spor = if_else(protein == prot, dups$phylo_spor[i], phylo_spor))
  
}

d.new <- 
  d.phage_clstrs %>% 
  select(protClstr, phylo_spor) %>% 
  left_join(d.new, ., by = "protClstr")



# summary -----------------------------------------------------------------

d.new %>% 
  # filter(protClstr_rep) %>% 
  group_by(phylo_spor, sig.class, phylum  ) %>% 
  summarise(n=n())

spor_sigs <- d.new %>% 
  filter (phylo_spor) %>% 
  filter(sig.class == "sporulation")


# host genera
spor_sigs %>% 
  mutate(host_genus = str_remove(`host name`, " .*")) %>% 
  group_by(host_genus) %>% 
  summarise(n=n())

# viral family
spor_sigs %>% 
  mutate(vir_fam = 
           str_remove(`virus lineage`, ".*Caudovirales; ") %>% 
           str_remove(";.*")) %>% 
  group_by(vir_fam) %>% 
  summarise(n=n())


# overall viral spp.
d.new%>% 
  group_by(taxon) %>% 
  summarise(n_sig=n()) %>% 
  ungroup() %>% 
  group_by(n_sig) %>% 
  summarise(n = n())
