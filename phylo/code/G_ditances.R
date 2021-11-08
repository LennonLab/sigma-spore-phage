library(here)
library(tidyverse)
library(cowplot)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)


# import trees ------------------------------------------------------------

# Best ML tree with Felsenstein bootstrap (FBP) support values
rax.fbp <- read.tree(here("phylo/data/align-trim-tree",
                          "tree/ML_TBE_tree.raxml.supportFBP"))
d.rax.fbp <- as_tibble(rax.fbp)

# Best ML tree with Transfer bootstrap (TBE) support values
rax.tbe <- read.tree(here("phylo/data/align-trim-tree",
                          "tree/ML_TBE_tree.raxml.supportTBE"))
d.rax.tbe <- as_tibble(rax.tbe)

#join the supports into one df
d.rax <- d.rax.tbe %>%
  select(node, label) %>%
  left_join(d.rax.fbp, . , by = c("node"),suffix = c(".fbp", ".tbe")) %>%
  mutate(group = case_when( str_detect(label.fbp, "bacteria") ~ "bacteria",
                            str_detect(label.fbp, "phage") ~ "phage",
                            TRUE ~ "NA")) %>%
  mutate(label = if_else(is.na(group), "", label.fbp)) %>%
  mutate(protein.id = str_remove(label.fbp, "-.*"))



# Add metadata ------------------------------------------------------------

# load viral sigmas data from vog HMM analysis
load(here("vogdb/data/vog_sigma_clean_tigr.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
d.phage <- d.faa %>%
  filter(protein %in% d.rax$protein.id) %>% 
  rename(protein.id = protein) %>% 
  #adjust host phyla names to shorten
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum))
rm(d.faa)

# load data collected for bacteria from feature tables
d.bact <- read_csv( here("phylo/data/bacterial_features.csv"))
d.bact <- d.bact %>%
  filter(product_accession %in% d.rax$protein.id) %>% 
  rename(protein.id = product_accession, description=name) %>% 
  # remove strain designation for compactness
  mutate(check_sp = str_count(sp, "_")) %>% 
  mutate(sp = case_when(check_sp == 1 ~ sp,
                        check_sp > 1  ~ str_extract(sp, regex("^.*?_.*?_")))) %>% 
  mutate(sp = str_remove(sp, "_$")) %>% 
  select(-check_sp)


# add bacterial taxonomy data. collected by phylo/code/get_taxa.R
d.bac_tax <- read_csv(here("phylo","data","bacterial_taxonomy.csv"))

#add bacteria tigr classification
d.bac_tigr <- read_csv(here("phylo","data","bacterial_sigma_tigr.csv"))

d.bact <- 
  d.bac_tax %>% 
  left_join(d.bact, ., by = c("assembly" = "assembly_accession")) %>%
  left_join(., d.bac_tigr, by = c("protein.id" = "query_name")) 

# add meta data to tree tibble
d.meta <-
  d.phage %>% 
  select( protein.id, description, sp, host_phylum=phylum,tigr.hit, sig.class) %>% 
  mutate(symbol=NA) %>% 
  bind_rows(., d.bact %>% select( protein.id, description, host_phylum=phylum, tigr.hit, sp, symbol))

d.rax <- left_join(d.rax, d.meta, by = "protein.id")



# distance ----------------------------------------------------------------

dist.tips <- 
  d.rax %>% 
  as.treedata() %>% 
  as.phylo() %>% 
  ape::cophenetic.phylo()

# convert to data frame
dist.tips <-
  as_tibble(dist.tips, rownames = "tip1") %>% 
    pivot_longer(cols = where(is.numeric), 
                 names_to = "tip2", values_to = "distance") %>% 
  #delete diaganol (distance = 0)
  filter(tip1 != tip2)


# select tips of interest
sub <- d.rax %>% 
  filter(str_detect(sp, regex("Eldridge", ignore_case = T))|
           str_detect(sp, regex("SP-10", ignore_case = T))|
           str_detect(sp, regex("Goe3", ignore_case = T))|
           str_detect(sp, regex("Bacillus_subtilis", ignore_case = T))|
           str_detect(sp, regex("Bacillus_cereus", ignore_case = T))|
           str_detect(sp, regex("Clostridioides_difficile", ignore_case = T))) %>% 
  filter(sig.class == "sporulation" | group == "bacteria") %>% 
  filter(group == "phage" | str_detect(tigr.hit, regex("sig.*[EKFG]"))|
           str_detect(symbol, "rpoD") ) %>% 
  pull(label)

dist.tips <- dist.tips %>% 
  filter((tip1 %in% sub) & (tip2 %in% sub)) %>% 
  #delete duplicates and separate phage and bacteria by column
    # each comparison is duplicated - from full matrix
  filter(str_detect(tip1, "phage") & str_detect(tip2, "bacteria"))  


# names of proteins
# protein.names <- read_csv(here("vogdb/data/hmm_align","seq_names.csv"))
protein.names <-d.rax %>% 
  filter(label %in% sub) %>% 
  #add sigK to Cdiff symbol
  mutate(symbol = if_else(tigr.hit == "spore_sigmaK", "sigK", symbol)) %>% 
  # phage protein names
  mutate(name = case_when(protein.id == "YP_009274875.1" ~ "ELDg168",
                          protein.id == "YP_009274876.1" ~ "ELDg169",
                          protein.id == "YP_007003377.1" ~ "SP10",
                          protein.id == "YP_009832138.1" ~ "Goe3",
                          TRUE ~ sp)) %>% 
  # bacterial names
  mutate(name = case_when(name == "Bacillus_subtilis" ~ paste0("Bsub_",symbol),
                          name == "Bacillus_cereus" ~ paste0("Bcer_",symbol),
                          name == "Clostridioides_difficile" ~ paste0("Cdif_",symbol),
                          TRUE ~ name)) %>% 
  select(protein.id, name)
  
  
  
# replace names
dist.tips <- dist.tips %>% 
  mutate(tip1 = str_remove(tip1, "-.*")) %>% 
  mutate(tip2 = str_remove(tip2, "-.*")) %>% 
  left_join( . , protein.names, by = c("tip1" = "protein.id")) %>% 
  mutate(tip1 = name) %>% 
  rename(phage_gene = tip1) %>% 
  select(-name) %>% 
  left_join(., protein.names, by = c("tip2" = "protein.id")) %>% 
  mutate(tip2 = name) %>% 
  rename(bacterial_gene = tip2) %>% 
  select(-name)


# plot --------------------------------------------------------------------


p <- dist.tips%>% 
  mutate(bacterial_gene = str_replace(bacterial_gene, "rpoD", "sigA")) %>% 
  mutate(host_gene = str_extract(bacterial_gene, "sig.")) %>% 
  mutate(host_gene = factor(host_gene, levels = c("sigF","sigG","sigA","sigK","sigE"))) %>% 
  mutate(bacterial_gene = str_remove(bacterial_gene, "_.*") %>% 
           fct_relevel("Bsub", after = 0)) %>% 
  mutate(phage_gene = factor(phage_gene, levels = c("ELDg169","ELDg168","Goe3","SP10"))) %>% 
  
  ggplot(aes(bacterial_gene, fct_rev(phage_gene)))+
  geom_tile(aes(fill = distance))+
  geom_text(aes(label = distance %>% round(digits = 2)))+
  scale_fill_viridis_c(direction = 1,
                       guide = guide_colourbar(title = "Branch Length"))+
  facet_wrap(~host_gene)+
  ylab("Phage Genes")+
  xlab("Bacterial sp.")+
  theme_classic()+
  theme(legend.position = c(0.8,0.2))+
  panel_border(color = "black")

save_plot(here("phylo/plots","distance_clonedGenes.png"),p)


# sporulation vs distance -------------------------------------------------

# import sporulation data
d.spor <- read_csv(here("phylo/data/FCM_sporulation_induction.csv")) %>% 
  filter(strain != "Empty Vector")


d.both <- dist.tips %>% 
  filter(str_detect(bacterial_gene, "Bsub_sigG")) %>% 
  left_join(d.spor, . , by = c("strain" = "phage_gene")) %>% 
  filter(! str_detect(strain, "sigF")) %>% 
  mutate(distance = if_else(strain == "sigG", 0, distance))


p <- d.both %>% 
  ggplot(aes(distance, induction.ratio)) +
  geom_smooth(method  = "lm")+
  geom_point(aes(fill  = pnl),shape=21, size = 2)+
  theme_classic()+
  panel_border(color = "black")+
  scale_shape_manual(values = 21:25)+
  scale_y_continuous(labels = scales::percent, trans = "log10")+
  scale_fill_manual(values = c("grey","white"))+
  xlab("Phylogenetic distance (from B. subtilis sigG)")+
  ylab("Spore yield")+
  guides(fill = guide_legend(title = "source"))

save_plot(here("phylo/plots","distance_yield.png"),p)

cor.test(x = d.both$distance, d.both$induction.ratio)
