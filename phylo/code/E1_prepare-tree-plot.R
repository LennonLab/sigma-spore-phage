#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)
library(ggtreeExtra) #https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
library(ggnewscale)
library(ggstar)



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
  


# root at common ansector of ECF -----------------------------------------------------

### identify ECF nodes

# B. subtilis ECF 
ecf <- c("sigM", "sigV", "sigW","sigX", "sigY","sigZ","ylaC") 
ecf.nodes <- d.rax %>% 
  filter(str_detect(sp, "subtilis")) %>% 
  filter(symbol %in% ecf) %>% 
  pull(node)

# Add other bacterial ECF by description
ecf.nodes <- d.rax %>% 
  filter(str_detect(group, "bacteria")) %>% 
  filter(str_detect(description, regex("ecf", ignore_case = T)) |
           str_detect(description, regex("extracyto", ignore_case = T))) %>% 
  pull(node) %>% 
  c(., ecf.nodes) %>% 
  unique()

# as.treedata(d.rax) %>% 
#   ggtree()+
#   geom_point2(aes(subset = (node %in% ecf.nodes)), size=2, color = "red")+
#   geom_point2(aes(subset = (node %in% root_ecf)), size=2, color = "blue")+
# geom_point2(aes(subset = (node == 656)), size=2, color = "blue")


# new root node at most recent common ancestor of ECF
root_ecf <- MRCA(d.rax, ecf.nodes) %>% pull(node)

# assign root and add to metadata
# rax <- root(rax, node = root_ecf)
rax.fbp <- root(rax.fbp, node = root_ecf, resolve.root = T)
rax.tbe <- root(rax.tbe, node = root_ecf, resolve.root = T)

d.rax <- as_tibble(rax.tbe) %>%
  select(node, label) %>%
  left_join(as_tibble(rax.fbp), . , by = c("node"),suffix = c(".fbp", ".tbe")) %>%
  mutate(group = case_when( str_detect(label.fbp, "bacteria") ~ "bacteria",
                            str_detect(label.fbp, "phage") ~ "phage",
                            TRUE ~ "NA")) %>%
  mutate(label = if_else(is.na(group), "", label.fbp)) %>%
  mutate(protein.id = str_remove(label.fbp, "-.*")) %>%
  left_join(., d.meta, by = "protein.id")

# replot to verify re-rooting
d.rax %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  mutate(ecf = groupOTU(d.rax,ecf.nodes) %>% pull(group)) %>% 
  as.treedata() %>% 
  ggtree(branch.length = "none")+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = 2)+
  geom_fruit(geom = "geom_tile",
             mapping=aes(y=ID, fill=ecf),color = "transparent")+
  scale_fill_viridis_d()+
  new_scale_fill()+
  geom_fruit(geom = "geom_tile",
             mapping=aes(y=ID, fill=group),color = "transparent")+
  new_scale_fill()



# prepare plot layers -----------------------------------------------------


# >> branches -----------------------------------------------------------------

# make informative label
d.rax <- d.rax %>% 
  mutate(tip.label = case_when(group == "phage" ~ sp,
                               group == "bacteria" ~ paste(sp, symbol, sep = "_"))) %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) 


# function to mark phage only monophyletic clades
# go over all internal nodes, divide to clade and check if one side is phage only
# if yes mark as phage clade node

internal <- d.rax %>% 
  filter(node > as.phylo(.) %>% Ntip()) %>% 
  pull(node)

#assign phage tip nodes to clade, and all other as empty vector
d.rax <- d.rax %>% 
  
  mutate(clade = case_when( str_detect(label, "bacteria") ~ "bacteria",
                            str_detect(label, "phage") ~ "phage",
                            TRUE ~ "bacteria"))

for (i in internal){
  x <- groupClade(rax.fbp,.node=i)
  # x <- groupClade(rax,.node=i)
  
  # add groups to main tree
  d.x <- as_tibble(x) %>% 
    select(node, cur.split = group) %>% 
    left_join(d.rax,., by = "node")
  
  # get groups as charcater vectors
  g1 <- d.x %>% 
    filter(cur.split == 1) %>% 
    filter(! node %in% internal) %>% 
    pull(group)
  
  g2 <- d.x %>% 
    filter(cur.split == 0) %>% 
    filter(! node %in% internal) %>% 
    pull(group)
  
  # test if any of the groups is phage only
  if (length(g1)==0) next # avoid empty charcter returning TRUE
  if (length(g2)==0) next
  if (all(g1=="phage") | all(g2=="phage")){
    d.rax$clade[d.rax$node==i] <- "phage" # assign to phage only clade
  }
}


# >> tigr -----------------------------------------------------------------


#################
# add tigr classification
tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")
d.rax <- d.rax %>% 
  #classify by spore function
  mutate(tigr.type = case_when(tigr.hit %in% tigr.spore ~ tigr.hit,
                               tigr.hit == "no_hit" ~ "no hit",
                               TRUE ~ "other")) %>% 
  mutate(tigr.type = str_remove(tigr.type, "spore_")) %>% 
  mutate(tigr.type = str_replace(tigr.type, "sigma", "sig")) %>% 
  mutate(tigr.type = str_replace(tigr.type, "Sig", "sig")) %>% 
  mutate(tigr.type = fct_relevel(tigr.type, "sigF", "sigG", "sigK", "sigE", "other", "no hit")) 
  


d.rax <- d.rax %>% 
  mutate(tbe = if_else(is.na(tip.label), label.tbe %>% parse_number(), NULL)) %>%
  mutate(support = case_when( #is.na(UFboot) ~ "NA",
                              100*tbe>= 90 ~ ">90%",
                              # 100*tbe>= 80 ~ ">80%",
                              100*tbe>= 70 ~ ">70%",
                              # 100*tbe>= 50 ~ ">50%",
                              TRUE ~ NA_character_)) %>%
  mutate(support = fct_rev(support)) %>%
  
  mutate(circ1.lab = case_when(str_detect(sp, "Escherichia virus T4 ") ~ "T4",
                               str_detect(sp, "Bcp1") ~ "Bcp1",
                               str_detect(sp, "Fah") ~ "Fah",
                               str_detect(sp, "SPO1") ~ "SPO1",
                               str_detect(sp, "SP-10") ~ "SP10",
                               str_detect(sp, "Goe3") ~ "Goe3",
                               str_detect(sp, "Eldridge") ~ "ELD",
                               TRUE ~ bs.label)) %>% 
  # height for mapping phage vs bacteria
  mutate(clade.x = case_when(clade == "phage" ~ 0,
                             clade == "bacteria" ~ 0.02,
                             TRUE ~ 0))

  
# Tree to identify nodes that define clades
d.rax %>% 
  as.treedata(.) %>% 
    ggtree(aes(color = clade), branch.length = "none", 
         layout = "fan", open.angle=5, size=0.1)+
  scale_color_manual(values = c("grey", "blue"))+
  geom_tiplab(aes(label = node), offset = 0.1, size=2)+
  geom_tiplab(aes(label = tip.label), offset = 1.5, size=2)+
  geom_nodelab(aes(label = node), color = "pink", size=2)+
  geom_fruit(geom = "geom_text", 
             mapping = aes(label=bs.label), color="grey20", offset = 0.1)
ggsave(here("phylo","plots","sigma_circle_NODES.pdf"), height=15, width = 15)

mrca.sigA <- MRCA(d.rax, c(409,440,380)) %>% pull(node)
mrca.gp55 <- MRCA(d.rax, c(469,566)) %>% pull(node)
mrca.gp34 <- MRCA(d.rax, c(185,205)) %>% pull(node)
mrca.ecf <- MRCA(d.rax, c(177,54)) %>% pull(node)
mrca.sigBFG <- d.rax %>% 
  filter(bs.label %in% c("sigF", "sigG", "sigB")) %>%
  pull(node) %>% 
  MRCA(d.rax, .) %>% pull(node)

# base tree
p1 <-  
  as.treedata(d.rax) %>% 
  
  ggtree(aes(color = clade), #branch.length = "none",
         layout = "fan", open.angle=5, size=0.1)+
  scale_color_manual(values = c("grey", "blue"), 
                     guide = guide_legend(title = "source" ))+
# bootstrap support
  geom_nodepoint(aes(fill=support), colour = "transparent", size=0.5, shape = 21)+
  scale_fill_manual(values = c("red", "pink"), na.value = NA,
                    guide = guide_legend(override.aes = list(size=5)))
 

p2 <- p1 +
  # circle 1 :phage vs bacteria
  new_scale_fill()+
  new_scale_color()+
  geom_fruit(geom = "geom_tile",
             mapping=aes( fill=clade),
             width = 0.2, color = "transparent")+
  scale_color_manual(values = c("grey", "blue"), guide = guide_legend(title = "source" ))+
  scale_fill_manual(values = c("grey", "blue"), guide = guide_legend(title = "source" )) +
  
  # circle 2: Bacterial/ host taxonomy
  
  new_scale_fill()+
  new_scale_color()+
  geom_fruit(geom = "geom_tile",
             mapping=aes( fill=host_phylum),
             width = 0.2, color = "transparent",offset = 0.05)+
  scale_fill_discrete(guide = guide_legend(ncol = 2))+
  
  # circle 3: TIGR
  
  new_scale_fill()+
  new_scale_color()+
  geom_fruit(geom = "geom_tile",
             width = 0.2, mapping=aes(fill=tigr.type), offset = 0.05)+
  # B. subtilis tip labels
  geom_fruit(geom = "geom_text", size=2,
             mapping = aes(label=bs.label), color="grey20", offset = 0.15)+
  scale_fill_viridis_d(guide = guide_legend(ncol = 2))
  

# clade annotations 
p3 <- p2+
  geom_cladelabel(node=mrca.ecf, label = "ECF", fontsize = 2.5,
                  offset = 4, offset.text = 0.2 )+
  geom_cladelabel(node=mrca.gp55, label = "phage\n(gp55)", fontsize = 2.52,
                  offset = 3, offset.text = 0.8)+
  geom_cladelabel(node=mrca.gp34, label = "phage\n(gp34)", fontsize = 2.5,
                  offset = 2, offset.text = 0.8)+
  geom_cladelabel(node=mrca.sigA, label = "sigA", fontsize = 2.5,
                  offset = 2.2, offset.text = 0.2 )+ 
  geom_cladelabel(node=mrca.sigBFG, label = "sigB/F/G",fontsize = 2.5,
                  offset = 3.5, offset.text = 0.2 )


ggsave(filename = here("phylo","plots","sigma_circle_rooted.pdf"),
       plot = p3,#+theme(legend.position = "none"),
       height=6, width = 8)


# trim margins ------------------------------------------------------------
# https://yulab-smu.top/treedata-book/faq.html#circular-blank

# library(magick)
# library(pdftools)
#     
# x <- image_read_pdf(here("phylo","plots","sigma_circle_rooted.png"), density = 1200)
# y <- image_trim(x)
# 
# ggsave(filename = here("phylo","plots","sigma_circle_rooted.png"),
#        plot = image_ggplot(y),
#        height=6, width = 8)
  
  