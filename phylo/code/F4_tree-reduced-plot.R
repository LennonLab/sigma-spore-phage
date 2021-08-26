#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)
library(ggtreeExtra) #https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
library(ggnewscale)
library(ggstar)



# import trees ------------------------------------------------------------

# Best ML tree with Felsenstein bootstrap (FBP) support values 
rax.fbp <- read.tree(here("phylo/data/reduced_set_to_align",
                          "tree/ML_TBE_tree.raxml.supportFBP"))
d.rax.fbp <- as_tibble(rax.fbp)

# Best ML tree with Transfer bootstrap (TBE) support values
rax.tbe <- read.tree(here("phylo/data/reduced_set_to_align",
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




####load metadata for sequences ####

# load viral sigmas data from vog HMM analysis
load(here("vogdb/data/vog_sigma_clean_tigr.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
d.phage <- d.faa %>%
  filter(protein %in% d.rax$protein.id) %>% 
  rename(protein.id = protein)
rm(d.faa)

# load data collected for bacteria from fearure tables
d.bact<- read_csv( here("phylo/data/bacterial_features.csv"))
d.bact <- d.bact %>%
  filter(product_accession %in% d.rax$protein.id) %>% 
  rename(protein.id = product_accession, description=name)

# add meta data to tree tibble
d.meta <-
  d.phage %>% 
    select( protein.id, description, sp, tigr.hit, sig.class) %>% 
    mutate(symbol=NA) %>% 
  bind_rows(., d.bact %>% select( protein.id, description, sp, symbol))
  
d.rax <- left_join(d.rax, d.meta, by = "protein.id")
  

#### root at base of ECF ####

# find node
p <- d.rax %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node), color="red", size=3)

  ggsave(here("phylo","plots","reduced-sigma_nodeNUMS_unrooted.pdf"), 
         plot = p, height=10, width = 10)


# new root node
root_ecf <- 315

# assign root and add data
rax.fbp <- root(rax.fbp, node = root_ecf)
rax.tbe <- root(rax.tbe, node = root_ecf)

#join the supports into one df
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
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=tigr.hit, shape = group), size=2)+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = .1)

#### Plotting ####
# make informative label
d.rax <- d.rax %>% 
  mutate(tip.label = case_when(group == "phage" ~ sp,
                               group == "bacteria" ~ paste(sp, symbol, sep = "_"))) %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) 

#################
# function to mark phage only monophyletic clades
# go over all internal nodes, divide to clade and check if one side is phage only
# if yes mark as phage lade node

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


#################
# add tigr classification
tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")
d.rax <- d.rax %>% 
  #classify by spore function
  mutate(tigr.type = if_else(tigr.hit %in% tigr.spore,tigr.hit, "other" )) %>% 
  mutate(tigr.type = str_remove(tigr.type, "spore_")) %>% 
  mutate(tigr.type = str_replace(tigr.type, "sigma", "sig")) %>% 
  mutate(tigr.type = str_replace(tigr.type, "Sig", "sig")) %>% 
  mutate(tigr.type = fct_relevel(tigr.type, "sigF", "sigG", "sigK", "sigE", "other")) %>%
  mutate(tigr.type=if_else(group == "phage",tigr.type,NULL)) 
  

# # first tree
ggtree(as.treedata(d.rax) , aes(color = clade))+
  geom_tippoint(aes(shape=group), size=2, shape=20)+
  geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)

# # # unrooted trees
# p <-
#   ggtree(as.treedata(d.rax), layout = 'equal_angle')+
#   geom_tippoint(aes(color=group), size=2, shape=20)+
#   geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)
# 
# # ggsave(here("phylo","plots","sigma_reduced_unrooted.pdf"),p, height=10, width = 10)
# 
# 
# as.treedata(d.rax) %>% 
#   ggtree( layout = 'equal_angle', aes(color = clade))
# 
# as.treedata(d.rax) %>% 
#   ggtree( layout = 'radial', aes(color = clade))+
#   geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = 1)
# 
# p <-as.treedata(d.rax) %>% 
#   ggtree(aes(color = clade))+
#   # geom_tippoint(aes(fill=group), size=1, shape=20)+
#   geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)+
#   scale_color_manual(values = c("grey30", "blue"))
# 
# # ggsave(here("phylo","plots","sigma_reduced_rooted.pdf"),p, height=10, width = 10)


d.rax %>% 
  mutate(tbe = if_else(is.na(tip.label), label.tbe %>% parse_number(), NULL)) %>% 
  mutate(support = case_when( #is.na(UFboot) ~ "NA",
                              100*tbe>= 90 ~ ">90%",
                              # 100*tbe>= 80 ~ ">80%",
                              100*tbe>= 70 ~ ">70%",
                              100*tbe>= 50 ~ ">50%",
                              TRUE ~ NA_character_)) %>% 
  mutate(support = fct_rev(support)) %>% 
  
  mutate(circ1.lab = case_when(
                               str_detect(sp, "SP-10") ~ "SP10",
                               str_detect(sp, "Goe3") ~ "Goe3",
                               str_detect(sp, "Eldridge") ~ "ELD",
                               TRUE ~ bs.label)) %>% 
  as.treedata(.) %>% 
  ggtree(aes(color = clade), 
         layout = "fan", open.angle=10, size=0.5)+
  geom_fruit(geom = "geom_tile",
             mapping=aes(y=ID, fill=clade),color = "white",
             offset = 0)+
  geom_fruit(geom = geom_text,
             mapping = aes(y=ID, label=circ1.lab), color = "black",
             offset = 0, size=2)+
  # geom_fruit(geom = geom_text,
  #            mapping = aes(y=ID, label=phage.lab), color = "black",offset = .1)+
  new_scale_fill()+
  geom_fruit(geom = "geom_tile",
             mapping=aes(y=ID, fill=tigr.type), color = "white",
             offset = 0.18)+
  scale_fill_viridis_d()+
  new_scale_fill()+
  geom_nodepoint(aes(fill=support), colour = "transparent", size=1.5, shape = 21)+
  scale_fill_grey(na.translate =FALSE)+
  geom_tiplab(aes(label=tip.label), color="red",alpha = 0.5,
              size=1, offset = .1)

ggsave(here("phylo","plots","sigma_reduced_rooted.pdf"), height=10, width = 10)

  # 
  # geom_tippoint(aes(shape=tigr.type), size=2, shape = 20, shape = c(21,23))+
  # geom_nodepoint(aes(fill=support, shape = support), size=3)+
  # 
  # geom_tiplab(aes(label=tip.label, color = group), alpha = 0.5,
  #             size=2, offset = .1, show.legend=F)+
  # geom_tiplab(aes(label=bs.label), color="red",alpha = 0.5,
  #             size=2, offset = .1)+
  # 
  # scale_color_manual(values = c("grey30", "blue"))+
  # scale_shape_manual(values = c(NA,21,21,21,21))
  # 
  # ggsave(here("phylo","plots","sigma_reduced_rooted.pdf"), height=10, width = 10)
  # 

