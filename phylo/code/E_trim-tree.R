# setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

# In the analysis below the RAxML tree is used to identify and remove clades of bacterial only sigma factors
# This is done to reduce the number of sigma factor genes analyzed

#import tree
rax <- read.tree(list.files(path = here("phylo/data/align-trim-tree/tree"),
                            pattern = "bestTree", full.names = T))


# convert to table of tree data
d.rax <- as_tibble(rax) %>% 
  mutate(group = if_else( str_detect(label, "bacteria"), "bacteria", "phage")) %>% 
  mutate(protein.id = str_remove(label, "-.*"))


# add info on bacterial genes 
d.bact<- read_csv( here("phylo/data/bacterial_features.csv")) %>% 
  select(protein.id = product_accession, sp=sp, description=name, symbol)
 
# add bacterial metadata to tree tibble
d.rax <-
  left_join(d.rax, d.bact, by = "protein.id")

# # make informative label
d.rax <- d.rax %>%
  mutate(tip.label = case_when(str_detect(description, regex("ecf", ignore_case = T)) ~ "ECF",
                               str_detect(sp, regex("subtilis", ignore_case = T)) ~ symbol,
                               TRUE ~ ""))
tree <- as.treedata(d.rax)
# Plot tree
p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=5, shape=20, alpha=.5)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)

ggsave(here("phylo","plots","sigma_all_unrooted.pdf"),p, height=10, width = 10)

# there are two clades that can be reduced:
# > The ECF clade - has no phage sequences
# > the RpoD/sigA clade - has no phage sequences

# identify splitting nodes
p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node), size =1.5, color = "red") 

ggsave(here("phylo","plots","sigma_nodeNUMS_unrooted.pdf"),p, height=10, width = 10)


ecf_base_node <- 522
rpod_base_node <- 440


# split tree by nodes found above 
x <- groupClade(rax,.node=c(rpod_base_node, ecf_base_node))
# plot to check
ggtree(x, layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=5, shape=20, alpha = .5)


#need to keep split == 1
split_keep <- 1


# add split groups to main tree
d.x <- as_tibble(x) %>% 
  select(node, split = group)

d.rax <- full_join(d.rax, d.x) 

tree <- as.treedata(d.rax)

p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color = group, shape=split),
                size=3, alpha = .5)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)


ggsave(here("phylo","plots","sigma_TrimGroups_unrooted.pdf"),p, height=10, width = 10)




#split by group, and add sigma of phages, B. subtilis and . E. coli from removed groups, for reference.
d.keep <- d.rax %>% 
  filter(split == split_keep | 
           #group == "phage" |
           str_detect(sp, regex("subtilis", ignore_case = T))|
           str_detect(sp, regex("Escherichia_coli", ignore_case = T)))

keep <- d.keep %>%
  filter(!is.na(protein.id)) %>% 
  pull(protein.id)

#check removed
d.rax %>% 
  filter(! protein.id %in% keep ) %>% 
  pull(group) %>% table()
# 183 bacteria sequences were removed.
# no phage sequences

##### Save results #####
if (!dir.exists(here("phylo","data", "reduced_set_to_align"))){
  dir.create(here("phylo","data", "reduced_set_to_align"))
}
# filter multifasta of hsearch results for filtering
library(seqinr)
sigma_fa <- read.fasta(here("phylo","data/sigmas_to_align.faa"),
                       seqtype =  "AA", whole.header = TRUE)

header.protein <- names(sigma_fa) %>% 
  str_remove("-bacteria") %>% 
  str_remove("-phage")

keep_fa <- sigma_fa[c(which(header.protein %in% keep))]

write.fasta(sequences = getSequence(keep_fa),
            names = getName(keep_fa),
            file.out = here("phylo","data", "reduced_set_to_align","sigmas_to_align.faa"))

# make filtered versions of meta data tables

# load viral sigmas from vog HMM analysis
# load(here("vogdb/data/vog_sigma_clean_Whost.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
# d.faa <- d.faa %>%
#   filter(protein %in% keep)
# write_csv(d.faa %>% select(-seq), here("vogdb/data", "reduced_set_to_align","sigmas_to_align.csv"))

# d.sp<- read_csv( here("phylo/data/bacterial_features.csv"))
# d.b <- d.sp %>% 
#   filter(product_accession %in% keep)
# write_csv(d.sp, here("phylo","data", "reduced_set_to_align","bacterial_features.csv"))

