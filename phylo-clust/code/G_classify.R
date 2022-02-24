#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)
library(ggtreeExtra) #https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
library(ggnewscale)
library(cowplot)


# load data ---------------------------------------------------------------

load(here("phylo-clust/data/plotted_tree.Rdata"))
load(here("phylo-clust/data/faa_table_clustered.RData"))

# load viral sigmas data from vog HMM analysis
load(here("TIGR/data/vog_sigma_clean_tigr.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
# define sporulation-related sigma TIGRFAMs
tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")

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

# get all proteins in sigma B/F/G/E/K clade (phage and bacteria)
all.spore <- d.rax %>% 
   offspring(.node = mrca.sigSPORE, tiponly = T) %>% 
   pull(label)

# get all proteins in sigma B clade (phage and bacteria)
# This is not a sporulation sigma, nested in sporulation sigmas
sig.b <- d.rax %>% 
  offspring(.node = mrca.sigB, tiponly = T) %>% 
  pull(label)

# remove sigB from all spore sigmas to get sporulation sigmas
sig.fgek <- all.spore[!all.spore %in% sig.b]

# remove bacteria from above to get phage.
sig.fgek_phage <- 
  sig.fgek[str_detect(sig.fgek,"phage")] %>% 
  str_remove("-phage")

# get all phage proteins included in phylogeny
all.phage <- 
  d.rax %>% 
  filter(group == "phage") %>% 
  pull(label) %>% 
  str_remove("-phage")


# Get alignment duplicate list -------------------------------------------------------
  # these were removed automatically when constructing phylogeny
  # were marked in duplicates in trimmed alignment
log <- readLines(here("phylo-clust/data/align-trim-tree/check_msa/check-msa.raxml.log"))

log <- log[str_detect(log, "identical")]

dups <- str_extract_all(log,"(.P_[0-9]*..-(phage|bacteria))", simplify = T) %>% 
  as_tibble() %>% 
  rename(kept = V1, removed = V2)
  
# mark sporulation classisfication for duplicates
dups <- dups %>% mutate(phylo_spor = kept %in% sig.fgek_phage)


# Mark protein cluster classifications ------------------------------------
  # The assumption here is that clustered proteins
  # have the same phylogentic position

# get cluster numbers and thier represntative protein 
# (rep. proteins were used for phylogeny)
d.phage_clstrs <- 
  d.new %>% 
  filter(protClstr_rep) %>% 
  select(protein, protClstr)

# mark phylogentic sporulation classification
d.phage_clstrs <- 
  d.phage_clstrs %>% 
  mutate(phylo_spor = case_when(
    protein %in% sig.fgek_phage ~ TRUE, #phylo sporulation
    protein %in% all.phage ~ FALSE, # phylo non-sporulation
    TRUE ~ NA)) # protein rep. not in phylogeny (alignment duplicates)

# mark alignment duplicates (based on their matching proteins included in phylogeny)
for (i in 1:nrow(dups)){
  prot <- dups$removed[i] %>% 
    str_remove("-phage")
  
  d.phage_clstrs <- d.phage_clstrs %>% 
    mutate(phylo_spor = if_else(protein == prot, dups$phylo_spor[i], phylo_spor))
  
}

# apply phylo classification to full list of proteins
d.new <- 
  d.phage_clstrs %>% 
  select(protClstr, phylo_spor) %>% 
  left_join(d.new, ., by = "protClstr")

# combine phylogeny and TIGR predictions
d.new <- 
  d.new%>% 
  mutate(pred_spor = 
           case_when(
             (phylo_spor) & (sig.class == "sporulation") ~ "phylo+HMM",
             (!phylo_spor) & (sig.class != "sporulation") ~ "neither",
             (phylo_spor) & (sig.class != "sporulation") ~ "phylo",
             (!phylo_spor) & (sig.class == "sporulation") ~ "HMM",
           )) %>% 
  mutate(pred_spor = fct_relevel(pred_spor, "neither"))

# _________________----------------------------------------------
# Plot by phylum ---------------------------------------------


# summarize for labels
n.gene.phylum <- d.new %>%
  mutate(phylum=str_remove(phylum,"/.*"))%>%
  group_by(phylum)%>%
  summarise(n=n(),.groups = "drop") %>% 
  mutate(phyl.lab = paste0(phylum,"\n(n=",n,")"))

# summarize by spore prediction and host phylum
d.plot <- 
  d.new%>%
  group_by(phylum,pred_spor)%>%
  summarise(n=n(), .groups = "drop") %>% 
  group_by(phylum)%>%
  mutate(perc=n/sum(n)) %>% 
  
  mutate(phylum=str_remove(phylum,"/.*"))%>%
  mutate(phylum=fct_reorder(phylum,n)) %>% 
  left_join(., n.gene.phylum, by = "phylum")

p.phylum <- d.plot %>% 
  ggplot(aes(x=phyl.lab, y = perc,fill = pred_spor)) + 
  geom_bar(position="fill", stat="identity", color = "transparent", size = 0,width = 0.5) +
  
  xlab("Host Phylum") +
  ylab("Phage-encoded Sigma Factor Genes")+
  guides(fill = guide_legend("Predicted\nSporulation\nSigma\nFactor", reverse = T))+
  
  scale_y_continuous(labels=scales::percent, position = "left") +
  scale_fill_viridis_d(direction = -1) + 
  theme_classic(base_size = 8)+
  panel_border(color = "black")+
  coord_flip()

ggsave2(here("phylo-clust","plots","spor_sigma_HostPhylum.png"),
        p.phylum,#+theme(legend.position = "none"),
        width = 4,height = 2)


# _________________----------------------------------------------
# Focus on phages of Firmicutes --------------------------
d.nsig <-
d.new %>%
  filter(phylum == "Firmicutes" ) %>% 
  # arrange phage by n.sigmas
  mutate(`virus name` = fct_infreq(`virus name`)) %>%
  # assign a positional index to each sigma factor (ordered by prediction)
  group_by(taxon)%>%
  arrange(pred_spor)%>%
  mutate(sig.position=row_number())%>%
  ungroup()%>%
  mutate(`virus name`=str_remove(`virus name`, ".*phage "))%>%
  mutate(`virus name`=str_remove(`virus name`, ".*virus"))%>%
  mutate(`virus name`=as_factor(`virus name`)) %>%
  
  # extract host genus
  separate(family.etc, into = c("family","genus","species","strain"),sep=";") %>%
    # genus is also in viral sp name as first word
  mutate(genus2=str_extract(sp,regex(".*? ")) %>% trimws()) %>%
  mutate(genus = trimws(genus)) %>%
  mutate(genus.plot = if_else(is.na(genus), genus2, genus))



# > Plot ---------------------------------------------------------------
    # extract viral family
    # d.sp <-
    d.nsig <- d.nsig %>%
      mutate(viral.family=str_extract(`virus lineage`,
                                      regex("caudovirales;.*", ignore_case = T)))%>%
      mutate(viral.family=str_remove(viral.family,
                                     regex("caudovirales; ", ignore_case = T)))%>%
      mutate(viral.family=str_remove(viral.family,
                                     regex(";.*", ignore_case = T))) 
    lab_text_size=3

    l.plot <- list()

    for (vfam in unique(d.nsig$viral.family)){
l.plot[[vfam]] <-
      d.nsig %>%

      filter(viral.family ==vfam) %>%
      mutate(genus.plot = fct_infreq(genus.plot)) %>%
      # filter(! genus.plot %in% c("Bacillus", "Priestia")) %>%
      # filter(! genus.plot %in% c("Staphylococcus", "Enterococcus")) %>%

      
      ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
      geom_tile(aes(fill=pred_spor), color="black")+
      theme_classic()+
      panel_border(color = "black", size=0.1)+
      scale_fill_viridis_d(direction = -1, drop = FALSE)+
      facet_grid(genus.plot~., scales = "free", space = "free")+
      geom_text(aes(label = genus.plot), size=lab_text_size,
                x = Inf, y = Inf, hjust = 1.1, vjust = 1.05) +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = "none",
            axis.title.y = element_blank(),
            panel.spacing = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5))+
      xlab("Sigma Factor Gene")+
      expand_limits(x=5)+
  scale_x_continuous(breaks = c(1:3))+
  ggtitle(vfam)
    }

#count viral families for proportional plotting
n.fams <- d.nsig %>% group_by(viral.family) %>% summarise(n=n())

pA <- egg::ggarrange(l.plot$Myoviridae+
                       theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
               # ggplot()+ theme_void(),
               l.plot$Podoviridae+
                 theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
               # ggplot()+ theme_void(),
               l.plot$Siphoviridae,
          ncol = 1,heights = c(7,6,61))



ggsave2(here("phylo-clust/plots/","sigma_TIGR_content_Firmi.png"),
        plot = plot_grid(l.plot$Herelleviridae, pA , ncol = 2),
        width = 8, height = 11)

# summary -----------------------------------------------------------------

d.new %>% 
  # group_by(phylo_spor, sig.class, phylum  ) %>% 
  group_by(phylo_spor, sig.class) %>% 
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
