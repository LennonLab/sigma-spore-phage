#setwd("/N/u/danschw/Carbonate/GitHub/sigma-spore-phage")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)
library(ggtreeExtra) #https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html
library(ggnewscale)
library(cowplot)
library(viridisLite)


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


# plot  sigmas all VOG phages -------------------------------------------------

# d.sp from vogdb/C_phage host
vog.sp <- read_csv(here("vogdb","data","vog_phages_Whost.csv"))
 

  

# summarize sigma content per phage
new.sp <-
  d.new %>% 
  group_by(taxon) %>% 
  summarise(n_sigma = n(), pred_spor = str_c(pred_spor, collapse = ";")) %>% 
  mutate(has_sigma = TRUE, spor_sigma = str_detect(pred_spor,"phylo\\+HMM")) %>% 
  #remove genome cluster duplicates from spp list and add sigma data
  right_join(., filter(vog.sp, is_rep),by = c("taxon" = "tax.id" )) %>% 
  # mark phages without sigma as such
  mutate(has_sigma = if_else(is.na(has_sigma), FALSE, has_sigma),
         sigma = case_when(
           spor_sigma ~ "spore-like sigma",
           has_sigma ~ "w. sigma",
           TRUE ~ "no sigma"
         )) %>% 
  # phylum name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) 


# sum phages by host phylum
d.plot <- new.sp %>% 
  group_by(phylum, sigma) %>% 
  summarise(n_phages = n()) 

# add total of phages
d.plot <-
  d.plot %>% 
  group_by(phylum) %>% 
  summarise(n_total = sum(n_phages)) %>% 
  left_join(d.plot,.)


#organize for plotting
d.plot <- d.plot %>% 
  filter(n_total>10)
  

# order of phyla by phage number
phyl_fct <- d.plot %>% 
  select(phylum, n_total) %>% 
  distinct() %>% 
  arrange(n_total) %>% 
   pull(phylum)

#plotting order for sigmas
sig_fct <- c("no sigma","w. sigma","spore-like sigma")

p <- 
  d.plot %>% 
  mutate(phylum = factor(phylum, levels=phyl_fct)) %>% 
  mutate(sigma = factor(sigma, levels=sig_fct)) %>% 
  # plot
  ggplot(aes(phylum, n_phages ,fill = sigma %>% fct_rev() ))+
  geom_col()+
  coord_flip()+
  scale_fill_viridis_d()+
  theme_classic()+
  theme(legend.position = c(0.8,0.25),
        legend.key.size = unit(0.2, "cm"),
        legend.background = element_blank())+
  guides(fill = guide_legend(title = "Sigmas"))+
  xlab ("Host Phylum")+
  ylab ("No. phage genomes")

ggsave2(here("phylo-clust","plots","spor_sigma_HostPhylum.png"),
        p, #+theme(legend.position = "none"),
        width = 4,height = 2)

# _________________----------------------------------------------
# Focus on phages of Firmicutes --------------------------

d.nsig <-
d.new %>%
  filter(phylum == "Firmicutes" ) %>% 
  mutate(spor_sigma = str_detect(pred_spor,"phylo\\+HMM")) %>% 
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
  
  



# extract viral family
    d.nsig <- d.nsig %>%
      mutate(viral.family=str_extract(`virus lineage`,
                                      regex("caudovirales;.*", ignore_case = T)))%>%
      mutate(viral.family=str_remove(viral.family,
                                     regex("caudovirales; ", ignore_case = T)))%>%
      mutate(viral.family=str_remove(viral.family,
                                     regex(";.*", ignore_case = T))) 
    
# spore forming genera by Bergey's
    bergey_spore <- read_csv(here("phylo-clust/data/genus_spore.csv"))
    
  d.nsig <- 
    bergey_spore %>% 
      select(genus.plot = genus, host_spore_former = spore_former) %>% 
      left_join(d.nsig, .)
    
# > Plot ---------------------------------------------------------------
  
  
  lab_text_size=2.5

    l.plot <- list()

    for (vfam in unique(d.nsig$viral.family)){
l.plot[[vfam]] <-
      d.nsig %>%

      filter(viral.family ==vfam) %>%
      #order phages
      arrange(desc(host_spore_former), desc(sig.position),  genus.plot) %>% 
      mutate(genus.plot = fct_inorder(genus.plot)) %>%
      mutate(host_spore_former = 
               factor(host_spore_former, levels = c(TRUE,FALSE))) %>% 
  
      ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
      geom_tile(aes(fill=spor_sigma), color="black")+
      theme_classic()+
      panel_border(color = "black", size=0.1)+
      scale_fill_manual(values = viridis(2)[2:1])+
      facet_grid(genus.plot~., scales = "free", space = "free")+
      geom_text(aes(label = genus.plot, color = host_spore_former),
                size=lab_text_size, show.legend = F,
                x = Inf, y = Inf, hjust = 1.1, vjust = 1.05) +
      scale_color_manual(values = c("blue","red"))+
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = "none",
            axis.title.y = element_blank(),
            panel.spacing = unit(0, "lines"),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8))+
      xlab("Sigma Factor Gene")+
      expand_limits(x=5)+
  scale_x_continuous(breaks = c(1:3))+
  ggtitle(vfam)
    }

#count viral families for proportional plotting
n.fams <- d.nsig %>% group_by(viral.family) %>% summarise(n=n())

pA <- egg::ggarrange(l.plot$Siphoviridae+
                       theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
                     l.plot$Podoviridae+
                       theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
                     l.plot$Myoviridae+
                       theme(legend.position = "bottom"),

   
          ncol = 1,heights = c(61,6,10))



ggsave2(here("phylo-clust/plots/","sigma_TIGR_content_Firmi.png"),
        plot = plot_grid(l.plot$Herelleviridae, pA , ncol = 2),
        width = 8, height = 10)

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
