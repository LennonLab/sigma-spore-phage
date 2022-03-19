library(here)
library(tidyverse)
library(cowplot)
library(gtools)

# Import data -------------------------------------------------------------
# Collecting DE data for all loci in all treatments
# reading in the FDR corrected p-values for indicating significantly DEx'ed genes.
f <- 
  list.files(here("RNAseq","data","DEseq"), pattern = "txt",recursive = T, full.names = T)

f <- f[str_detect(f, "pDR110")]
d <- read_delim(f, delim = "\t")

# select controls
d <- 
  d %>% 
  select(id,
         colnames(d)[colnames(d) %>% str_detect("pDR110.*r[123]")])
      


# * arrange data -------------------------------------------------
# locus tags
# list matching delta6 locus_tags with 168:  

d6.spor.genes <- 
  read_csv(here("RNAseq/data/annotations/delta6_Spor_Annot.csv"))

#removing RNA genes
d6.spor.cds <- 
  d6.spor.genes %>%
  filter(str_detect(locus_tag.168,"BSU"))%>%
  filter(!str_detect(locus_tag.168,"RNA"))

# sporulation genes by subtiwiki
d6.spor.cds <- d6.spor.cds %>% 
  mutate(sw.spore  = str_detect(replace_na(category2," "), regex("sporulation", ignore_case = T)))


d.long <-
  # filter only cds
  right_join(d, d6.spor.cds, by=c("id"="locus_tag.d6")) %>% 
  pivot_longer(cols = (ends_with("_r1")|ends_with("_r2")|ends_with("_r3"))) %>% 
  separate(name, into = c("induced","iptg","rep")) %>% 
  rename(reads = value) %>% 
  #normalize to gene length
  mutate(normalized_reads = reads+1 / (end-start))


# plot --------------------------------------------------------------------

p <- d.long %>%
  rename(sporul.gene = sw.spore) %>% 
  ggplot(aes(log10(normalized_reads)))+
  geom_density(aes(fill = sporul.gene, color = sporul.gene), alpha = 0.6)+
  geom_density()+
  scale_fill_viridis_d(direction = -1)+
  scale_color_viridis_d(direction = -1)+
  facet_grid(rep ~ iptg)+
  theme_classic()+
  theme(legend.position = "bottom")+
  panel_border(color = "black")


ggsave(here("RNAseq/plots/control_reads.png"), plot = p,
       width = 6, height = 4)



# END ---------------------------------------------------------------------


