library(here)
library(tidyverse)
library(cowplot)
library(gtools)


# Collecting DE data for all loci in all treatments
# reading in the FDR corrected p-values for indicating significantly DEx'ed genes.
f <- 
  list.files(here("RNAseq","data","DEseq"), pattern = "txt",recursive = T, full.names = T)
d <- read_delim(f[1], delim = "\t")
d <-   select(d,id)
d2 <- d
for (file in f){
  
  de <- read_delim(file, delim = "\t")
  d <- full_join(d, select(de, id, contains("IHW")))
  d2<- full_join(d2, select(de, id, contains("Fold_Change")))
}

d <- d%>%
  rename_at(vars(-id), ~str_remove(.,"_.*"))
p.val <- d
d2 <- d2%>%
  rename_at(vars(-id), ~str_remove(.,"_.*"))
#fold change
fc <- d2
#clean up
rm(de,d,d2)


#------------------
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


d.all <- 
  full_join(fc, p.val, by = "id", suffix = c("_fc", "_p")) %>% 
  # filter only cds
  right_join(., d6.spor.cds, by=c("id"="locus_tag.d6")) %>% 
  pivot_longer(cols = (ends_with("_fc") | ends_with("_p"))) %>% 
  separate(name, into = c("induced","val")) %>% 
  pivot_wider(names_from = "val", values_from = "value" ) %>% 
  # change genes that were not assigned a p-value for DE to 1 (no change)
  mutate(p=if_else(is.na(p),1,p))%>%
  filter(induced != "pDR110")

# label 20 most significant genes of each treatment
gene_labs <- d.all %>%
  # mutate(induced = paste0("P[IPTG]-", induced)) %>% 
  filter(p<0.05) %>% 
  group_by(induced) %>%
  slice_max( n = 10, order_by =  fc)

# sum Dexed
gene_dexed <- d.all %>%
  # define logical vector of DE based on p-value
  mutate(up.down = case_when(
    (p<0.05) & (fc > 2) ~ "up", 
    (p<0.05) & (fc < 0.5)~ "down",
    TRUE ~ "unchaged")) %>%
  group_by(induced, sw.spore, up.down) %>%
  summarise(n_changed = n()) %>%
  mutate(pnl=case_when(induced %in% c("sigF","sigG") ~ "host",
                       TRUE ~ "phage") )  %>% 
  mutate(strip = paste0(pnl,": ",induced)) %>% 
mutate(strip = fct_relevel(strip, "phage: ELDg169", after = 2)) %>% 
  mutate(sw.spore = if_else(sw.spore, "spor. genes", "other") %>% 
           as_factor() %>% fct_rev()) %>% 
  mutate(x_pos = case_when(up.down=="up"~ 7,
                           up.down=="down"~ -7,
                           TRUE ~ 0)) %>% 
  mutate(y_pos = case_when(up.down=="up"~ 200,
                           up.down=="down"~ 200,
                           TRUE ~ 150))


# # triangle characters
# c.up <- sprintf("\u25B3")
# c.down <- sprintf("\u25BD")

p <-  d.all %>%
  mutate(up.down = case_when(
    (p<0.05) & (fc > 2) ~ "up", 
    (p<0.05) & (fc < 0.5)~ "down",
    TRUE ~ "unchaged") %>% 
      as_factor() %>% fct_relevel("up","down")) %>%
  mutate(pnl=case_when(induced %in% c("sigF","sigG") ~ "host",
                       TRUE ~ "phage") ) %>%  
  mutate(strip = paste0(pnl,": ",induced)) %>% 
  mutate(strip = fct_relevel(strip, "phage: ELDg169", after = 2)) %>%
  mutate(sw.spore = if_else(sw.spore, "spor. genes", "other") %>% 
           as_factor() %>% fct_rev()) %>% 
  
  ggplot(aes(log2(fc), -log10(p)))+
  geom_rect(xmin=-log2(2), xmax=log2(2), ymin=-Inf, ymax=50,
            fill = "grey90", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = 2)+
  # geom_text(data = gene_dexed, aes(label = paste0("\n ",c.up,up,"\n ",c.down ,down)), 
            # x= -Inf, y=Inf, vjust= 1, hjust = 0, size = 3)+
  geom_point(aes(color = up.down), size=.5, alpha = 0.5, shape = 20)+
  geom_text(data = gene_dexed, 
            aes(label = n_changed, color = up.down, x= x_pos, y = y_pos), 
            vjust= 0.5, hjust = 0.5, size = 2.5,  show.legend = F)+
  # ggrepel::geom_text_repel(data = gene_labs, aes(label = gene), max.overlaps = 20)+
  facet_grid(sw.spore~strip)+
  scale_colour_viridis_d(direction = 1)+
  theme_classic()+
  panel_border(color = "black")+
  xlab(expression(log[2]~FC)) + 
  ylab(expression(-log[10]~P~value)) + 
  theme(legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size = 12))+
  labs(color = "expression change")+
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4, alpha = 1)))

ggsave(here("RNAseq/plots/volcano_plot2.png"), plot = p,
       width = 7, height = 4)

save(p, file = here("RNAseq/plots/volcano-plot.Rdata"))


# END ---------------------------------------------------------------------


