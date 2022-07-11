library(here)
library(tidyverse)
library(cowplot)
library(gtools)

# Import data -------------------------------------------------------------
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

# sporulation gene enrich ---------------------
# Are sporulation genes enriched in the sample of deferentially expressed genes?
# the FDR corrected p-values for indicating significantly DExd genes.

# sporulation genes by subtiwiki
d6.spor.cds <- d6.spor.cds %>% 
  mutate(sw.spore  = str_detect(replace_na(category2," "), regex("sporulation", ignore_case = T)))


# define up/sown regulated
d.all <- d.all %>%
  mutate(up.down = case_when(
    (p<0.05) & (fc > 2) ~ "up", 
    (p<0.05) & (fc < 0.5)~ "down",
    TRUE ~ "unchanged"))

# record p values
p.val.spore <- tibble()

# summarize by induced gene
p.val.spore <-
  d.all %>% 
  group_by(induced, up.down, sw.spore) %>% 
  summarise(n = n()) %>% 
  # mutate(total = sum(n)) %>% 
  ungroup() %>% 
  mutate(sw.spore = if_else(sw.spore, "spor", "other"),
         gene_cat = paste(sw.spore, up.down, sep = ".")) %>% 
  select(-up.down, -sw.spore) %>% 
  pivot_wider(names_from = gene_cat, values_from = n, values_fill = 0) %>% 
  # define hg parameters
  mutate(
    # total sporulation genes
    m = spor.down + spor.up + spor.unchanged,
    # non-sporulation genes
    n = other.down + other.up + other.unchanged,
    # total number of upregulated genes
    k.up = spor.up + other.up,
    # total number of downregulated genes
    k.down = spor.down + other.down
  ) %>%
  # hg test
  mutate(
    hg.up = signif(phyper(spor.up, m, n, k.up, lower.tail =F),5),
    hg.down = signif(phyper(spor.down, m, n, k.down, lower.tail =F),5) 
  )

# Adjust Pvalue for multiple testing
p.val.spore <- 
  p.val.spore %>% 
  # get single column of al P values
  pivot_longer(cols = starts_with("hg."), 
               names_to = "direction",
               values_to = "p") %>% 
    mutate(direction = str_replace(direction, "hg", "Padj")) %>% 
  # P adjustment
    mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  # add back to main table of stats
    select(induced, direction, p.adj) %>% 
    pivot_wider(names_from = direction, values_from = p.adj) %>% 
    left_join(p.val.spore, ., by = "induced")

# export analysis test results
# write_csv(p.val.spore,
#           here("RNAseq/data/sporulation_gene_enrichment.csv"))

# format for sul. data
p.down <- select(p.val.spore, -ends_with("up")) %>% mutate(direction = "down")
p.up <- select(p.val.spore, -ends_with("down")) %>% mutate(direction = "up")

colnames(p.up) <-  str_replace (colnames(p.up), ".up", "\\.DExed")
colnames(p.down) <-  str_replace (colnames(p.down), ".down", "\\.DExed")

p.sup <- bind_rows(p.up, p.down) %>%
  arrange(induced) %>%
  relocate(direction, .after = 1)

write_csv(p.sup,
          here("RNAseq/data/sporulation_gene_enrichment.csv"))

# Plotting --------------------------------------------------------------

# gene count summary and stats
gene_dexed <- 
  p.val.spore %>% 
  select(induced, starts_with("spor"),starts_with("other"), starts_with("Padj")) %>% 
  pivot_longer(cols = c(starts_with("spor"),starts_with("other")), 
               names_to = "g_cat", values_to = "n") %>% 
  separate(g_cat, into = c("sw.spore", "up.down"), sep = "\\.") %>% 
  # label for plot with n and adjusted P-value
  mutate(n_lab = case_when(
    (sw.spore == "spor" & 
       up.down == "up") ~ paste(n, stars.pval(Padj.up), sep = " "),
    (sw.spore == "spor" & 
       up.down == "down") ~ paste(n, stars.pval(Padj.down), sep = " "),
    TRUE ~ as.character(n)
  )) %>% 

# adjustments for plotting
  mutate(pnl=case_when(induced %in% c("sigF","sigG") ~ "host",
                       TRUE ~ "phage") )  %>% 
  mutate(strip = paste0(pnl,": ",induced)) %>% 
  mutate(sw.spore = str_replace(sw.spore, "spor", "spor. genes") %>% 
           as_factor() %>% fct_rev()) %>% 
  mutate(up.down = as_factor(up.down) %>% fct_relevel("up","down")) %>% 
  mutate(x_pos = case_when(up.down=="up"~ 7,
                           up.down=="down"~ -7,
                           TRUE ~ 0)) %>% 
  mutate(y_pos = case_when(up.down=="up"~ 275,
                           up.down=="down"~ 275,
                           TRUE ~ 215)) %>% 
  mutate(induced = fct_relevel(induced, "ELDg169"))

# prep for plot
p <-  d.all %>%
  mutate(up.down = case_when(
    (p<0.05) & (fc > 2) ~ "up", 
    (p<0.05) & (fc < 0.5)~ "down",
    TRUE ~ "unchanged") %>% 
      as_factor() %>% fct_relevel("up","down")) %>%
  mutate(pnl=case_when(induced %in% c("sigF","sigG") ~ "host",
                       TRUE ~ "phage") ) %>%  
  mutate(induced = fct_relevel(induced, "ELDg169")) %>%
  mutate(sw.spore = if_else(sw.spore, "spor. genes", "other") %>% 
           as_factor() %>% fct_rev()) %>% 
  # plot
  ggplot(aes(log2(fc), -log10(p)))+
  geom_rect(xmin=-log2(2), xmax=log2(2), ymin=-Inf, ymax=210,
            fill = "transparent", color = "grey90", alpha = 0.5)+
  geom_hline(yintercept = -log10(0.05), color = "grey90", linetype = 2)+
  geom_point(aes(color = up.down), size=.5, alpha = 0.5, shape = 20)+
  geom_label(data = gene_dexed, 
            aes(label = n_lab, color = up.down, x= x_pos, y = y_pos), 
            vjust= 0.5, hjust = 0.5, size = 2.5,  show.legend = F,
            fill = "white", label.size = NA)+
  facet_grid(sw.spore~pnl + induced)+
  scale_colour_viridis_d(direction = 1)+
  theme_classic()+
  panel_border(color = "black")+
  xlab(expression(log[2]~FC)) + 
  ylab(expression(-log[10]~P~value)) + 
  ylim(0,280)+
  # scale_x_continuous(breaks = c(-10,-1,1,10), labels = c(-10,-2,2,10))+
  theme(legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size = 12),
        strip.background = element_blank())+
  labs(color = "expression change")+
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 4, alpha = 1)))

ggsave(here("RNAseq/plots/volcano_plot.png"), plot = p,
       width = 6, height = 4)

# save(p, file = here("RNAseq/plots/volcano-plot.Rdata"))


# END ---------------------------------------------------------------------


