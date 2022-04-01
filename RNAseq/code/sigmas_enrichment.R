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


# sigH regulon -----------------------------------
d.reg <- 
  read_csv(here("RNAseq/data/annotations/SW_regulations.csv"))

sigmas <- d.reg %>% 
  filter(mode == "sigma factor") %>% 
  group_by(regulon, regulator ) %>% 
  summarise(n = n()) %>% 
  filter(str_detect(regulon, "Sig|YlaC")) %>% 
  pull(regulator)

all.pvals <- tibble()


for(sig in sigmas){
  
  if(sig %in% c("SigL", "YlaC")) next
  # the above two caus error: "object 'sig.down' not found"
  
  d.sig <- 
    d.reg %>% 
    filter(str_detect(regulon, regex(sig, ignore_case = T)))
  
  # add delat6 locus tags
  d.sig <- 
    d6.spor.cds %>% 
    select(locus_tag.168, locus_tag.d6) %>% 
    right_join(., d.sig, by = c(locus_tag.168 = "locus tag")) %>% 
    filter(!is.na(locus_tag.d6)) 
  
  d.all <-  
    full_join(fc, p.val, by = "id", suffix = c("_fc", "_p")) %>%
    # filter only cds
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6")) %>%
    # mark sig regulated genes
    mutate(sig.reg = id %in% d.sig$locus_tag.d6 ) %>%
    pivot_longer(cols = (ends_with("_fc") | ends_with("_p"))) %>% 
    separate(name, into = c("induced","val")) %>% 
    pivot_wider(names_from = "val", values_from = "value" ) %>% 
    # change genes that were not assigned a p-value for DE to 1 (no change)
    mutate(p=if_else(is.na(p),1,p))%>%
    filter(induced != "pDR110") 
  
  
  
  
  # sig-regulated gene enrich ---------------------
  # Are sig-regulated genes enriched in the sample of deferentially expressed genes?
  # the FDR corrected p-values for indicating significantly DExd genes.
  
  
  # define up/sown regulated
  d.all <- d.all %>%
    mutate(up.down = case_when(
      (p<0.05) & (fc > 2) ~ "up", 
      (p<0.05) & (fc < 0.5)~ "down",
      TRUE ~ "unchanged"))
  
  # record p values
  p.val.sig <- tibble()
  
  # summarize by induced gene
  p.val.sig <-
    d.all %>% 
    group_by(induced, up.down, sig.reg) %>% 
    summarise(n = n()) %>% 
    # mutate(total = sum(n)) %>% 
    ungroup() %>% 
    mutate(sig.reg = if_else(sig.reg, "sig", "other"),
           gene_cat = paste(sig.reg, up.down, sep = ".")) %>% 
    select(-up.down, -sig.reg) %>% 
    pivot_wider(names_from = gene_cat, values_from = n, values_fill = 0) %>% 
    # define hg parameters
    mutate(
      # total sig genes
      m = sig.down + sig.up + sig.unchanged,
      # non-sig genes
      n = other.down + other.up + other.unchanged,
      # total number of upregulated genes
      k.up = sig.up + other.up,
      # total number of downregulated genes
      k.down = sig.down + other.down
    ) %>%
    # hg test
    mutate(
      hg.up = signif(phyper(sig.up, m, n, k.up, lower.tail =F),5),
      hg.down = signif(phyper(sig.down, m, n, k.down, lower.tail =F),5) 
    )
  
  all.pvals <- 
    p.val.sig %>% 
    mutate(cur_sig = sig) %>% 
    relocate(cur_sig) %>% 
    bind_rows(all.pvals, .)
}





# ------------------------------------------


# Adjust Pvalue for multiple testing
all.pvals <-
  all.pvals %>% 
  # get single column of al P values
  pivot_longer(cols = starts_with("hg."), 
               names_to = "direction",
               values_to = "p") %>% 
    mutate(direction = str_replace(direction, "hg", "Padj")) %>% 
  # P adjustment
    mutate(p.adj = p.adjust(p, method = "BH")) %>% 
  # add back to main table of stats
    select(cur_sig, induced, direction, p.adj) %>% 
    pivot_wider(names_from = direction, values_from = p.adj) %>% 
    left_join(all.pvals, .,  by = c("cur_sig" ,"induced"))

# # export analysis test results
# write_csv(p.val.spore,
#           here("RNAseq/data/sigH_gene_enrichment.csv"))



# plot --------------------------------------------------------------------

p <- all.pvals %>% 
  
  pivot_longer(cols = c(starts_with("Padj")),
               names_to = "direction", values_to = "p") %>% 
  mutate(astr = if_else(p<=0.05, stars.pval(p), ""),
         direction = if_else(str_detect(direction,"up"),
                             "upregulaed", "downregulated") %>% fct_rev()) %>%
  mutate(sig.grp = case_when(
    str_detect(cur_sig, "E|F|G|K") ~ "sporulation",
    str_detect(cur_sig, "A|B|D|H") ~ "vegetative",
    str_detect(cur_sig, "M|W|X|V|Y") ~ "ECF",
    TRUE ~ "other"
  ) %>% factor(levels = c("sporulation","vegetative", "ECF", "other" ))) %>% 
  ggplot(aes(y=induced, x=cur_sig))+
  geom_tile(aes(fill = -log10(p)), color = "black")+
  geom_text(aes(label = astr), size = 5, vjust=0.8)+
  facet_grid(direction ~ sig.grp, scales = "free_x", space = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  xlab("regulon tested")+
  theme_classic()+
  panel_border(color = "black")


ggsave(here("RNAseq/plots/sigmas_enrichment.png"), plot = p,
       width = 8, height = 4)


# plot 2 ------------------------------------------------------------------

p <- all.pvals %>% 
  mutate(frac_up = sig.up/m,
         frac_down = sig.down/m) %>% 
  pivot_longer(cols = c(starts_with("frac")),
               names_to = "direction", values_to = "frac") %>% 
  # mutate(astr = if_else(p<=0.05, stars.pval(p), ""),
  #        direction = if_else(str_detect(direction,"up"),
  #                            "upregulaed", "downregulated") %>% fct_rev()) %>%
  mutate(sig.grp = case_when(
    str_detect(cur_sig, "E|F|G|K") ~ "sporulation",
    str_detect(cur_sig, "A|B|D|H") ~ "vegetative",
    str_detect(cur_sig, "M|W|X|V|Y") ~ "ECF",
    TRUE ~ "other"
  ) %>% factor(levels = c("sporulation","vegetative", "ECF", "other" ))) %>% 
  ggplot(aes(y=induced, x=cur_sig))+
  geom_tile(aes(fill = frac), color = "black")+
  # geom_text(aes(label = astr), size = 5, vjust=0.8)+
  facet_grid(direction ~ sig.grp, scales = "free_x", space = "free")+
  scale_fill_gradient(low = "white", high = "red")+
  xlab("regulon tested")+
  theme_classic()+
  panel_border(color = "black")

p
# ggsave(here("RNAseq/plots/sigmas_enrichment.png"), plot = p,
#        width = 8, height = 4)



# plot 3 ------------------------------------------------------------------
d.fc <- tibble()
all.fc <- tibble()
for(sig in sigmas){
  
  # if(sig %in% c("SigL", "YlaC")) next
  # the above two caus error: "object 'sig.down' not found"
  
  d.sig <- 
    d.reg %>% 
    filter(str_detect(regulon, regex(sig, ignore_case = T)))
  
  # add delat6 locus tags
  d.sig <- 
    d6.spor.cds %>% 
    select(locus_tag.168, locus_tag.d6) %>% 
    right_join(., d.sig, by = c(locus_tag.168 = "locus tag")) %>% 
    filter(!is.na(locus_tag.d6)) 
  
  d.all <-  
    full_join(fc, p.val, by = "id", suffix = c("_fc", "_p")) %>%
    # filter only cds
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6")) %>%
    # mark sig regulated genes
    mutate(sig.reg = id %in% d.sig$locus_tag.d6 ) %>%
    pivot_longer(cols = (ends_with("_fc") | ends_with("_p"))) %>% 
    separate(name, into = c("induced","val")) %>% 
    pivot_wider(names_from = "val", values_from = "value" ) %>% 
    # change genes that were not assigned a p-value for DE to 1 (no change)
    mutate(p=if_else(is.na(p),1,p))
    # filter(induced != "pDR110") 
  
  d.fc <- d.all %>% 
    
    filter(sig.reg) %>% 
    mutate(fc.log2 = log2(fc)) %>% 
    group_by(induced) %>% 
    summarise(m.fc = mean(fc.log2, na.rm = T), 
              sd.fc = sd(fc.log2, na.rm = T), 
              n=n()) %>% 
    mutate(regulon = sig) %>% 
      relocate(regulon) %>% 
      bind_rows(d.fc,.)
  
  all.fc <-
    d.all %>% 
    filter(sig.reg) %>% 
    mutate(regulon = sig) %>% 
    select(induced, regulon, id, fc,p ) %>% 
    
    bind_rows(all.fc, .)
  
  
}

# diverging colors 
cl <-  viridisLite::viridis(2)
cl <- c(cl[2], "#FFFFFF", cl[1])

p <- d.fc %>% 
  mutate(sig.grp = case_when(
    str_detect(regulon, "E|F|G|K") ~ "sporulation",
    str_detect(regulon, "A|B|D|H") ~ "vegetative",
    str_detect(regulon, "M|W|X|V|Y") ~ "ECF",
    TRUE ~ "other"
  ) %>% factor(levels = c("sporulation","vegetative", "ECF", "other" ))) %>% 
  ggplot(aes(y=induced, x=regulon))+
  geom_tile(aes(fill = m.fc), color = "black")+
  facet_grid(. ~ sig.grp, scales = "free_x", space = "free")+
  # scale_fill_gradient2()+
  # scale_fill_gradientn(colours = colorspace::diverge_hcl(palette = "Blue-Red2",3),
  #                      values = c(0,0.21,1) )+
  scale_fill_gradientn(colours =  cl,
                       values = c(0,0.21,1) )+
  xlab("regulon tested")+
  theme_classic()+
  panel_border(color = "black")

p

p <- d.fc %>% 
  mutate(sig.grp = case_when(
    str_detect(regulon, "E|F|G|K") ~ "sporulation",
    str_detect(regulon, "A|B|D|H") ~ "vegetative",
    str_detect(regulon, "M|W|X|V|Y") ~ "ECF",
    TRUE ~ "other"
  ) %>% factor(levels = c("sporulation","vegetative", "ECF", "other" ))) %>% 
  filter(sig.grp %in% c("sporulation","vegetative")) %>% 
  ggplot(aes(y=m.fc, x=regulon))+
  geom_hline(yintercept = 0, size=0.2)+
  geom_col(fill = "grey", color = "black", position = "dodge", width = 0.5)+
  geom_errorbar(aes(ymin = m.fc-sd.fc, ymax = m.fc+sd.fc), width = 0.2)+
  facet_grid(induced ~ sig.grp, scales = "free_x", space = "free")+
  # xlab("regulon tested")+
  theme_classic()+
  panel_border(color = "black")

p
  

ggsave(here("RNAseq/plots/sigmas_fc.png"), plot = p,
       width = 6, height = 6)

# ppp ---------------------------------------------------------------------


p <- all.fc %>% 
  mutate(sig.grp = case_when(
    str_detect(regulon, "E|F|G|K") ~ "sporulation",
    str_detect(regulon, "A|B|D|H") ~ "vegetative",
    str_detect(regulon, "M|W|X|V|Y") ~ "ECF",
    TRUE ~ "other"
  ) %>% 
    factor(levels = c("sporulation","vegetative", "ECF", "other" ))) %>% 
  filter(sig.grp %in% c("sporulation","vegetative")) %>% 
  mutate(pnl=case_when(induced %in% c("sigF","sigG") ~ "host",
                       induced=="pDR110" ~ "ctrl.",
                       TRUE ~ "phage") ) %>% 
  mutate(induced = fct_relevel(induced, "ELDg169")) %>%
  ggplot(aes(y=log2(fc), x=regulon))+
  geom_hline(yintercept = 0, size=0.2)+
  geom_violin(fill = "grey", color = "black", width=0.5)+
  # geom_col(fill = "grey", color = "black", position = "dodge", width = 0.5)+
  # geom_errorbar(aes(ymin = m.fc-sd.fc, ymax = m.fc+sd.fc), width = 0.2)+
  facet_grid(pnl + induced ~ sig.grp, scales = "free_x", space = "free")+
  # xlab("regulon tested")+
  theme_classic()+
  panel_border(color = "black")

p


ggsave(here("RNAseq/plots/sigmas_fc.png"), plot = p,
       width = 6, height = 6)
# END ---------------------------------------------------------------------


