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



# END ---------------------------------------------------------------------


