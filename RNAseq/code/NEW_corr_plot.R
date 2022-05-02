library(here)
library(tidyverse)
library(cowplot)
library(gtools)
library(scales)
library(broom)


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


# compare phage- to host- derived genes  ----------------------------------

sig.host.v <- c("sigF","sigG")
sig.phage.v <- c("ELDg168","ELDg169","Goe3","SP10")

#organize data for host-phage comparisons
d.hp <- d.all %>% 
  filter(!is.na(fc)) %>% 
  select(id, sw.spore, induced,fc) %>% 
  pivot_wider(names_from = induced, values_from = fc) %>% 
  pivot_longer(cols = sig.phage.v, names_to = "induced.phage", values_to = "fc.phage")  
  
#test correlation
d.cor.test <-
  bind_rows(
    
    d.hp %>% 
      group_by(induced.phage) %>% 
      summarise(
        cor.test(log2(sigF), log2(fc.phage), method = "spearman") %>%
          tidy() %>% mutate(induced.host = "sigF")),
    
    d.hp %>% 
      group_by(induced.phage) %>% 
      summarise(
        cor.test(log2(sigG), log2(fc.phage), method = "spearman") %>%
          tidy() %>% mutate(induced.host = "sigG"))
  )
# adjust correlation p value for multiple testing
d.cor.test <- d.cor.test %>% 
  mutate(adj.p.BH = p.adjust(p.value, method = "BH")) %>% 
  relocate(adj.p.BH, .after = "p.value") %>% 
  relocate(induced.host)

write_csv(d.cor.test, here("RNAseq/data/FC_spearman.csv"))





# Plot --------------------------------------------------------------------

# arrangements for plotting
brx <- c( -10,0,10)

d.cor.test <- 
  d.cor.test %>% 
  mutate(p.lab = paste('rho', "==",signif(estimate,3))) %>% 
  mutate(induced.phage = fct_relevel(induced.phage,"ELDg169","ELDg168","Goe3","SP10"))

# plot
p <- d.hp %>% 
  # arrangements for plotting
  mutate(induced.phage = fct_relevel(induced.phage,"ELDg169","ELDg168","Goe3","SP10"),
         sw.spore = if_else(sw.spore, "spor. genes", "other") %>% fct_rev()) %>% 
  #plot
  ggplot(aes(x=log2(sigF), y=log2(fc.phage)))+
  geom_hline(yintercept = 0, color="grey40")+
  geom_vline(xintercept = 0, color="grey40")+
  geom_point(aes(color = sw.spore),shape=20, alpha = 0.5, size = 0.8)+
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black")+
  # add cor test result
  geom_text(data = d.cor.test %>% filter(induced.host=="sigF"),
            aes(label = p.lab),  parse = TRUE,
            x = Inf, y = -Inf, hjust = 1.1, vjust = -0.2, lineheight = .5 )+
  scale_y_continuous(breaks = brx)+
  scale_x_continuous(breaks = brx)+

  xlab(bquote(log[2]~FC~(italic("sigF")))) +
  ylab(bquote(log[2]~FC~("phage gene"))) +
  
  facet_wrap(~induced.phage, ncol = 1)+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  scale_color_grey()+
  theme(legend.position="bottom",
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b=2)))+
  guides(color = guide_legend(title = NULL, override.aes = list(size=5, alpha = 1)))

ggsave(here("RNAseq/plots/correlations-sigF.png"),
       plot = p, width = 2, height = 6)



# sigG --------------------------------------------------------------------

# plot
p <- d.hp %>% 
  # arrangements for plotting
  mutate(induced.phage = fct_relevel(induced.phage,"ELDg169","ELDg168","Goe3","SP10"),
         sw.spore = if_else(sw.spore, "spor. genes", "other") %>% fct_rev()) %>% 
  #plot
  ggplot(aes(x=log2(sigF), y=log2(fc.phage)))+
  geom_hline(yintercept = 0, color="grey40")+
  geom_vline(xintercept = 0, color="grey40")+
  geom_point(aes(color = sw.spore),shape=20, alpha = 0.5, size = 0.8)+
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black")+
  # add cor test result
  geom_text(data = d.cor.test %>% filter(induced.host=="sigG"),
            aes(label = p.lab),  parse = TRUE,
            x = Inf, y = -Inf, hjust = 1.1, vjust = -0.2, lineheight = .5 )+
  scale_y_continuous(breaks = brx)+
  scale_x_continuous(breaks = brx)+
  
  xlab(bquote(log[2]~FC~(italic("sigG")))) +
  ylab(bquote(log[2]~FC~("phage gene"))) +
  
  facet_wrap(~induced.phage, ncol = 1)+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  scale_color_grey()+
  theme(legend.position="bottom",
        legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_text(margin = margin(b=2)))+
  guides(color = guide_legend(title = NULL, override.aes = list(size=5, alpha = 1)))

ggsave(here("RNAseq/plots/correlations-sigG.png"),
       plot = p, width = 2, height = 6)

