library(here)
library(tidyverse)
library(cowplot)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)


# import PID ------------------------------------------------------------
# > parse results -----------------------------------------------------------


pid_sigBFG <- read.table(here("synthesis/data/hmm_align", "pid_sigBFG.tsv")) %>% 
  as_tibble()
colnames(pid_sigBFG) <- c("seqname1","seqname2","perc_id","nid","denomid","perc_match","nmatch","denommatch")

# names of proteins
protein.names <- read_csv(here("synthesis/data/hmm_align","seq_names.csv"))

# replace names
pid_sigBFG <- left_join(pid_sigBFG, protein.names, by = c("seqname1" = "protein")) %>% 
  mutate(seqname1 = name) %>% 
  rename(phage1 = phage) %>% 
  select(-name) %>% 
  left_join(., protein.names, by = c("seqname2" = "protein")) %>% 
  mutate(seqname2 = name) %>% 
  rename(phage2 = phage) %>% 
  select(-name)


#organize data: phage vs host

# remove phage-phage comparisons
pid_sigBFG <- pid_sigBFG %>% 
  filter(str_detect(seqname1, "sigF"))


# sporulation vs PID -------------------------------------------------

# import sporulation data
d.spor <- read_csv(here("synthesis/data/FCM_sporulation_induction.csv")) %>% 
  filter(strain != "Empty Vector")

# summarise 
d.spor <- d.spor %>% 
  group_by(strain) %>% 
  summarise(yield = mean(induction.ratio),n=n(), v = sd(induction.ratio)/sqrt(n))


d.both <- pid_sigBFG %>% 
  rename(strain = seqname2) %>% 
  mutate(strain = str_remove(strain, "Bsub_")) %>% 
  left_join(d.spor, . , by = "strain") %>% 
  mutate(perc_id = if_else(strain == "sigF", 100, perc_id))



# linear model
# linearMod <- lm(yield ~ perc_id, data = d.both)
# sum.lm <- summary(linearMod)
# rsq <- sum.lm$adj.r.squared %>% signif(2)
# pval <- sum.lm$coefficients["perc_id","Pr(>|t|)"]%>% signif(2)

# get r for plot
r <- cor(d.both$perc_id, d.both$yield, method = "spearman") %>% 
  signif(3) 
p <- cor.test(d.both$perc_id, d.both$yield, 
              method = "spearman")$p.value %>% 
  signif(3)

p1 <- d.both %>% 
  ggplot(aes(perc_id, yield)) +
  # geom_abline(intercept = linearMod$coefficients[1],
  #             slope = linearMod$coefficients[2],
  #             linetype = 2)+
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black", 
              se = F, linetype = 2)+
  geom_errorbar(aes(ymin  = yield-v, ymax = yield + v))+
  geom_point(aes(fill  = strain),shape=21, size = 2)+
  geom_text(label = paste('rho', "==",r), parse = TRUE,
            x=Inf, y = Inf, hjust = 1.1, vjust =1.5)+
  geom_text(label = paste("\nP =",p), parse = FALSE,
           x=Inf, y = Inf, hjust = 1.1, vjust =1.5)+
  theme_classic()+
  panel_border(color = "black")+
  scale_shape_manual(values = 21:25)+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.25))+
  scale_fill_viridis_d()+
  xlab("Percent ID (compared to sigF)")+
  ylab("Spore yield")+
  guides(fill = guide_legend(title = "cloned\nsigma"))




# PID vs DEG spor gene ----------------------------------------------------

deg <- read_csv("RNAseq/data/sporulation_gene_enrichment.csv")

d.both2 <- pid_sigBFG %>% 
  rename(strain = seqname2) %>% 
  mutate(strain = str_remove(strain, "Bsub_")) %>% 
  right_join(.,deg %>% filter(direction=="up"),
             by = c("strain" = "induced")) %>% 
  mutate(perc_id = if_else(strain == "sigF", 100, perc_id))

# # linear model
# linearMod <- lm(spor.DEG ~ perc_id, data = d.both2)
# sum.lm <- summary(linearMod)
# rsq <- sum.lm$adj.r.squared %>% signif(2)
# pval <- sum.lm$coefficients["perc_id","Pr(>|t|)"]%>% signif(2)

# get r for plot
r <- cor(d.both2$perc_id, d.both2$spor.DExed, method = "spearman") %>% 
  signif(3)
p <- cor.test(d.both2$perc_id, d.both2$spor.DExed, 
              method = "spearman")$p.value %>% 
  signif(3)

p2 <- d.both2 %>% 
  ggplot(aes(perc_id, spor.DExed)) +
  geom_smooth(method = 'lm', formula = 'y ~ x', color = "black", 
              se = F, linetype = 2)+
  # geom_abline(intercept = linearMod$coefficients[1],
  #             slope = linearMod$coefficients[2],
  #             linetype = 2)+
  geom_point(aes(fill  = strain),shape=21, size = 2)+
  geom_text(label = paste('rho', "==",r), parse = TRUE, #paste0("R2=", rsq, "\nP = ",pval ),
            x=-Inf, hjust = -0.1, y = Inf, vjust = 1.5)+
  geom_text(label = paste("\nP =",p), parse = FALSE,
            x=-Inf, hjust = -0.1, y = Inf, vjust = 1.5)+
  theme_classic()+
  panel_border(color = "black")+
  scale_shape_manual(values = 21:25)+
  scale_x_continuous( breaks = seq(20,100,20), limits = c(20,100))+
  scale_fill_viridis_d()+
  xlab("Percent ID (compared to sigF)")+
  ylab("Upregulated sporulation genes")+
  guides(fill = guide_legend(title = "cloned\nsigma"))


save_plot(here("synthesis/plots","ID_yield_DEG.png"),
          plot_grid(p1+theme(legend.position = "none"),
                    NULL, p2,
                    rel_widths = c(1,0.1,1.4),
                    labels = c("c","","d"), nrow = 1),
          base_height = 3, base_width = 7)
