# summarise by viral family
tmp <- d.sp%>%
  # phylum name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  # extract viral family
  mutate(viral.family=str_extract(`virus lineage`,
                                  regex("caudovirales;.*", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex("caudovirales; ", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex(";.*", ignore_case = T)))%>%
  mutate(viral.family = if_else(is.na(viral.family), "non-Caudovirales*", viral.family)) %>% 
  mutate(viral.family = fct_relevel(viral.family, "non-Caudovirales*", after = 0)) %>% 
filter(! is.na(viral.family)) %>%
  filter(! is.na(phylum))
  # group_by(viral.family, phylum)%>%
  # mutate(has_sigma = n.sigma >0,
  #        multi_sigma = n.sigma >1) %>% 
  # summarise(n=n(),
  #           w.sigma=sum(has_sigma),
  #           perc_sigma=100*sum(has_sigma)/n(),
  #           w.multi=sum(multi_sigma),
  #           perc_multi=100*sum(multi_sigma)/n())


p.phylum <-  tmp %>% 
  group_by(phylum)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>% 
  summarise(n=n(),
            w.sigma=sum(has_sigma),
            perc_sigma=100*sum(has_sigma)/n(),
            w.multi=sum(multi_sigma),
            perc_multi=100*sum(multi_sigma)/n()) %>% 
  
  filter(n>1) %>% 
  ggplot(aes(phylum))+
  geom_col(aes(y=n), position=position_dodge(preserve = "single"), alpha=0.5, color="grey80")+
  # geom_col(aes(y=w.multi, fill=viral.family), position=position_dodge(preserve = "single"))+
  geom_col(aes(y=w.sigma), position=position_dodge(preserve = "single"), alpha=1, color="black")+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,1000))+
  xlab("Host Phylum")+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  coord_flip(expand = F)+
  scale_fill_viridis_d()+
  scale_x_discrete()

p.phylum

p.Vfam <-  tmp %>% 
  group_by(viral.family)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>% 
  summarise(n=n(),
            w.sigma=sum(has_sigma),
            perc_sigma=100*sum(has_sigma)/n(),
            w.multi=sum(multi_sigma),
            perc_multi=100*sum(multi_sigma)/n()) %>% 
  
  filter(n>1) %>% 
  ggplot(aes(viral.family))+
  geom_col(aes(y=n), position=position_dodge(preserve = "single"), alpha=0.5, color="grey80")+
  # geom_col(aes(y=w.multi, fill=viral.family), position=position_dodge(preserve = "single"))+
  geom_col(aes(y=w.sigma), position=position_dodge(preserve = "single"), alpha=1, color="black")+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,2000))+
  xlab("Viral Family")+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  coord_flip(expand = F)+
  scale_fill_viridis_d()

p.Vfam

p.both <-  tmp %>% 
  group_by(viral.family, phylum)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>% 
  summarise(n=n(),
            w.sigma=sum(has_sigma),
            perc_sigma=100*sum(has_sigma)/n(),
            w.multi=sum(multi_sigma),
            perc_multi=100*sum(multi_sigma)/n()) %>% 

  filter(n>1) %>% 
  mutate(viral.family = fct_infreq(viral.family) ) %>% 
  ggplot(aes(viral.family))+
  geom_col(aes(y=n, fill=viral.family), position=position_dodge(preserve = "single"), alpha=0.5, color="grey80")+
  # geom_col(aes(y=w.multi, fill=viral.family), position=position_dodge(preserve = "single"))+
  geom_col(aes(y=w.sigma, fill=viral.family), position=position_dodge(preserve = "single"), alpha=1, color="black")+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,1000))+
  xlab("Viral Family")+
  theme_classic(base_size = 12)+
  panel_border(color = "black")+
  facet_wrap(~phylum, nrow = 2)+
  coord_cartesian(expand = F)+
  scale_fill_viridis_d()+
  scale_x_discrete(limits=rev)+
  theme(axis.text.x = element_blank(),
        legend.position = c(0.9,0.25),
        legend.text = element_text(size = 8),
        legend.key.height = unit(1,"mm"))+
  guides(fill =  guide_legend(title = "Viral Family", 
                              reverse = TRUE,
                              override.aes = list(alpha = 1)))


p.both

top_row <- plot_grid(p.phylum, NULL, p.Vfam, nrow = 1,
                     labels = c("a","","b"),rel_widths = c(1,0.1,1))
all.3 <- plot_grid(top_row, p.both, ncol = 1, labels = c("","c"), 
                   rel_heights = c(1,1.5))



ggsave2(here("vogdb","figures","viral_family_host_phylum.png"),
        plot = ggdraw(all.3) +
               theme(plot.background = element_rect(fill="white", color = NA)),
        width = 8,height = 6)



# t1 <- aov(perc_sigma ~ viral.family * phylum, tmp)
# summary.aov(t1)
# 
# t2 <- aov(perc_multi ~ viral.family * phylum, tmp)
# summary.aov(t2)



# tmp <- d.sp %>% filter(! str_detect(d.sp$`virus lineage`,regex("caudovirales", ignore_case = T)))
