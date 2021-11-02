library(tidyverse)
library(here)
library(cowplot)
library(ggh4x)
library(viridisLite)

#load data frame of phage sigma factors generated in (B)
load(file=here("vogdb","data","vog_sigma_clean.RData"))

d.faa <- d.faa%>%
  mutate(taxon=as.numeric(taxon))


# import all viruses used to assemble VOGs
vog.sp <-  read_tsv(here("vogdb","vogdb_downloads","vog.species.list"), trim_ws = T) %>% 
  as_tibble( .name_repair = "universal")

#import virus-host data
# downloaded from ftp://ftp.genome.jp/pub/db/virushostdb/ (24/Nov/2020)
vh.db <- read_tsv(here("vogdb","vogdb_downloads","virushostdb.tsv"))


# virus duplicates in VHDB data -------------------------------------------
# These reflect multiple hosts
duplicated(vh.db$`virus tax id`)%>%sum() # 3461
vh.db%>%
  filter(str_detect(`host lineage`,regex("bacteria",ignore_case = T)))%>%
  group_by(`virus tax id`, `virus name`)%>%
  summarise(n=n())%>%
  ggplot(aes(x=n))+
  geom_histogram()+
  scale_y_log10()+
  ggtitle("VHDB host number for phages" )


# match hosts to vog spp
d.sp <- left_join(vog.sp, vh.db, by = c("tax.id" = "virus tax id")) %>% 
  # focus on viruses of bacteria
  filter(str_detect(`host lineage`, "Bacteria"))


# Here I will break up the host lineage data so I verify that hosts are of similar taxonomy
# the host lineage is not uniform having mostly 4-9 levels
# str_count(d.sp$`host lineage`,";")%>%range(na.rm = T)
# str_count(d.sp$`host lineage`,";")%>%hist()
# An example looks like this: 
# Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales
# Which corresponds to 
# Domain; Phylum; class; order. The next level down would be family

#but some have an intermediate rank between domain and phylum
# specifically FCB group AND Terrabacteria group 
# I will remove those for consistency

d.sp <- d.sp%>%
  mutate(`host lineage`=str_remove(`host lineage`," FCB group;"))%>%
  mutate(`host lineage`=str_remove(`host lineage`," Terrabacteria group;"))

d.sp <- d.sp%>%
  separate(`host lineage`, sep="; ",
           into = c("domain", "phylum", "class", "order", "family.etc"), 
           extra = "merge")



# Remove duplicate hosts --------------------------------------------------


d.sp %>% 
  group_by(tax.id, `virus name`)%>%
  summarise(n=n())%>%
  arrange(desc(n)) %>% 
  ggplot(aes(x=n))+
  geom_histogram()+
  scale_y_log10()


# Back to phage duplicates due to multiple hosts
dup <- d.sp%>%
  group_by(tax.id, `virus name`)%>%
  summarise(n.host=n())%>%
  filter(n.host>1)

#there are 246 phages with duplicate hosts
# how different are the hosts?
# Are they from different orders?

dup.order <- d.sp%>%
  group_by(tax.id, `virus name`, domain, phylum, class, order)%>%
  summarise(n.host=n())%>%
  filter(n.host>1) %>% 
  group_by(tax.id, `virus name`, domain, phylum, class) %>% 
  summarise(n.host.order=n())%>%
  filter(n.host.order>1)

#there is only one such phage: PRD1. This has also hosts from different orders within the gammaproteobacteria.
# otherwise, all hosts for each phage are the same down to order

# Look at host species
dup.hosts <- 
  d.sp %>% 
  filter(tax.id %in% dup$tax.id) %>% 
  # filter(phylum=="Firmicutes") %>% 
  group_by(tax.id) %>% 
  mutate(host.num = paste0("host_",row_number())) %>% 
  ungroup() %>% 
  select(tax.id, `virus name`, host.num, `host name`) %>% 
  pivot_wider(names_from = host.num, values_from = `host name`)
# in Firmicutes phages multiple hosts are mostly from the same species or 
# species complex (antracis, cereus). 
# SPbeta infects subtilis and pumilus, in same genus.

#remove host duplicate
#add host number to assist and to remember
d.sp <- dup%>%
  select(tax.id,n.host)%>%
  left_join(d.sp, ., by="tax.id")
#add 1 host for all the rest
d.sp$n.host[is.na(d.sp$n.host)] <- 1

d.sp <- d.sp%>%
   filter(!duplicated(tax.id))


# add host to sigma factor data

d.faa <- left_join(d.faa,d.sp,by=c("taxon"="tax.id"))


# # save data
# write.csv(here("vogdb","data","vog_sigma_clean_Whost.csv"))
# save(d.faa,file = here("vogdb","data","vog_sigma_clean_Whost.RData"))
# # load(file = here("vogdb","data","vog_sigma_clean_Whost.RData"))





# sigma factor gene number per phage --------------------------------------

d.sp <-
  d.faa%>%
  group_by(taxon)%>%
  summarise(n.sigma=n()) %>%
  left_join(d.sp, ., by = c("tax.id" = "taxon"))

# add 0 to phages without sigma genes
d.sp$n.sigma[is.na(d.sp$n.sigma)] <- 0

#How many unique viruses have sigma factors?
length(unique(d.faa$taxon))
#471

# Phylum independent summary of sigma factors/genome
d.sp%>%
  group_by(tax.id)%>%
  group_by(n.sigma)%>%
  summarise(n.genomes=n())%>%
  mutate(perc=100*n.genomes/sum(n.genomes))

# n.sigma n.genomes   perc
# 0      2962   86.3
# 1       397   11.6
# 2        50    1.46
# 3        24     0.699


# summarise by phylum
d.sp%>%
  group_by(phylum)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>%
  summarise(n=n(),
            perc_sigma=100*sum(has_sigma)/n(),
            perc_multi=100*sum(multi_sigma)/n())

#   phylum                                  n perc_sigma perc_multi
#   <chr>                               <int>      <dbl>      <dbl>
# 1 Actinobacteria                        797       1.88      0
# 2 Aquificae                               1       0         0
# 3 Bacteroidetes/Chlorobi group           76       2.63      0
# 4 Cyanobacteria/Melainabacteria group   102      64.7       3.92
# 5 Deinococcus-Thermus                     8       0         0
# 6 Firmicutes                            745      18.9       7.79
# 7 Proteobacteria                       1683      14.7       0.713
# 8 PVC group                               6       0         0
# 9 Spirochaetes                            3       0         0
# 10 Tenericutes                             9       0         0
# 11 NA                                      3       0         0


d.sp%>%
  filter(n.sigma>0) %>%
  group_by(phylum)%>%
  summarise(n= n(), avg_Nsig=mean(n.sigma), sd=sd(n.sigma),
            perc_multi = sum(n.sigma>1)/n(),
            perc_single = sum(n.sigma==1)/n())


# _____________--------------------
# Sigma presence plots ----------------------------------------------------

# summarise by viral family
d.presence <- d.sp%>%
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
  filter(! is.na(viral.family))  %>%
  filter(! is.na(phylum))


# > plot by host phyla ----------------------------------------------------

p.phylum <-  d.presence %>% 
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
  geom_col(aes(y=w.sigma), position=position_dodge(preserve = "single"),size=0, width = 0.6, fill="grey10")+
  geom_col(aes(y=n), position=position_dodge(preserve = "single"),
           fill = alpha("grey20", 0.2), color="black", size = 0.7, width=0.6)+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,2000))+
  xlab("Host Phylum")+
  theme_classic(base_size = 12)+
  panel_border(color = "black", size = 1)+
  coord_flip(expand = F)+
  scale_fill_viridis_d()+
  scale_x_discrete()

# > plot by viral family -------------------------------------------------

p.Vfam <-  d.presence %>% 
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
  geom_col(aes(y=w.sigma), position=position_dodge(preserve = "single"),size=0, width = 0.6, fill="grey10")+
  geom_col(aes(y=n), position=position_dodge(preserve = "single"),
           fill = alpha("grey20", 0.2), color="black", size = 0.7, width=0.6)+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,2000))+
  xlab("Viral Family")+
  theme_classic(base_size = 12)+
  panel_border(color = "black", size = 1)+
  coord_flip(expand = F)+
  scale_fill_viridis_d()

# > plot by host phyla AND viral family --------------------------------


p.both.data <-  d.presence %>% 
  mutate(phylum = str_replace(phylum, "Deinococcus-Thermus","Deino.-Thermus")) %>% 
  group_by(viral.family, phylum,)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>% 
  summarise(n=n(),
            w.sigma=sum(has_sigma),
            perc_sigma=100*sum(has_sigma)/n(),
            w.multi=sum(multi_sigma),
            perc_multi=100*sum(multi_sigma)/n())
  
# phyla that have any sigma factor
phyla.keep <- p.both.data %>%
  group_by(phylum) %>%
  summarise(sig=sum(w.sigma)) %>% 
  filter(sig>2) %>% 
  pull(phylum)

# viral families that have any sigma factor
viral.keep <- p.both.data %>%
  group_by(viral.family) %>%
  summarise(sig=sum(w.sigma)) %>% 
  filter(sig>1) %>% 
  pull(viral.family) %>% as.character()

#plot
p.both <- 
  p.both.data %>% 
  filter(phylum %in% phyla.keep) %>% 
  filter(viral.family %in% viral.keep) %>% 
  #removing a single phage that is unclassified
  filter(viral.family != "unclassified Caudovirales") %>% 
  mutate(viral.family = fct_infreq(viral.family) ) %>% 
  ggplot(aes(viral.family))+
  geom_col(aes(y=w.sigma), position=position_dodge(preserve = "single"),size=0, width = 0.6, fill="grey10")+
  geom_col(aes(y=n), position=position_dodge(preserve = "single"),
           fill = alpha("grey20", 0.2), color="black", size = 0.7, width=0.6)+
  ylab("No. Phage Genomes") +
  scale_y_log10(limits = c(1,2000))+
  xlab("Viral Family")+
  theme_classic(base_size = 12)+
  panel_border(color = "black", size = 0.7)+
  # facet_wrap(~phylum, nrow = 1)+
  facet_nested_wrap("Host Phylum" + phylum~., nest_line = TRUE, nrow = 1)+
  coord_flip(expand = F)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+

  theme(#axis.text.y = element_blank(),
        legend.position = c(0.9,0.22),
        legend.text = element_text(size = 8),
        legend.key.height = unit(1,"mm"),
        strip.background = element_blank())+
        # strip.background = element_rect(size = 0.7))+
  guides(fill =  guide_legend(title = "Viral Family", 
                              reverse = TRUE),
         color =  guide_legend(title = "Viral Family", 
                              reverse = TRUE))


# n sigma panel -----------------------------------------------------------

#summarise data
n.sig.data <- 
  d.sp %>% 
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  group_by(phylum, n.sigma) %>% 
  summarise(n=n(), .groups = "drop") %>% 
  group_by(phylum) %>% 
  summarise( n.sigma= n.sigma, n = n, total = sum(n), .groups = "drop") %>% 
  mutate(pnl = "Host Phylum")

n.sig.all <-   d.sp %>% 
  group_by(n.sigma) %>% 
  summarise(n=n(), .groups = "drop", ) %>% 
  summarise( n.sigma= n.sigma, n = n, total = sum(n), .groups = "drop") %>% 
  mutate(pnl = "All Hosts", phylum = "Bacteria") %>% 
  bind_rows(.,filter(n.sig.data, phylum %in% phyla.keep)) %>% 
  mutate(perc = n/total)



n.labs <- n.sig.all %>% 
  group_by(pnl,phylum) %>% 
  summarise(n=sum(n)) %>% 
  mutate(lab = paste("n=", n))


p.multi <-  
  n.sig.all %>% 
  ggplot( aes(n.sigma, group = phylum)) + 
  geom_col(aes(y = perc),  fill = "grey10", color="black", size = 0.7, width=0.6) + 
  geom_text(data = n.labs, aes(label=lab), x=2,y=0.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap(pnl + phylum~., nest_line = TRUE, nrow = 1)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")+
  theme(strip.background = element_blank())


# > combine plots ----------------------------------------------------------


top_row <- plot_grid(p.phylum, NULL, p.Vfam, nrow = 1,
                     labels = c("a","","b"),rel_widths = c(1,0.2,1))
# all.3 <- plot_grid(top_row, NULL, p.both, ncol = 1, labels = c("","c"), 
#                    rel_heights = c(1,0.1,1)) 
middle_row <- plot_grid(NULL, p.both, nrow = 1,rel_widths = c(0.1,1))

all.4 <- plot_grid(top_row, middle_row,p.multi, ncol = 1, labels = c("","c", "d")) 

ggsave2(here("vogdb","figures","viral_family_host_phylum.png"),
        plot = ggdraw(all.4) +
          theme(plot.background = element_rect(fill="white", color = NA)),
        width = 8,height = 8)

# > stats for sigma presence --------------------------------------------

# >> Host phyla -----------------------------------------------------------


# make contingency table 
contin.t1 <- 
  d.sp%>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  # remove singeltons
  group_by(phylum) %>% filter( n() > 1 ) %>% ungroup() %>% 
  mutate(has.sigma = n.sigma>0) %>% 
  select(phylum, has.sigma)%>%
  table()

m1 <- chisq.test(contin.t1, simulate.p.value = TRUE, B = 1e6)
m2 <- fisher.test(contin.t1, simulate.p.value = TRUE, B = 1e6)

write.csv(contin.t1, here("vogdb/data/","Hphyla_sigma_contingency.csv"))

# >> viral family -------------------------------------------------------------



# make contingency table 
contin.t2 <- 
  d.sp%>%
  #extract viral family
  mutate(viral.family=str_extract(`virus lineage`,
                                  regex("caudovirales;.*", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex("caudovirales; ", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex(";.*", ignore_case = T))) %>%
  mutate(viral.family = if_else(is.na(viral.family), "non-Caudovirales", viral.family)) %>% 
  # remove singeltons
  group_by(viral.family) %>% filter( n() > 1 ) %>% ungroup() %>%
  
  mutate(has.sigma = n.sigma>0) %>% 
  select(viral.family, has.sigma)%>%
  table()

m1 <- chisq.test(contin.t2, simulate.p.value = TRUE, B = 1e6)
m2 <- fisher.test(contin.t2, simulate.p.value = TRUE, B = 1e6)

write.csv(contin.t2, here("vogdb/data/","VIRfam_sigma_contingency.csv"))

# _____________--------------------


# Sigma number plots ------------------------------------------------------


# > Plot by phylum ----------------------------------------------------------


# filter for phyla with less than 20 phages
phyla_rm <-  d.sp %>% 
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  group_by(phylum) %>% 
  summarise(n=n()) %>% 
  filter(n<20) %>% 
  pull(phylum)


# summary table for N by phylum
n.labs <- d.sp%>%
  filter(! phylum %in% phyla_rm) %>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  group_by(phylum)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))



p.phylum <- 
  d.sp%>%
  filter(! phylum %in% phyla_rm) %>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  ggplot( aes(n.sigma, group = phylum)) + 
  geom_bar(aes(y = ..prop..), stat="count", fill = viridis(2)[1], color="black") + 
  geom_text(data=n.labs, aes(label=lab), x=3,y=0.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap("Host Phylum" + phylum~., nest_line = TRUE, ncol = 1)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")

# > Plot by Firmicute genera ------------------------------------------------


# extracting genus data for  firmicutes
firmi <- d.sp%>%
  filter(str_detect(phylum,regex("firmicutes", ignore_case = T)))%>%
  # filter(!str_detect(order,regex("lacto", ignore_case = T)))%>%
  separate(family.etc, into = c("family","genus","species","strain"),sep=";")

firmi$.species.name[is.na(firmi$genus)]
#   "Streptococcus phage MM1" "Clostridium phage phiCDHM11"  ... 
# genus is also in viral sp name as first word
firmi <- firmi%>%
  mutate(genus2=str_extract(.species.name,regex(".*? ")) %>% trimws()) %>% 
  mutate(genus = trimws(genus)) %>% 
  mutate(genus.plot = if_else(is.na(genus), genus2, genus))

firmi%>%
  group_by(genus)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>% 
  summarise(n=n(),
            perc_sigma=100*sum(has_sigma)/n(),
            perc_multi=100*sum(multi_sigma)/n())

firmi%>%
  filter(n.sigma>0) %>%
  group_by(genus)%>%
  summarise(n= n(), avg_Nsig=mean(n.sigma), sd=sd(n.sigma),
            perc_multi = 100*sum(n.sigma>1)/n(),
            perc_single = 100*sum(n.sigma==1)/n())

# firmi %>% 
#   select(genus, genus2) %>% 
#   filter(genus != genus2)

# summary table for N by order
n.labs <- firmi%>%
  mutate(genus.plot = str_replace_all(genus.plot, "bacillus", "bacil.")) %>% 
  mutate(genus.plot = str_replace_all(genus.plot, "bacterium", "bact.")) %>% 
  group_by(tax.id,genus.plot)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(genus.plot)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) %>% 
  filter(n>5)

p.genus <-
firmi%>%
  mutate(genus.plot = str_replace_all(genus.plot, "bacillus", "bacil.")) %>% 
  mutate(genus.plot = str_replace_all(genus.plot, "bacterium", "bact.")) %>% 
  filter(genus.plot %in% n.labs$genus.plot) %>% 
  group_by(genus.plot, n.sigma) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>% 
  
  ggplot( aes(n.sigma, perc)) + 
  geom_col(fill = viridis(2)[2], color = "black") + 
  geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap(~"Host genus (Firmicutes)" + genus.plot, 
                    scales = "fixed", nrow = 5)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")


  # ggsave2(here("vogdb","figures","sigma_HostGenus_firmicutes.png"),
  #         plot = p.genus+
  #           ggtitle("phage sigma factor content by host genus",
  #                   "Firmicute infecting phages"),
  #         width = 7,height = 7)




# > Plot by viral type ------------------------------------------------------



# First I will break up the virus lineage data so
# the host lineage is not uniform having 2-9 levels
str_count(firmi$`virus lineage`,";")%>%range()
str_count(firmi$`virus lineage`,";")%>%hist()
# the last looks like this: 
# Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes;
#  Caudovirales; Siphoviridae; unclassified Siphoviridae

# All but 1 belong to the tailed phages (Caudovirales)
# the single exception is Salisaeta icosahedral phage 1
# indeed this is " a tailless bacteriophage" (PMID: 22509017)
  tmp <- (!str_detect(d.faa$`virus lineage`,"Caudovirales"))%>%which()
  d.faa$`virus name`[tmp]
# I will only look at the family level below tailed phage (order Caudovirales)
firmi <- firmi%>%
  mutate(viral.family=str_extract(`virus lineage`,
                                  regex("caudovirales;.*", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                  regex("caudovirales; ", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex(";.*", ignore_case = T)))

# d.sp$`virus name`[is.na(d.sp$viral.family)]
# 110/3433 left without family assignments


# # make contingency table 
# contin.t <-
#   firmi%>%
#   #summarise sigma factors per phage
#   group_by(tax.id,viral.family)%>%
#   summarise(n.sig=n())%>%
#   ungroup()%>%
#   select(-tax.id)%>%
#   table()
# 
# m1 <- chisq.test(contin.t)
# m2 <- fisher.test(contin.t,workspace = 6e8)


# summary table for N by viral.family
n.labs <- firmi %>% 
  filter(str_detect(genus.plot,"Bacillus"))%>%
  filter(!is.na(viral.family)) %>% 
  group_by(viral.family)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) %>% 
  filter(n>1)



d.vir <-
  firmi%>%
  filter(str_detect(genus.plot,"Bacillus"))%>%
  filter(viral.family %in% n.labs$viral.family) %>% 
  group_by(viral.family, n.sigma ,phylum)%>%
    summarise(n.genomes=n())%>%
  ungroup()%>%
  group_by(viral.family)%>%
    mutate(perc=n.genomes/sum(n.genomes))%>%

  filter(! phylum %in% phyla_rm) %>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  
  mutate(viral.family=fct_reorder(viral.family,n.genomes))

p.vir <- d.vir %>% 
  ggplot( aes(n.sigma,perc)) + 
  geom_col(fill = "green4") + 
  geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
  scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap(~"Viral Family (Bacillus host)" + viral.family, nrow = 5)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")
  

  # ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
  #         plot = p.vir+
  #           ggtitle("phage sigma factor content by viral family"),
  #         width = 10,height = 6)

# Who are the siphoviruses with 3 sigma factors?
firmi%>%
  filter(viral.family=="Siphoviridae")%>%
  filter(n.sigma>2)%>%
  arrange(desc(n.sigma)) %>% 
  select(`virus name` , n.sigma)



# `virus name`                   n.sigma
# <chr>                            <int>
# 1 Bacillus phage Slash                 3
# 2 Bacillus phage Stahl                 3
# 3 Bacillus phage Staley                3
# 4 Bacillus phage Stills                3
# 5 Bacillus phage BceA1                 2
# 6 Brevibacillus phage Jenst            2
# 7 Cyanophage PSS2                      2
# 8 Staphylococcus phage SpaA1           2
# 9 Synechococcus phage S-CBS2           2

# They are all Bacillus phages(n=4)
# with 2 sigmas there are other hosts


# >> combine plots -----------------------------------------------------------


# right_col <- 
#   plot_grid(p.vir, NULL,
#             rel_heights = c(4.5, 1),
#             ncol = 1)

p <-  plot_grid(p.phylum,p.genus,#right_col,
          rel_widths = c(1.5, 3), #1.3),
          nrow = 1, labels = LETTERS)
  ggsave2(here("vogdb","figures","sigma_taxonomy.png"),
          plot = p, width = 12,height = 10)


# > Herelle host genera -----------------------------------------------------
  d.herelle <- filter(firmi, phage.nonphage == "phage") %>% 
    filter(str_detect(`virus lineage`, "Herelle"))
  
  
  # summary table for N by viral.family
  n.labs <- d.herelle %>% 
    filter(!is.na(viral.family)) %>% 
    group_by(genus.plot, viral.family)%>%
    summarise(n=n())%>%
    mutate(lab=paste0("n=",n)) %>% 
    filter(n>1) %>%
    ungroup()
  
  
  
  d.vir <-
    d.herelle%>%
    filter(genus.plot %in% n.labs$genus.plot) %>% 
    filter(viral.family %in% n.labs$viral.family) %>% 
    group_by(viral.family, n.sigma ,genus.plot)%>%
    summarise(n.genomes=n())%>%
    ungroup()%>%
    group_by(genus.plot)%>%
    mutate(perc=n.genomes/sum(n.genomes))
  
  p.herelle <- d.vir %>% 
    ggplot( aes(n.sigma,perc)) + 
    geom_col(fill = viridis(2)[1], color = "black") + 
    geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    facet_wrap(~ genus.plot )+
    # facet_nested_wrap(~"Viral Family (Bacillus host)" + viral.family, nrow = 5)+
    theme_classic(base_size = 13)+
    panel_border(color = "black")
  
  
  # ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
  #         plot = p.vir+
  #           ggtitle("phage sigma factor content by viral family"),
  #         width = 10,height = 6)
  
  # > sipho host genera -----------------------------------------------------
  d.sipho <- filter(firmi, phage.nonphage == "phage") %>% 
    filter(str_detect(`virus lineage`, "Sipho")) %>% 
    # merge Clostridium and Clostridioides host genera
    mutate(genus.plot = str_replace(genus.plot,"Clostridioides","Clostridium"))
    # Oren, A., & Rupnik, M. (2018). 
    # Clostridium difficile and Clostridioides difficile: two validly published and correct names. 
    # Anaerobe, 52, 125-126.
  
  
  # summary table for N by viral.family
  n.labs <- d.sipho %>% 
    filter(!is.na(viral.family)) %>% 
    group_by(genus.plot, viral.family)%>%
    summarise(n=n())%>%
    mutate(lab=paste0("n=",n)) %>% 
    filter(n>1) %>%
    ungroup()
  
  
  
  d.vir <-
    d.sipho%>%
    filter(genus.plot %in% n.labs$genus.plot) %>% 
    filter(viral.family %in% n.labs$viral.family) %>% 
    group_by(viral.family, n.sigma ,genus.plot)%>%
    summarise(n.genomes=n())%>%
    ungroup()%>%
    group_by(genus.plot)%>%
    mutate(perc=n.genomes/sum(n.genomes))
  
  p.sipho <- d.vir %>% 
    ggplot( aes(n.sigma,perc)) + 
    geom_col(fill = viridis(2)[2], color = "black") + 
    geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    facet_wrap(~ genus.plot )+
    # facet_nested_wrap(~"Viral Family (Bacillus host)" + viral.family, nrow = 5)+
    theme_classic(base_size = 13)+
    panel_border(color = "black")
  
  
  # ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
  #         plot = p.vir+
  #           ggtitle("phage sigma factor content by viral family"),
  #         width = 10,height = 6)
  
  # >> combine plots -----------------------------------------------------------
  
  right_col <-
    plot_grid(p.herelle, NULL,
              rel_widths = c(1, 0.5),
              ncol = 2)
  
  p <-  plot_grid(right_col,p.sipho,
                  rel_heights =  c(0.6, 1),
                  nrow = 2, labels = LETTERS)
  
  p <- ggdraw(p) + 
    theme(plot.background = element_rect(fill="white", color = NA))
  
  ggsave2(here("vogdb","figures","viral_family_Nsig.png"),
          plot = p, width = 7,height = 8)
  
  
# > Viruses of all phyla ----------------------------------------------------
  d.phage <- filter(d.sp, phage.nonphage == "phage")
  # First I will break up the virus lineage data so
  # the host lineage is not uniform having 2-9 levels
  str_count(d.phage$`virus lineage`,";")%>%range()
  str_count(d.phage$`virus lineage`,";")%>%hist()
  # the last looks like this: 
  # "Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; 
  # Caudoviricetes; Caudovirales; Myoviridae; Phikzvirus; unclassified Phikzvirus"
  str_detect(d.phage$`virus lineage`, 'Caudovirales')[tmp] %>% sum()
  
  tmp <- d.phage %>% 
    filter(!str_detect(`virus lineage`, 'Caudovirales'))
  # All but 97 belong to the tailed phages (Caudovirales)
  
  d.phage <- d.phage %>% 
    filter(str_detect(`virus lineage`, 'Caudovirales'))
  # I will only look at the family level below tailed phage (order Caudovirales)
  d.phage <- d.phage%>%
    mutate(viral.family=str_extract(`virus lineage`,
                                    regex("caudovirales;.*", ignore_case = T)))%>%
    mutate(viral.family=str_remove(viral.family,
                                   regex("caudovirales; ", ignore_case = T)))%>%
    mutate(viral.family=str_remove(viral.family,
                                   regex(";.*", ignore_case = T)))
  
  # summary table for N by viral.family
  n.labs <- d.phage %>% 
    filter(!is.na(viral.family)) %>% 
    filter(n.sigma>0) %>% 
    group_by(phylum, viral.family)%>%
    summarise(n=n())%>%
    mutate(lab=paste0("n=",n)) %>% 
    filter(n>1) %>% 
    # name adjust
    mutate(phylum =
             case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                       str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                       TRUE ~ phylum)) %>% 
    ungroup()
  
  
  
  d.vir <-
    d.phage%>%
    # name adjust
    mutate(phylum =
             case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                       str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                       TRUE ~ phylum)) %>% 
    filter(n.sigma>0) %>% 
    filter(phylum %in% n.labs$phylum) %>% 
    filter(viral.family %in% n.labs$viral.family) %>% 
    group_by(viral.family, n.sigma ,phylum)%>%
    summarise(n.genomes=n())%>%
    ungroup()%>%
    group_by(viral.family)%>%
    mutate(perc=n.genomes/sum(n.genomes))%>%
    
    filter(! phylum %in% phyla_rm) %>%
    
    mutate(viral.family=fct_reorder(viral.family,n.genomes))
  
  p.vir <- d.vir %>% 
    ggplot( aes(n.sigma,perc)) + 
    geom_col(fill = "green4") + 
    geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    facet_grid(viral.family ~ phylum )+
    # facet_nested_wrap(~"Viral Family (Bacillus host)" + viral.family, nrow = 5)+
    theme_classic(base_size = 13)+
    panel_border(color = "black")
  
  
  # ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
  #         plot = p.vir+
  #           ggtitle("phage sigma factor content by viral family"),
  #         width = 10,height = 6)
  
# > Summary plot - MAIN ------------------------------------------------------

# Phylum independent summary of sigma factors/genome
  nsig_all <- d.sp%>%
    group_by(tax.id)%>%
    group_by(n.sigma)%>%
    summarise(n.genomes=n())%>%
    mutate(perc=n.genomes/sum(n.genomes))
  
  nsig_firmi <- d.sp%>%
    filter(phylum=="Firmicutes") %>% 
    group_by(tax.id)%>%
    group_by(n.sigma)%>%
    summarise(n.genomes=n())%>%
    mutate(perc=n.genomes/sum(n.genomes))  
  
  nsig_bacil <- firmi%>%
    filter(genus.plot=="Bacillus") %>% 
    group_by(tax.id)%>%
    group_by(n.sigma)%>%
    summarise(n.genomes=n())%>%
    mutate(perc=n.genomes/sum(n.genomes)) 
  
  
 dplot <- 
    bind_rows(
      nsig_all %>% mutate(host_tax = "Bacteria"),
      nsig_firmi %>% mutate(host_tax = "Firmicutes"),
      nsig_bacil %>% mutate(host_tax = "Bacillus")
     ) %>% 
    mutate(host_tax = fct_relevel(host_tax, c("Bacteria", "Firmicutes", "Bacillus")))

  p.tax <- dplot %>% 
    ggplot(aes(n.sigma, perc, fill = host_tax))+
    geom_col(color = "black", position = position_dodge2(padding = 0.4), width = 0.7)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage Genomes") +
    xlab("Sigma Factors per Genome")+
    theme_classic()+
    panel_border(color = "black")+
    # scale_fill_grey(start = 0.8,end = 0.2, name = "Host taxa")+
    scale_fill_viridis_d( name = "Host Taxa", direction = -1)+
    theme(legend.position = c(1,1),
          legend.justification = c(1,1),
          legend.background = element_blank(),
          legend.key.size = unit(0.1, 'in'))
          
  # ggsave(filename = here("vogdb","figures","nSigma_sum.png"),
  #       plot = p.tax, width = 3,height = 3)
  # summary table for N by viral.family
  n.labs <- firmi %>% 
    filter(str_detect(genus.plot,"Bacillus"))%>%
    filter(!is.na(viral.family)) %>% 
    group_by(viral.family)%>%
    summarise(n=n())%>%
    mutate(lab=paste0("n=",n)) %>% 
    filter(n>1)
  
  
  
  d.vir <-
    firmi%>%
    filter(str_detect(genus.plot,"Bacillus"))%>%
    filter(viral.family %in% n.labs$viral.family) %>% 
    group_by(viral.family, n.sigma ,phylum)%>%
    summarise(n.genomes=n())%>%
    ungroup()%>%
    group_by(viral.family)%>%
    mutate(perc=n.genomes/sum(n.genomes))%>%
    
    filter(! phylum %in% phyla_rm) %>%
    # name adjust
    mutate(phylum =
             case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                       str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                       TRUE ~ phylum)) %>% 
    
    mutate(viral.family=fct_reorder(viral.family,n.genomes))
  
  p.vir2 <- d.vir %>% 
    mutate(viral.family = fct_relevel(viral.family, c("Podoviridae", "Myoviridae",
                                                  "Siphoviridae", "Herelleviridae"))) %>% 
    ungroup() %>% 
    select(viral.family, n.sigma, perc) %>% 
    complete(n.sigma, viral.family, fill = list(perc=0)) %>% 
    ggplot(aes(n.sigma, perc, fill = viral.family))+
    geom_col(color = "black", position = position_dodge2(padding = 0.4), width = 0.7)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage Genomes") +
    xlab("Sigma Factors per Genome")+
    theme_classic()+
    panel_border(color = "black")+
    # scale_fill_grey(start = 0.8,end = 0.2, name = "Host taxa")+
    scale_fill_viridis_d( name = "Viral Family", direction = -1)+
    theme(legend.position = c(1,1),
          legend.justification = c(1,1),
          legend.background = element_blank(),
          legend.key.size = unit(0.1, 'in'))
  
  # ggsave(filename = here("vogdb","figures","nSigma_virFam.png"),
  #        plot = p.vir2, width = 3,height = 3)
  
  plot_grid(p.tax,p.vir2, labels = letters) %>% 
save_plot(filename = here("vogdb","figures","nSigmaAB.png"),
          base_width = 6 , base_height = 3)
  