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
vog.sp <-  read_tsv(here("vogdb/data","vogdb_downloads","vog.species.list"), trim_ws = T) %>% 
  as_tibble( .name_repair = "universal") %>% 
  # keep only phages
  filter(phage.nonphage == "phage")

#import virus-host data
download.file(url = "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv",
              destfile = here("vogdb/data", "virushostdb.tsv"))
# downloaded on (6/JAN/2022)
vh.db <- read_tsv(here("vogdb/data","virushostdb.tsv"))


# virus duplicates in VHDB data -------------------------------------------
# These reflect multiple hosts
duplicated(vh.db$`virus tax id`)%>%sum() # 4042
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
  # focus on viruses of bacteria (remove Archaea)
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

#there are 275 phages with duplicate hosts
# how different are the hosts?
# Are they from different orders?

dup.order <- d.sp%>%
  group_by(tax.id, `virus name`, domain, phylum, class, order)%>%
  summarise(n.host=n())%>%
  filter(n.host>1) %>% 
  group_by(tax.id, `virus name`, domain, phylum, class) %>% 
  summarise(n.host.order=n())%>%
  filter(n.host.order>1)

#there is only one such phage: PRD1. 
# This has also hosts from different orders within the gammaproteobacteria.
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
# In many cases hosts are differen species within genus.
# E.g. SPbeta infects subtilis and pumilus, in same genus of Bacillus.
# An exception I noted is 'Thermus phage phi OH2' that has listed 
# hosts from two different phyla:  Deinococcus-Thermus and Firmicutes.



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
write_csv(select(d.faa,-seq), here("vogdb","data","vog_sigma_clean_Whost.csv"))
save(d.faa,file = here("vogdb","data","vog_sigma_clean_Whost.RData"))
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
#617

# Phylum independent summary of sigma factors/genome
d.sp%>%
  group_by(tax.id)%>%
  group_by(n.sigma)%>%
  summarise(n.genomes=n())%>%
  mutate(perc=100*n.genomes/sum(n.genomes))

# n.sigma n.genomes   perc
# 0       3770        86.4  
# 1        518        11.9  
# 2         47         1.08 
# 3         30         0.687


# summarise by phylum
d.sp%>%
  group_by(phylum)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>%
  summarise(n=n(),
            perc_sigma=100*sum(has_sigma)/n(),
            perc_multi=100*sum(multi_sigma)/n())

# phylum                                  n perc_sigma perc_multi
# <chr>                               <int>      <dbl>      <dbl>
# 1 Actinobacteria                       1314      0.609      0    
# 2 Aquificae                               1      0          0    
# 3 Bacteroidetes/Chlorobi group           83      2.41       0    
# 4 Cyanobacteria/Melainabacteria group   101     65.3        1.98 
# 5 Deinococcus-Thermus                     5      0          0    
# 6 Firmicutes                            865     21.6        7.28 
# 7 Fusobacteria                            1      0          0    
# 8 Proteobacteria                       1980     16.8        0.606
# 9 PVC group                               6      0          0    
# 10 Spirochaetes                            3      0          0    
# 11 Tenericutes                             3      0          0    
# 12 NA                                      3      0          0 


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
  scale_y_log10(limits = c(1,3000))+
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
  scale_y_log10(limits = c(1,3000))+
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
  scale_y_log10(limits = c(1,3000))+
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


# > n sigma panel -----------------------------------------------------------

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
  mutate(lab = paste0("n=", n))


p.multi <-  
  n.sig.all %>% 
  ggplot( aes(n.sigma, group = phylum)) + 
  geom_col(aes(y = perc),  fill = "grey80", color="black", size = 0.7, width=0.6) + 
  geom_text(data = n.labs, aes(label=lab), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+
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


all.4 <- plot_grid(top_row, p.both,p.multi, ncol = 1, labels = c("","c", "d")) 

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
# Plot by Firmicute genera ------------------------------------------------


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
  mutate(genus.plot = str_replace_all(genus.plot, "bacillus", "bacill.")) %>%
  mutate(genus.plot = str_replace_all(genus.plot, "bacterium", "bact.")) %>%
  mutate(genus.plot = fct_infreq(genus.plot)) %>% 
  group_by(tax.id,genus.plot)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(genus.plot)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) 

p.genus <-
  firmi%>%
  mutate(genus.plot = str_replace_all(genus.plot, "bacillus", "bacill.")) %>%
  mutate(genus.plot = str_replace_all(genus.plot, "bacterium", "bact.")) %>%
  filter(genus.plot %in% n.labs$genus.plot) %>% 
  mutate(genus.plot = fct_infreq(genus.plot)) %>% 
  group_by(genus.plot, n.sigma) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>% 
  
  ggplot( aes(n.sigma, perc)) + 
  geom_col(fill = "grey80", color = "black") +
  geom_text(data = n.labs  ,aes(label=lab), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~genus.plot, scales = "fixed", ncol = 4)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")+
  theme(strip.background = element_blank())

ggsave2(here("vogdb","figures","sigma_taxonomy.png"),
        plot = plot_grid(p.genus, labels = "a"), width = 7,height = 8)

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

# firmi$`virus name`[is.na(firmi$viral.family)]
# 5/865 left without family assignments

# all Firmicutes phages ---------------------------------------------------

n.labs <- firmi %>% 
  # filter(str_detect(genus.plot,"Bacillus"))%>%
  filter(!is.na(viral.family)) %>% 
  group_by(viral.family)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) %>% 
  filter(n>1)

d.vir <-
  firmi%>%
  # filter(str_detect(genus.plot,"Bacillus"))%>%
  filter(viral.family %in% n.labs$viral.family) %>% 
  group_by(viral.family, n.sigma ,phylum)%>%
  summarise(n.genomes=n())%>%
  ungroup()%>%
  group_by(viral.family)%>%
  mutate(perc=n.genomes/sum(n.genomes))%>%
  
  # name adjust
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
  theme_classic(base_size = 13)+
  panel_border(color = "black")+
  # scale_fill_grey(start = 0.8,end = 0.2, name = "Host taxa")+
  scale_fill_viridis_d( name = "Viral Family", direction = 1)+
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key.size = unit(0.15, 'in'))

# > Herelle host genera -----------------------------------------------------
d.herelle <- filter(firmi, phage.nonphage == "phage") %>% 
  filter(str_detect(`virus lineage`, "Herelle"))


# summary table for N by viral.family
n.labs <- d.herelle %>% 
  filter(!is.na(viral.family)) %>% 
  mutate(genus.plot = fct_infreq(genus.plot)) %>%
  group_by(genus.plot, viral.family)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) %>% 
  filter(n>1) %>%
  ungroup()



d.vir <-
  d.herelle%>%
  filter(genus.plot %in% n.labs$genus.plot) %>% 
  mutate(genus.plot = fct_infreq(genus.plot)) %>%
  filter(viral.family %in% n.labs$viral.family) %>% 
  group_by(viral.family, n.sigma ,genus.plot)%>%
  summarise(n.genomes=n())%>%
  ungroup()%>%
  group_by(genus.plot)%>%
  mutate(perc=n.genomes/sum(n.genomes))

p.herelle <- d.vir %>% 
  ggplot( aes(n.sigma,perc)) + 
  geom_col(fill = viridis(7)[4], color = "black") + 
  geom_text(data = n.labs, aes(label=lab), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+  
  scale_y_continuous(labels=scales::percent, limits = c(0,1),  breaks = c(0,.5,1)) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~ genus.plot,ncol = 2 )+
  # facet_nested_wrap(~"Host genus (Firmicutes)" + genus.plot, nrow = 2)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))


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
  mutate(genus.plot = fct_infreq(genus.plot)) %>%
  group_by(genus.plot, viral.family)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n)) %>% 
  filter(n>1) %>%
  ungroup()



d.vir <-
  d.sipho%>%
  filter(genus.plot %in% n.labs$genus.plot) %>% 
  mutate(genus.plot = fct_infreq(genus.plot)) %>%
  filter(viral.family %in% n.labs$viral.family) %>% 
  group_by(viral.family, n.sigma ,genus.plot)%>%
  summarise(n.genomes=n())%>%
  ungroup()%>%
  group_by(genus.plot)%>%
  mutate(perc=n.genomes/sum(n.genomes))



p.sipho <- d.vir %>% 
  ggplot( aes(n.sigma,perc)) + 
  geom_col(fill = viridis(7)[3], color = "black") + 
  geom_text(data = n.labs, aes(label=lab), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5)+  
  scale_y_continuous(labels=scales::percent, limits = c(0,1), breaks = c(0,.5,1)) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~ genus.plot )+
  # facet_nested_wrap(~"Host genus (Firmicutes)" + genus.plot, nrow = 4)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 10))


# > combine plots -----------------------------------------------------------

row1 <-  plot_grid(p.vir2,NULL, p.herelle,labels = c("b","","c"),
                   rel_widths =  c(1.1,0.1, 1),nrow = 1)
# row2 <-  plot_grid(p.herelle,NULL,
#                    rel_widths =  c(1, 0.5),nrow = 1)

p <- plot_grid(row1, NULL, p.sipho, labels = c("","","d"),
               ncol = 1, rel_heights =  c(1,0.05,1.1))

ggsave2(here("vogdb","figures","viral_family_Nsig.png"),
        plot = p, width = 7,height = 8)

