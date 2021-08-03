library(tidyverse)
library(here)
library(cowplot)
library(ggh4x)

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

##############################################
#  virus duplicates in VHDB data
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



#### remove duplicate hosts ####

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



#### sigma factor gene number per phage ####
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


##########################
# Plot by phylum
##########################


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
  geom_bar(aes(y = ..prop..), stat="count", fill = "blue4") + 
  geom_text(data=n.labs, aes(label=lab), x=3,y=0.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap("Host Phylum" + phylum~., nest_line = TRUE, ncol = 1)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")

p.phylum+
  ggsave2(here("vogdb","figures","sigma_HostPhylum.png"),
          width = 4,height = 8)


# make contingency table 
contin.t <- 
  d.sp%>%
  filter(! phylum %in% phyla_rm) %>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  select(phylum, n.sigma)%>%
  table()

m1 <- chisq.test(contin.t, simulate.p.value = TRUE, B = 1e6)
m2 <- fisher.test(contin.t, simulate.p.value = TRUE, B = 1e6)

##########################
# Plot by order
##########################

# summary table for N by order
n.labs <- d.sp%>%
  filter(! phylum %in% phyla_rm) %>%
  filter (!is.na(order)) %>% 
  group_by(order)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))%>%
  mutate(order=str_replace(order,"/","\n"))%>%
  mutate(order=fct_reorder(order,n))



d.sp %>%
  filter(! phylum %in% phyla_rm) %>%
  filter (!is.na(order)) %>% 
  ggplot( aes(n.sigma, group = order)) + 
  geom_bar(aes(y = ..prop..), stat="count") + 
  geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~order)+
  theme_classic()+
  panel_border()+
  ggtitle("phage sigma factor content by host order")+
  ggsave2(here("vogdb","figures","sigma_HostOrder.png"),
          width = 10,height = 6)

##########################
# Plot by genus within firmicutes
##########################

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
  geom_col(fill = "red4") + 
  geom_text(data=n.labs, aes(label=lab), x=2.75,y=.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested_wrap(~"Host genus (Firmicutes)" + genus.plot, scales = "fixed", nrow = 5)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")

p.genus+
  ggtitle("phage sigma factor content by host genus",
          "Firmicute infecting phages")+
  ggsave2(here("vogdb","figures","sigma_HostGenus_firmicutes.png"),
          width = 7,height = 7)



##########################
# Plot by viral type
##########################

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
  # tmp <- (!str_detect(d.faa$`virus lineage`,"Caudovirales"))%>%which()
  # d.faa$`virus name`[tmp]
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
  


p.vir+
  ggtitle("phage sigma factor content by viral family")+
  ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
          width = 10,height = 6)

# Who are the siphoviruses with 3 sigma factors?
d.sp%>%
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

#combine plots


right_col <- 
  plot_grid(p.vir, NULL,
            rel_heights = c(4.5, 1),
            ncol = 1)

plot_grid(p.phylum,p.genus,right_col,
          rel_widths = c(1.5, 3, 1.3),
          nrow = 1, labels = LETTERS)
  ggsave2(here("vogdb","figures","sigma_taxonomy.png"),
          width = 12,height = 10)


  
###############
# Summary plot
###############
  
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
  
  
  
  nsig_all %>% 
    ggplot(aes(n.sigma, perc))+
    geom_col(color = "grey30", fill = "white")+
    geom_col(data = nsig_firmi,fill = "grey", color ="black",
             alpha=.5, width = 0.6, position = position_nudge(x=0.15))+
    geom_col(data = nsig_bacil, fill = "black", color ="black",
             alpha=.8, width = 0.3, position = position_nudge(x=0.3))+
    geom_label(label = paste0("All phages: n=", sum(nsig_all$n.genomes)),
               x=2,y=.9, color = "black", fill = "white", hjust = 0)+
    geom_label(label = paste0("Firmicute phages: n=", sum(nsig_firmi$n.genomes)), 
               x=2,y=.85, fill = "grey80", color = "black",hjust = 0)+
    geom_label(label = paste0("Bacillus phages: n=", sum(nsig_bacil$n.genomes)), 
               x=2,y=.8, fill = "grey20", color = "white",hjust = 0)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    theme_classic()+
    panel_border(color = "black")
  # ggsave(here("vogdb","figures","nSigma_allVbacil.png"),
  #        width = 6,height = 6)
  
  dplot <- 
    bind_rows(
      nsig_all %>% mutate(host_tax = "Bacteria"),
      nsig_firmi %>% mutate(host_tax = "Firmicutes"),
      nsig_bacil %>% mutate(host_tax = "Bacillus")
     ) %>% 
    mutate(host_tax = fct_relevel(host_tax, c("Bacteria", "Firmicutes", "Bacillus")))

  p.tax <- dplot %>% 
    ggplot(aes(n.sigma, perc, fill = host_tax))+
    geom_col(color = "black", position = position_dodge2(padding = 0.2), width = 0.8)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    theme_classic()+
    panel_border(color = "black")+
    # scale_fill_grey(start = 0.8,end = 0.2, name = "Host taxa")+
    scale_fill_viridis_d( name = "Host taxa", direction = -1)+
    theme(legend.position = c(0.75,0.8),
          legend.key.size = unit(0.15, 'in'))
  ggsave(filename = here("vogdb","figures","nSigma_sum.png"),
        plot = p.tax, width = 3,height = 3)

  p.vir2 <- d.vir %>% 
    mutate(viral.family = fct_relevel(viral.family, c("Podoviridae", "Myoviridae",
                                                  "Siphoviridae", "Herelleviridae"))) %>% 
    ungroup() %>% 
    select(viral.family, n.sigma, perc) %>% 
    complete(n.sigma, viral.family, fill = list(perc=0)) %>% 
    ggplot(aes(n.sigma, perc, fill = viral.family))+
    geom_col(color = "black", position = position_dodge2(padding = 0.2), width = 0.8)+
    scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
    ylab("Phage genomes") +
    xlab("Sigma factors per genome")+
    theme_classic()+
    panel_border(color = "black")+
    # scale_fill_grey(start = 0.8,end = 0.2, name = "Host taxa")+
    scale_fill_viridis_d( name = "viral family", direction = -1)+
    theme(legend.position = c(0.7,0.8),
          legend.key.size = unit(0.1, 'in'))
  
  ggsave(filename = here("vogdb","figures","nSigma_virFam.png"),
         plot = p.vir2, width = 3,height = 3)
  
  plot_grid(p.tax,p.vir2, labels = letters) %>% 
save_plot(filename = here("vogdb","figures","nSigmaAB.png"),
          base_width = 6 , base_height = 3)
  