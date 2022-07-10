library(tidyverse)
library(here)


#load data frame of phage sigma factors generated in (B)
load(file=here("vogdb","data","vog_sigma_clean.RData"))

d.faa <- d.faa%>%
  mutate(taxon=as.numeric(taxon))


# import all viruses used to assemble VOGs
# vog.sp <-  read_tsv(here("vogdb/data","vogdb_downloads","vog.species.list"), trim_ws = T) %>% 
#   as_tibble( .name_repair = "universal") %>% 
#   # keep only phages
#   filter(phage.nonphage == "phage")
vog.sp <- read_csv(here("vogdb/data/vog-members_cd-clstr.csv"))

#import virus-host data
# download.file(url = "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv",
#               destfile = here("vogdb/data", "virushostdb.tsv"))
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

# remove phages with no host at phylum level (there are 3 such rows)
d.sp <-
  d.sp %>% filter(!is.na(phylum))

# extract viral family
d.sp <- 
  d.sp %>% 
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

# there are three phages with viral family "unclassified"
d.sp %>% filter(str_detect(viral.family,"unclassified"))
# Nitrincola phage 1M3-16 - excluding, no information beside genome sequence
# Streptococcus phage 315.5 - excluding ("putative prophage ...There is no evidence that it can be propagated as a functional phage.")
# Escherichia phage Henu7 - NCBI taxonomy now ahs this :    Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Caudovirales; Drexlerviridae; Tempevirinae; Henuseptimavirus; Escherichia virus Henu7
  # do assigning family Drexlerviridae
d.sp <- d.sp %>% 
  mutate(viral.family = if_else(str_detect(.species.name,"Henu7"),
                                "Drexlerviridae",viral.family %>% as.character())) %>% 
  # filter(str_detect(.species.name,"Henu7")) %>% view()
  filter(!str_detect(viral.family,"unclassified"))

# add host to sigma factor data

d.faa <- left_join(d.faa,d.sp,by=c("taxon"="tax.id"))

# filter out sequences from phages that were excluded above
# (clustering analysis and host matching)
# these are mostly crAss phages (no host)
d.faa <-
  d.faa %>% 
  filter(!is.na(phylum)) %>% 
  filter(!is.na(n_clstr))

# for the rest of the analysis we will also remove vOTU clustering duplicates
clstr_dups <- 
  d.faa %>% filter (!is_rep)
# this removes 14 sigma factor genes
d.faa <- 
  d.faa %>% filter (is_rep)
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
#582

# Phylum independent summary of sigma factors/genome
d.sp%>%
  group_by(tax.id)%>%
  group_by(n.sigma)%>%
  summarise(n.genomes=n())%>%
  mutate(perc=100*n.genomes/sum(n.genomes))

# n.sigma n.genomes   perc
# 1       0      3768 86.6  
# 2       1       507 11.7  
# 3       2        46  1.06 
# 4       3        29  0.667


# summarise by phylum
d.sp%>%
  group_by(phylum)%>%
  mutate(has_sigma = n.sigma >0,
         multi_sigma = n.sigma >1) %>%
  summarise(n=n(),
            perc_sigma=100*sum(has_sigma)/n(),
            perc_multi=100*sum(multi_sigma)/n())

# phylum                                  n perc_sigma perc_multi
# 1 Actinobacteria                       1314      0.609      0    
# 2 Aquificae                               1      0          0    
# 3 Bacteroidetes/Chlorobi group           83      2.41       0    
# 4 Cyanobacteria/Melainabacteria group   100     63          2    
# 5 Deinococcus-Thermus                     5      0          0    
# 6 Firmicutes                            864     20.7        7.06 
# 7 Fusobacteria                            1      0          0    
# 8 Proteobacteria                       1970     16.8        0.609
# 9 PVC group                               6      0          0    
# 10 Spirochaetes                            3      0          0    
# 11 Tenericutes                             3      0          0 

